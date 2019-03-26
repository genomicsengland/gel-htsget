import os
import sys
import json
import yaml
from tempfile import NamedTemporaryFile
from tornado import ioloop, web
from pyCGA.opencgarestclients import OpenCGAClient
from pysam import AlignmentFile, VariantFile

_READS_FORMAT = ['BAM']
_VARIANTS_FORMAT = ['VCF']


class OpenCGAConfigHandler:
    def __init__(self, config_fpath):
        self.config_fpath = config_fpath
        self.host = None
        self.version = None
        self._set_up_config_params()

    def _set_up_config_params(self):
        try:
            config_fhand = open(self.config_fpath, 'r')
        except Exception:
            msg = 'Unable to read config file "' + self.config_fpath + '"'
            raise IOError(msg)

        if self.config_fpath.endswith('.yml') or self.config_fpath.endswith('.yaml'):
            config_dict = yaml.safe_load(config_fhand)
        elif self.config_fpath.endswith('.json'):
            config_dict = json.loads(config_fhand.read())
        else:
            msg = 'File path must end in "yml", "yaml" or "json"'
            raise IOError(msg)

        config_fhand.close()

        host = config_dict['rest']['hosts'][0]
        if not (host.startswith('http://') or host.startswith('https://')):
            self.host = 'http://' + host
        self.version = config_dict['version']

    def get_configuration(self):
        return {'version': self.version, 'rest': {'hosts': [self.host]}}


class VCFTypeConfig:
    def __init__(self, config_fpath):
        self.config_fpath = config_fpath
        self.vcf_types = {}
        self._get_vcf_types()

    def _get_vcf_types(self):
        try:
            config_fhand = open(self.config_fpath, 'r')
        except Exception:
            msg = 'Unable to read config file "' + self.config_fpath + '"'
            raise IOError(msg)

        for line in config_fhand:
            study, type_, regex = line.rstrip('\n').split()
            self.vcf_types.setdefault(study, {}).update({type_: regex})

        config_fhand.close()

    def get_vcf_types(self):
        return self.vcf_types


class BaseHandler(web.RequestHandler):
    def __init__(self, application, request, **kwargs):
        super(BaseHandler, self).__init__(application, request)
        self.oc = None
        self.opencga_config = kwargs['opencga_config']

    def options(self):
        method = self.request.headers.get('Access-Control-Request-Method', '')
        if method and method == 'GET':
            self.set_header('Access-Control-Max-Age', 2592000)
            headers = self.request.headers.get('Access-Control-Request-Headers', '')
            if headers:
                self.set_header('Access-Control-Allow-Headers', headers)
        self.set_status(204)
        self.finish()

    def set_default_headers(self):
        origin = self.request.headers.get('Origin', '')
        if origin:
            self.set_header('Access-Control-Allow-Origin', origin)
        self.set_header('Content-Type', 'application/vnd.ga4gh.htsget.v1.1.1+json; charset=utf-8')
        self.set_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')

    def write_error(self, status_code, **kwargs):
        # self.set_header('Content-Type', 'application/json')
        try:
            error, msg = self._reason.split('::')
        except Exception, e:
            error, msg = self._reason, e.message

        self.finish(json.dumps({
            'htsget': {
                'error': error,
                'message': msg,
            }
        }))

    @staticmethod
    def _get_session_id(auth_token):
        try:
            bearer, sid = auth_token.strip().split(' ')
            assert bearer == 'Bearer'
        except (ValueError, AssertionError):
            raise ValueError('Invalid Authorization token')
        return sid

    def _opencga_login(self, username=None, password=None, auth_token=None):
        if not (username and (password or auth_token)):
            msg = 'Username/password or username/token not provided'
            raise web.HTTPError(reason='InvalidAuthentication::' + msg,
                                status_code=401)

        try:
            if username and password:
                self.oc = OpenCGAClient(configuration=self.opencga_config,
                                        user=username,
                                        pwd=password,
                                        auto_refresh=False)
            elif username and auth_token:
                session_id = self._get_session_id(auth_token)
                self.oc = OpenCGAClient(configuration=self.opencga_config,
                                        user=username,
                                        session_id=session_id,
                                        auto_refresh=False)
        except Exception, e:
            try:
                msg = json.loads(e.message)['error']
            except (ValueError, KeyError):
                msg = str(e)
            raise web.HTTPError(reason='InvalidAuthentication::' + msg,
                                status_code=401)

    def _get_file_info(self, study, file_id=None, sample=None, path=None,
                       bioformat=None, name=None, **kwargs):
        include = 'id,uri'
        if file_id and not (name or path):
            file_info = self.oc.files.search(
                study=study, id=file_id, include=include, **kwargs
            ).get()
        elif sample and name and not (file_id or path):
            file_info = self.oc.files.search(
                study=study, samples=sample, bioformat=bioformat,
                name=name, include=include, **kwargs
            ).get()
        elif path and not (file_id or sample):
            file_info = self.oc.files.search(
                study=study, bioformat=bioformat, path=path, include=include,
                **kwargs
            ).get()
        else:
            file_info = None

        if not file_info:
            msg = 'NotFound::No file found for sample "{}" in study "{}"'
            raise web.HTTPError(reason=msg.format(sample, study),
                                status_code=404)
        elif len(file_info) > 1:
            msg = 'NotFound::Multiple files found for sample "{}" in study "{}": "{}"'
            raise web.HTTPError(reason=msg.format(sample, study, file_info),
                                status_code=404)
        else:
            return file_info[0]

    def _get_args(self, endpoint):

        assert endpoint in ['reads', 'variants', 'data']

        args = {arg: self.request.arguments[arg][0]
                for arg in self.request.arguments}

        ftype = args.get('format', None)
        ref = args.get('referenceName', None)
        start = args.get('start', None)
        end = args.get('end', None)

        # UnsupportedFormat error if the requested format is not supported
        if not ftype:
            args['format'] = 'BAM' if endpoint == 'reads' else 'VCF'
        else:
            if (endpoint == 'reads' and ftype not in _READS_FORMAT) or \
                    (endpoint == 'variants' and ftype not in _VARIANTS_FORMAT):
                msg = 'UnsupportedFormat::Requested format is not supported'
                raise web.HTTPError(reason=msg, status_code=400)

        # InvalidInput error if start but no reference or referenceName is "*"
        if (start or end) and (ref == '*' or not ref):
            msg = ('InvalidInput::Coordinates are specified and either no'
                   ' reference is specified or referenceName is "*"')
            raise web.HTTPError(reason=msg, status_code=400)

        # InvalidRange error if start is greater than end
        if start and end and int(start) > int(end):
            msg = 'InvalidRange::Start greater than end'
            raise web.HTTPError(reason=msg, status_code=400)

        return args

    def _check_ref(self, endpoint, study, file_id, reference):

        assert endpoint in ['reads', 'variants', 'data']

        if endpoint == 'reads':
            references = [
                ref['sequenceName']
                for ref in self.oc.files.search(
                    study=study, id=file_id
                ).get()[0]['attributes']['alignmentHeader']['sequenceDiccionary']
            ]
        elif endpoint == 'variants':
            references = [
                i['id'] for i in self.oc.files.search(
                    study=study, id=file_id
                ).get()[0]['attributes']['variantFileMetadata']['header'][
                    'complexLines']
                if i['key'] == 'contig'
            ]
        else:
            references = []

        if reference not in references:
            raise web.HTTPError(
                reason='NotFound::Requested reference does not exist',
                status_code=404
            )

    @staticmethod
    def _get_coordinates(args):
        # Start: 0-based, inclusive. End: 0-based exclusive.
        ref = str(args['referenceName']) if 'referenceName' in args else None
        start = int(args['start']) + 1 if 'start' in args else None
        end = int(args['end']) if 'end' in args else None
        return ref, start, end

    @staticmethod
    def _create_response(ref, start, end, study, host, format_, file_id, token,
                         username):
        # Creating the response URL to retrieve data
        url = '{host}/data?format={format}&study={study_id}&fileId={file_id}'
        url += '&referenceName={reference_name}' if ref else ''
        url += '&start={start}' if start else ''
        url += '&end={end}' if end else ''
        url = url.format(host='http://' + host,
                         format=format_,
                         study_id=study,
                         file_id=file_id,
                         reference_name=ref,
                         start=start,
                         end=end)

        response = {'htsget': {'format': format_, 'urls': []}}
        chunks = [{'url': url,
                       'headers': {'Authorization': token,
                                   'username': username}
                       }]
        response['htsget']['urls'] += chunks

        return response


class LoginHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(LoginHandler, self).__init__(application, request, **kwargs)

    def post(self):
        username = self.get_argument('username', '')
        password = self.get_argument('password', '')

        if not username:
            raise web.HTTPError(reason='InvalidInput::No username provided',
                                status_code=400)
        elif not password:
            raise web.HTTPError(reason='InvalidInput::No password provided',
                                status_code=400)
        else:
            self._opencga_login(username=username, password=password)
            self.write({'sessionId': self.oc.session_id})


class RefreshTokenHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(RefreshTokenHandler, self).__init__(application, request, **kwargs)

    def post(self):
        username = self.get_argument('username', '')
        token = self.get_argument('token', '')

        if not username:
            raise web.HTTPError(reason='InvalidInput::No username provided',
                                status_code=400)
        elif not token:
            raise web.HTTPError(reason='InvalidInput::No token provided',
                                status_code=400)
        else:
            self._opencga_login(username=username, auth_token=token)
            self.write({'sessionId': self.oc.session_id})


class ReadsHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(ReadsHandler, self).__init__(application, request, **kwargs)
        self.endpoint = 'reads'
        self.auth_token = self.request.headers.get('Authorization', '')
        self.username = self.request.headers.get('username', '')

    def get(self, study, sample):
        # Login into OpenCGA
        self._opencga_login(username=self.username, auth_token=self.auth_token)

        # Getting arguments
        args = self._get_args(self.endpoint)

        # Getting BAM ID
        bam_id = self._get_file_info(
            study=study, sample=sample, bioformat='ALIGNMENT', name='~.bam$',
            **{'attributes.gelStatus': 'READY'}
        )['id']

        # Getting coordinates
        ref, start, end = self._get_coordinates(args)
        self._check_ref(self.endpoint, study, bam_id, ref)

        # Creating the response URL to retrieve data
        response = self._create_response(
            ref, start, end, study, self.request.host, args['format'], bam_id,
            self.auth_token, self.username
        )
        self.write(response)


class ReadsByPathHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(ReadsByPathHandler, self).__init__(application, request, **kwargs)
        self.endpoint = 'reads'
        self.auth_token = self.request.headers.get('Authorization', '')
        self.username = self.request.headers.get('username', '')

    def get(self, study, path):
        # Login into OpenCGA
        self._opencga_login(username=self.username, auth_token=self.auth_token)

        # Getting arguments
        args = self._get_args(self.endpoint)

        # Getting BAM ID
        bam_id = self._get_file_info(
            study=study, path=path.replace(':', '/'), bioformat='ALIGNMENT',
            **{'attributes.gelStatus': 'READY'}
        )['id']

        # Getting coordinates
        ref, start, end = self._get_coordinates(args)
        self._check_ref(self.endpoint, study, bam_id, ref)

        # Creating the response URL to retrieve data
        response = self._create_response(
            ref, start, end, study, self.request.host, args['format'], bam_id,
            self.auth_token, self.username
        )
        self.write(response)


class VariantsHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(VariantsHandler, self).__init__(application, request, **kwargs)
        self.endpoint = 'variants'
        self.auth_token = self.request.headers.get('Authorization', '')
        self.username = self.request.headers.get('username', '')
        self.extra_options = kwargs

    def get(self, study, vcf_type, sample):
        # Login into OpenCGA
        self._opencga_login(username=self.username, auth_token=self.auth_token)

        # Getting arguments
        args = self._get_args(self.endpoint)

        # Getting VCF ID
        if 'vcf_types' not in self.extra_options:
            msg = 'InvalidInput::Configuration file with VCF types required'
            raise web.HTTPError(reason=msg, status_code=400)

        if study not in self.extra_options['vcf_types']:
            msg = 'InvalidInput::Study not defined in the VCF types configuration file'
            raise web.HTTPError(reason=msg, status_code=400)

        if vcf_type not in self.extra_options['vcf_types'][study]:
            msg = 'InvalidInput::VCF type not defined for the study "{}" in the VCF types configuration file'
            raise web.HTTPError(reason=msg.format(study), status_code=400)

        regex = self.extra_options['vcf_types'][study][vcf_type]
        vcf_id = self._get_file_info(
            study=study, sample=sample, bioformat='VARIANT', name=regex,
            **{'attributes.gelStatus': 'READY'}
        )['id']

        # Getting coordinates
        ref, start, end = self._get_coordinates(args)
        self._check_ref(self.endpoint, study, vcf_id, ref)

        # Creating the response URL to retrieve data
        response = self._create_response(
            ref, start, end, study, self.request.host, args['format'], vcf_id,
            self.auth_token, self.username
        )
        self.write(response)


class VariantsByPathHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(VariantsByPathHandler, self).__init__(application, request,
                                                    **kwargs)
        self.endpoint = 'variants'
        self.auth_token = self.request.headers.get('Authorization', '')
        self.username = self.request.headers.get('username', '')

    def get(self, study, path):
        # Login into OpenCGA
        self._opencga_login(username=self.username, auth_token=self.auth_token)

        # Getting arguments
        args = self._get_args(self.endpoint)

        # Getting VCF ID
        vcf_id = self._get_file_info(
            study=study, path=path.replace(':', '/'), bioformat='VARIANT',
            **{'attributes.gelStatus': 'READY'}
        )['id']

        # Getting coordinates
        ref, start, end = self._get_coordinates(args)
        self._check_ref(self.endpoint, study, vcf_id, ref)

        # Creating the response URL to retrieve data
        response = self._create_response(
            ref, start, end, study, self.request.host, args['format'], vcf_id,
            self.auth_token, self.username
        )
        self.write(response)


class DataHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(DataHandler, self).__init__(application, request, **kwargs)
        self.endpoint = 'data'
        self.auth_token = self.request.headers.get('Authorization', '')
        self.username = self.request.headers.get('username', '')

        self.ntf = NamedTemporaryFile(prefix='htsget', suffix='', dir='/tmp',
                                      mode='wb', delete=False)
        self._buf_size = 1000000000

    def get(self):
        # Login into OpenCGA
        self._opencga_login(username=self.username, auth_token=self.auth_token)

        # Getting arguments
        args = self._get_args(self.endpoint)

        # Getting file uri
        file_uri = self._get_file_info(
            study=args['study'], file_id=args['fileId']
        )['uri'].replace('file://', '')

        ref, start, end = self._get_coordinates(args)

        # Getting file data
        if args['format'] in _READS_FORMAT:
            self.set_header('Content-Type', 'application/vnd.ga4gh.bam')
            samfile = AlignmentFile(file_uri, 'rb')
            bamfile = AlignmentFile(self.ntf.name, header=samfile.header, mode='wb')
            for read in samfile.fetch(ref, start, end):
                bamfile.write(read)
            samfile.close()
            bamfile.close()
        elif args['format'] in _VARIANTS_FORMAT:
            self.set_header('Content-Type', 'application/vnd.ga4gh.vcf')
            vcffile_in = VariantFile(file_uri, 'r')
            vcffile_out = VariantFile(self.ntf.name, header=vcffile_in.header, mode='w')
            for read in vcffile_in.fetch(ref, start, end):
                vcffile_out.write(read)
            vcffile_in.close()
            vcffile_out.close()

        # Data blocks should not exceed ~1GB
        with open(self.ntf.name, 'rb') as f:
            while True:
                data = f.read(self._buf_size)
                if not data:
                    break
                self.write(data)

        os.remove(self.ntf.name)


class Htsget:
    def __init__(self, opencga_config, vcf_types_config=None):

        opencga_config = OpenCGAConfigHandler(opencga_config).get_configuration()

        vcf_types = None
        if vcf_types_config is not None:
            vcf_types = VCFTypeConfig(vcf_types_config).get_vcf_types()

        self.application = web.Application([
            (r'/user/login', LoginHandler,
             dict(opencga_config=opencga_config)),
            (r'/user/refresh-token', RefreshTokenHandler,
             dict(opencga_config=opencga_config)),
            (r'/reads/(?!bypath)(.*)/(.*)', ReadsHandler,
             dict(opencga_config=opencga_config)),
            (r'/reads/bypath/(.*)/(.*)', ReadsByPathHandler,
             dict(opencga_config=opencga_config)),
            (r'/variants/(?!bypath)(.*)/(.*)/(.*)', VariantsHandler,
             dict(opencga_config=opencga_config, vcf_types=vcf_types) if vcf_types
             else dict(opencga_config=opencga_config)),
            (r'/variants/bypath/(.*)/(.*)', VariantsByPathHandler,
             dict(opencga_config=opencga_config)),
            # TODO data endpoint should be https
            (r'/data', DataHandler,
             dict(opencga_config=opencga_config))
        ])
        self.application.listen(8888)

    @staticmethod
    def start():
        ioloop.IOLoop.instance().start()

    @staticmethod
    def stop():
        ioloop.IOLoop.instance().stop()


def main():
    htsget = Htsget(opencga_config='./opencga.json',
                    vcf_types_config='./vcf_types.tsv')
    htsget.start()


if __name__ == '__main__':
    sys.exit(main())
