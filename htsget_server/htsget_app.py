import sys
import json
import yaml
import requests
from tempfile import NamedTemporaryFile

from tornado import ioloop, web
from pyCGA.opencgarestclients import OpenCGAClient
from pysam import AlignmentFile


def _get_coordinates(args):
    # Start: 0-based, inclusive. End: 0-based exclusive.
    ref = str(args['referenceName']) if 'referenceName' in args else None
    start = int(args['start']) + 1 if 'start' in args else None
    end = int(args['end']) if 'end' in args else None
    return ref, start, end


class ConfigHandler:
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

        host = config_dict['rest']['hosts'][0]
        if not (host.startswith('http://') or host.startswith('https://')):
            self.host = 'http://' + host
        self.version = config_dict['version']

    def get_configuration(self):
        return {'version': self.version, 'rest': {'hosts': [self.host]}}


class BaseHandler(web.RequestHandler):
    def __init__(self, application, request, **kwargs):
        super(BaseHandler, self).__init__(application, request)
        self.oc = None
        self.config = kwargs['config']

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
        try:
            if username and password:
                self.oc = OpenCGAClient(configuration=self.config,
                                        user=username,
                                        pwd=password,
                                        auto_refresh=False)
            elif username and auth_token:
                session_id = self._get_session_id(auth_token)
                self.oc = OpenCGAClient(configuration=self.config,
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

    def _get_file_info(self, study, file_id=None, sample=None, bioformat=None,
                       name=None, **kwargs):
        if file_id:
            file_info = self.oc.files.search(
                study=study, id=file_id, include='id,uri', **kwargs
            ).get()
        else:
            file_info = self.oc.files.search(
                study=study, samples=sample, bioformat=bioformat,
                name=name, include='id,uri', **kwargs
            ).get()

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

    def _download_file(self, file_id, study, auth_token, prefix='htsget',
                      suffix='', _dir='/tmp', mode='wb', delete=True):
        r = requests.get(
            '{h}/webservices/rest/{v}/utils/ranges/{f}?study={s}&sid={sid}'.format(
                h=self.config['rest']['hosts'][0], v=self.config['version'],
                f=file_id, s=study, sid=self._get_session_id(auth_token)
            )
        )
        ntf = NamedTemporaryFile(prefix=prefix, suffix=suffix, dir=_dir,
                                 mode=mode, delete=delete)
        ntf.write(r.content)
        ntf.close()
        return ntf.name

    def _get_args(self):
        args = {arg: self.request.arguments[arg][0]
                for arg in self.request.arguments}

        ftype = args.get('format', None)
        ref = args.get('referenceName', None)
        start = args.get('start', None)
        end = args.get('end', None)

        # UnsupportedFormat error if the requested format is not supported
        class_name = self.__class__.__name__
        if not ftype:
            args['format'] = 'BAM' if class_name == 'ReadsHandler' else 'VCF'
        else:
            if (class_name == 'ReadsHandler' and ftype not in ['BAM']) or \
                    (class_name == 'VariantsHandler' and ftype not in ['VCF']):
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
        self.auth_token = self.request.headers.get('Authorization', '')
        self.username = self.request.headers.get('username', '')

    def _check_ref(self, study, bam_id, reference):
        references = [
            ref['sequenceName']
            for ref in self.oc.files.search(
                study=study, id=bam_id
            ).get()[0]['attributes']['alignmentHeader']['sequenceDiccionary']
        ]
        if reference not in references:
            raise web.HTTPError(
                reason='NotFound::Requested reference does not exist',
                status_code=404
            )

    def get(self, study, sample):
        # Login into OpenCGA
        self._opencga_login(username='bertha', auth_token=self.auth_token)

        # Getting arguments
        args = self._get_args()

        # Getting BAM ID
        bam_id = self._get_file_info(
            study=study, sample=sample, bioformat='ALIGNMENT', name='~.bam$',
            **{'attributes.gelStatus': 'READY'}
        )['id']

        # Getting coordinates
        ref, start, end = _get_coordinates(args)
        self._check_ref(study, bam_id, ref)

        # Creating the response URL to retrieve data
        url = '{host}/data?format={format}&study={study_id}&fileId={file_id}'
        url += '&referenceName={reference_name}' if ref else ''
        url += '&start={start}' if start else ''
        url += '&end={end}' if end else ''
        url = url.format(host='http://' + self.request.host,
                         format=args['format'],
                         study_id=study,
                         file_id=bam_id,
                         reference_name=ref,
                         start=start,
                         end=end)

        # Generating URLs
        response = {'htsget': {'format': 'BAM', 'urls': []}}
        bam_chunks = [{'url': url,
                       'headers': {'Authorization': self.auth_token,
                                   'username': self.username}
                       }]
        response['htsget']['urls'] = response['htsget']['urls'] + bam_chunks

        self.write(response)


class DataHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(DataHandler, self).__init__(application, request, **kwargs)
        self.auth_token = self.request.headers.get('Authorization', '')
        self.username = self.request.headers.get('username', '')

    def get(self):
        # Login into OpenCGA
        self._opencga_login(username='bertha', auth_token=self.auth_token)

        # Getting arguments
        args = self._get_args()

        # Getting BAM ID
        file_uri = self._get_file_info(
            study=args['study'], file_id=args['fileId']
        )['uri'].replace('file://', '')

        ref, start, end = _get_coordinates(args)
        print ref, start, end

        # TODO REMOVE DEBUG LINES:
        print file_uri
        file_uri = '/home/dgil/dev/gel-htsget/htsget_server/LP3001241-DNA_D06.mini.bam'

        samfile = AlignmentFile(file_uri, 'rb')
        response = {'result': '\n'.join(
            [str(samfile.header).rstrip('\n')] +
            [read.to_string() for read in samfile.fetch(reference=ref,
                                                        start=start,
                                                        end=end)]
        )}
        self.write(response)


class Htsget:
    def __init__(self, config):

        config = ConfigHandler(config).get_configuration()

        self.application = web.Application([
            (r'/user/login', LoginHandler, dict(config=config)),
            (r'/user/refresh-token', RefreshTokenHandler, dict(config=config)),
            (r'/reads/(.*)/(.*)', ReadsHandler, dict(config=config)),
            (r'/data', DataHandler, dict(config=config))
        ])
        self.application.listen(8888)

    @staticmethod
    def start():
        ioloop.IOLoop.instance().start()

    @staticmethod
    def stop():
        ioloop.IOLoop.instance().stop()


def main():
    htsget = Htsget(config='./opencga.json')
    htsget.start()


if __name__ == '__main__':
    sys.exit(main())
