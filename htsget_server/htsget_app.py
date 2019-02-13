import os
import sys
import json
import yaml
import requests
from tempfile import NamedTemporaryFile
from tornado import ioloop, web

from index import Index
from pyCGA.opencgarestclients import OpenCGAClient


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
        return {"version": self.version, "rest": {"hosts": [self.host]}}


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

    def opencga_login(self, username=None, password=None, auth_token=None):
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

    def get_file_info(self, study, sample, bioformat, name, **kwargs):
        file_info = self.oc.files.search(
            study=study, samples=sample, bioformat=bioformat, name=name,
            include='id,uri', **kwargs
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

    def download_file(self, file_id, study, auth_token, prefix='htsget',
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

    def check_arguments(self):
        format = self.request.arguments.get('format', None)
        ref = self.request.arguments.get('referenceName', None)
        start = self.request.arguments.get('start', None)
        end = self.request.arguments.get('end', None)

        class_name = self.__class__.__name__
        if format and class_name == 'ReadsHandler' and format[0] not in ['BAM']:
            msg = 'InvalidInput::Requested format is not supported'
            raise web.HTTPError(reason=msg, status_code=400)
        if format and class_name == 'VariantsHandler' and format[0] not in ['VCF']:
            msg = 'InvalidInput::Requested format is not supported'
            raise web.HTTPError(reason=msg, status_code=400)
        if (start or end) and (ref == '*' or not ref):
            msg = 'InvalidInput::Coordinates are specified and either no reference is specified or referenceName is "*"'
            raise web.HTTPError(reason=msg, status_code=400)
        if start and end and int(start[0]) > int(end[0]):
            msg = 'InvalidRange::Start greater than end'
            raise web.HTTPError(reason=msg, status_code=400)


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
            self.opencga_login(username=username, password=password)
            self.write({"sessionId": self.oc.session_id})


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
            self.opencga_login(username=username, auth_token=token)
            self.write({"sessionId": self.oc.session_id})


class ReadsHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(ReadsHandler, self).__init__(application, request, **kwargs)
        self.auth_token = self.request.headers.get('Authorization', '')
        self.username = self.request.headers.get('username', '')

    @staticmethod
    def _get_header_length(bai_fpath):
        index = Index.from_file(filename=bai_fpath)
        return index.header_length()

    def _get_refs(self, study, bam_id):
        return [
            ref['sequenceName']
            for ref in self.oc.files.search(
                study=study, id=bam_id
            ).get()[0]['attributes']['alignmentHeader']['sequenceDiccionary']
        ]

    @staticmethod
    def _get_bytes(bai_fpath, refs, ref, start, end):
        index = Index.from_file(filename=bai_fpath)
        try:
            reference = refs.index(ref)
        except Exception:
            raise web.HTTPError(
                reason='NotFound::Requested reference does not exist',
                status_code=404
            )

        return [(i[0], i[0] + i[1])
                for i in index.region_offset_iter(ref=reference,
                                                  beg=start, end=end)]

    def _get_coordinates(self):
        ref = self.request.arguments.get('referenceName', None)
        start = self.request.arguments.get('start', None)
        end = self.request.arguments.get('end', None)
        return str(ref[0]), int(start[0]), int(end[0])

    def get(self, study, sample):
        self.check_arguments()

        self.opencga_login(username=self.username, auth_token=self.auth_token)

        response = {"htsget": {"format": "BAM", "urls": []}}

        # Getting BAM and BAI IDs
        bam_id = self.get_file_info(
            study=study, sample=sample, bioformat='ALIGNMENT', name='~.bam$',
            **{'attributes.gelStatus': 'READY'}
        )['id']
        bai_id = self.get_file_info(
            study=study, sample=sample, bioformat='ALIGNMENT', name='~.bai$',
            **{'attributes.gelStatus': 'READY'}
        )['id']

        # BAM chunks byte range
        refs = self._get_refs(study, bam_id)
        ref, start, end = self._get_coordinates()
        bai_fpath = self.download_file(
            bai_id, study, self.auth_token, prefix='htsget-{}'.format(sample),
            suffix='.bai', delete=False
        )
        bam_header_length = self._get_header_length(bai_fpath)
        byte_range_list = self._get_bytes(bai_fpath, refs, ref, start, end)
        os.unlink(bai_fpath)

        # OpenCGA URL
        opencga_url = '{host}/webservices/rest/{version}/utils/ranges/{file_id}?study={study_id}'
        bam_body_url = opencga_url.format(host=self.config['rest']['hosts'][0],
                                          version=self.config['version'],
                                          file_id=bam_id,
                                          study_id=study)

        # Generating URLs
        bam_chunks = [
            {
                "url": bam_body_url,
                "headers": {
                    "Authorization": self.auth_token,
                    "Range": 'bytes={}-{}'.format(byte_range[0], byte_range[1])
                }
            } for byte_range in [(0, bam_header_length+65536)] + byte_range_list
            # } for byte_range in [(0, bam_header_length+19558)] + byte_range_list
        ] + [
            {'url': 'data:application/vnd.ga4gh.bam;base64,H4sIBAAAAAAA/wYAQkMCABsAAwAAAAAAAAAAAA=='}  # EOF
        ]
        response['htsget']['urls'] = response['htsget']['urls'] + bam_chunks

        self.write(response)


class Htsget:
    def __init__(self, config):

        config = ConfigHandler(config).get_configuration()

        self.application = web.Application([
            (r"/user/login", LoginHandler, dict(config=config)),
            (r"/user/refresh-token", RefreshTokenHandler, dict(config=config)),
            (r"/reads/(.*)/(.*)", ReadsHandler, dict(config=config))
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


if __name__ == "__main__":
    sys.exit(main())
