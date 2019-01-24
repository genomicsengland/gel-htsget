import sys
import json
import yaml

from tornado import ioloop, web

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

    def write_error(self, status_code, **kwargs):
        self.set_header('Content-Type', 'application/json')
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
    def get_session_id(auth_token):
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
                                        pwd=password)
            elif auth_token:
                session_id = self.get_session_id(auth_token)
                self.oc = OpenCGAClient(configuration=self.config,
                                        session_id=session_id)
        except Exception, e:
            try:
                msg = json.loads(e.message)['error']
            except (ValueError, KeyError):
                msg = str(e)
            raise web.HTTPError(reason='InvalidAuthentication::' + msg,
                                status_code=401)

    def get_file_id(self, study, sample, bioformat, name, **kwargs):
        file_info = self.oc.files.search(
            study=study, samples=sample, bioformat=bioformat, name=name,
            include='id', **kwargs
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
            return file_info[0]['id']


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


class ReadsHandler(BaseHandler):
    def __init__(self, application, request, **kwargs):
        super(ReadsHandler, self).__init__(application, request, **kwargs)
        self.auth_token = self.request.headers.get('Authorization', '')

    def get(self, study, sample):
        self.opencga_login(auth_token=self.auth_token)

        file_id = self.get_file_id(
            study=study, sample=sample, bioformat='ALIGNMENT', name='~.bam$',
            **{'attributes.gelStatus': 'READY'}
        )

        opencga_url = '{host}/webservices/rest/{version}/utils/ranges/{file_id}?study={study_id}'
        url = opencga_url.format(host=self.config['rest']['hosts'][0],
                                 version=self.config['version'],
                                 file_id=file_id,
                                 study_id=study)
        print url

        auth = self.auth_token
        bam_header_url = 'XXX'
        bam_body_url = 'YYY'
        byte_range = 'Range ' + 'DDDD'

        response = {
            'ID': sample,
            'OPTIONS': self.request.arguments,
            "htsget": {
                "format": "BAM",
                "urls": [
                    {
                        "url": bam_header_url,
                        "headers": {
                            "Authorization": auth,
                            "Range": byte_range
                        }
                    },
                    {
                        "url": bam_body_url,
                        "headers": {
                            "Authorization": auth,
                            "Range": byte_range
                        }
                    }
                ]
            }
        }
        self.write(response)


class Htsget:
    def __init__(self, config):

        config = ConfigHandler(config).get_configuration()

        self.application = web.Application([
            (r"/user/login", LoginHandler, dict(config=config)),
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
    htsget = Htsget(config='./htsget/opencga.json')
    htsget.start()


if __name__ == "__main__":
    sys.exit(main())
