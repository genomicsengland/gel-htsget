import os

from setuptools import setup

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

dir_path = os.path.dirname(os.path.realpath(__file__))

with open(os.path.join(dir_path, './VERSION'), 'r') as version_file:
    version = str(version_file.readline()).strip()

reqs = [
    'pyCGA==1.3.0',
    'tornado==5.1.1',
    'pyyaml>=4.2b1',
    'pysam==0.15.2'
]


setup(
    name='htsget',
    version=version,
    packages=['htsget_server'],
    url='',
    license='',
    include_package_data=True,
    author='antonior,dapregi',
    author_email='antonio.rueda-martin@genomicsengland.co.uk,daniel.perez-gil@genomicsengland.co.uk',
    description='Genomics England htsget',
    install_requires=reqs
)
