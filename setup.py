from setuptools import setup

reqs = [
    'pyCGA==1.3.1',
    'tornado==5.1.1',
    'pyyaml>=4.2b1',
    'pysam==0.15.2'
]


setup(
    name='htsget',
    version='0.1.0',
    packages=['htsget_server'],
    url='',
    license='',
    include_package_data=True,
    author='antonior,dapregi',
    author_email='antonio.rueda-martin@genomicsengland.co.uk,daniel.perez-gil@genomicsengland.co.uk',
    description='Genomics England htsget',
    install_requires=reqs
)
