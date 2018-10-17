from setuptools import setup, find_packages
from codecs import open
from os import path

__version__ = '0.0.1'

here = path.abspath(path.dirname(__file__))


DESCRIPTION = 'Collection of object-oriented entities related by the central dogma of biology.',
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with open(path.join(here, 'requirements.txt'), encoding='utf-8') as f:
    all_reqs = f.read().split('\n')

install_requires = [x.strip() for x in all_reqs if 'git+' not in x]
dependency_links = [x.strip().replace('git+', '')
                    for x in all_reqs if x.startswith('git+')]

setup(
    name='dogma',
    version=__version__,
    description=DESCRIPTION,
    long_description=long_description,
    license='MIT',
    author='Griffin Clausen',
    author_email='griffinclausen@gmail.com',
    keywords='',
    entry_points='''
        [console_scripts]
        genetic_code=dogma.genetic_code:test
    ''',
    packages=find_packages(exclude=['docs', 'tests*']),
    include_package_data=True,
    install_requires=install_requires,
    dependency_links=dependency_links
)