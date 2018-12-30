from setuptools import setup, find_packages
import io
from os import path

__version__ = '1.0.1'

here = path.abspath(path.dirname(__file__))


DESCRIPTION = """
Collection of object-oriented entities related by the central dogma of biology.
""".strip()

with io.open('README.rst', 'rt', encoding='utf8') as f:
    readme = f.read()

with io.open(path.join(here, 'requirements.txt'), encoding='utf-8') as f:
    all_reqs = f.read().split('\n')

install_requires = [x.strip() for x in all_reqs if 'git+' not in x]
dependency_links = [x.strip().replace('git+', '')
                    for x in all_reqs if x.startswith('git+')]

setup(
    name='dogma',
    version=__version__,
    description=DESCRIPTION,
    long_description=readme,
    url='https://github.com/griffinclausen/dogma',
    license='MIT',
    author='Griffin Clausen',
    author_email='dogma.maintainer@gmail.com',
    keywords='',
    entry_points='''
        [console_scripts]
        genetic_code=dogma.genetic_codes:test
    ''',
    packages=find_packages(exclude=['docs', 'tests*']),
    include_package_data=True,
    install_requires=[
        'pandas'
        ],
    dependency_links=dependency_links
)
