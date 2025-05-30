#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev

import io
import os
import sys
from shutil import rmtree
from glob import glob
from pathlib import Path
from setuptools import find_packages, setup, Command
try:
    from pybind11.setup_helpers import Pybind11Extension
except ImportError:
    from setuptools import Extension as Pybind11Extension
    
def glob_fix(package_name, glob):
    # this assumes setup.py lives in the folder that contains the package
    package_path = Path(f'./{package_name}').resolve()
    return [str(path.relative_to(package_path)) 
            for path in package_path.glob(glob)]

# Package meta-data.
NAME = 'homolig'
DESCRIPTION = 'A simple tool to identify sequence similarity in immune receptor repertoires.'
URL = 'https://github.com/edavis71/homolig'
EMAIL = 'agirgis3@jhu.edu'
AUTHOR = 'Alexander Girgis'
REQUIRES_PYTHON = '>=3.8.0, <3.13.0'
VERSION = '1.0.0'

# What packages are required for this module to be executed?
#REQUIRED = [
#     # os, sys, and itertools are modules included with Python and should not be listed
#     'numpy==1.22.4', 'pandas==1.4.3', 'more-itertools==8.10.0', 'bio==1.3.9', 'anndata==0.8.0', 'tqdm==4.64.0', #'scipy==1.9.0',
#     'pybind11==2.10.0','scikit-learn==1.1.1'
#]

#REQUIRED = [
#     # os, sys, and itertools are modules included with Python and should not be listed
#     'numpy', 'pandas', 'more_itertools',
#     'Bio', 'anndata', 'tqdm', 'scipy','pybind11>=2.2'
#]

#Note: Specifying exact versions here is probably more stringent than necessary, but these versions have been validated to work with python 3.11. 
REQUIRED = [ 
'anndata==0.8.0',
'array_api_compat==1.9.1',
'bio==1.3.9',
'biopython==1.79',
'biothings-client==0.2.6',
'certifi==2024.8.30',
'charset-normalizer==3.4.0',
'exceptiongroup==1.2.2',
'h5py==3.7.0',
'homolig==1.0.0',
'idna==3.10',
'llvmlite==0.39.0',
'matplotlib==3.5.2',
'more-itertools==10.5.0',
'mygene==3.2.2',
'natsort==8.4.0',
'numpy==1.22.4',
'packaging==24.1',
'pandas==1.4.3',
'platformdirs==4.3.6',
'pybind11==2.13.6',
'python-dateutil==2.9.0.post0',
'pytz==2024.2',
'requests==2.32.3',
'scanpy==1.9.1',
'scikit-learn==1.1.1',
'scipy==1.13.1',
'six==1.16.0',
'tqdm==4.66.6',
'tzdata==2024.2',
'urllib3==2.2.3'

]


# What packages are optional?
EXTRAS = {
    # 'fancy feature': ['django'],
}

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION


class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')

        sys.exit()

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

ext_modules = [
    Pybind11Extension(
        "homoligcpp",
        ['homoligcpp/aamatrix.cpp',
        'homoligcpp/aavector.cpp',
        'homoligcpp/homolig.cpp',
        'homoligcpp/homolig_bind.cpp'],
        include_dirs=[
            get_pybind_include(),
            get_pybind_include(user=True),
            'homoligcpp/'
        ],
        language="c++",
        cxx_std=17),
]

# Where the magic happens:
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],

    install_requires=REQUIRED,
    extras_require=EXTRAS,
    #include_package_data=True,
    package_data={'homolig': ['data/align_matrices/*', 
                  *glob_fix('homolig', 'data/fastas/**/*'),
                  'data/imgt_genedb_full.csv',
                  'data/mapper.csv']},
    license='MIT',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        #'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy'
    ],
    # $ setup.py publish support.
    cmdclass={
        'upload': UploadCommand,
    },
    ext_modules=ext_modules,
    setup_requires=['pybind11>=2.2']
)
