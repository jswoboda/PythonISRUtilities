#!/usr/bin/env python
req = ['nose','numpy','scipy','tables','pathlib2']

from setuptools import setup, find_packages

config = {
    'description': 'ISR constants',
    'author': 'John Swoboda',
    'url': 'https://github.com/jswoboda/PythonISRUtilities.git',
    'version': '1',
    'install_requires': req,
    'packages': find_packages(),
    'name':     'isrutilities',
    'python_requires': '>=2.7',
}

setup(**config)
