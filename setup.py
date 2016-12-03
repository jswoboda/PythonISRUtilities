#!/usr/bin/env python
import subprocess
from setuptools import setup

try:
    subprocess.call(['conda','install','--yes','--file','requirements.txt'])
except Exception as e:
    pass


config = {
    'description': 'ISR constants',
    'author': 'John Swoboda',
    'url': 'https://github.com/jswoboda/PythonISRUtilities.git',
    'version': '1',
    'install_requires': [],
    'packages': ['isrutilities'],
    'name':      'isrutilities'
}

setup(**config)
