#!/usr/bin/env python

from setuptools import setup, find_packages
import pathlib

# Taken directly from
# https://github.com/pypa/sampleproject/blob/main/setup.py
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')


setup(
    name='smsn',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='v1.0.0',
    url='https://github.com/EMeyerLab/SMSN',
    author='DELEVOYE Guillaume',
    package_dir={'': 'src'},  # Optional
    packages=find_packages(where='src'),
    include_package_data=True,
    license=open('LICENSE.txt').read(),
    entry_points={'console_scripts': [
        "smsn=smsn.launchers.smsn_launcher:main",
    ]},
    setup_requires=['pytest-runner']
)

