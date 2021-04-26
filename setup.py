#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='smsn',
    version='v0.2-alpha',
    author='DELEVOYE Guillaume',
    license=open('LICENSE.txt').read(),
    packages=find_packages("."),
    python_requires='3.7',
    package_data={'smsn': ['resources/']},
    entry_points={'console_scripts': [
        "smsn = smsn.launchers.smsn_launcher:main",
    ]},
)
