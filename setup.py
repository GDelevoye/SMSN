#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='smsn',
    version='0.1',
    author='DELEVOYE Guillaume',
    license=open('LICENSE.txt').read(),
    packages=find_packages("."),
    python_requires='>=3.8',
    package_data={'smsn': ['resources/']},
    entry_points={'console_scripts': [
        "smsn = smsn.launchers.smsn_launcher:main",
    ]},
)
