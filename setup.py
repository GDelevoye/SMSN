#!/usr/bin/env python

from setuptools import setup, Extension, find_packages
import os
import sys

setup(
    name='smsn',
    version='0.1',
    author='DELEVOYE Guillaume',
    license=open('LICENSES.txt').read(),
    packages=find_packages("."),
    python_requires='>=3.8',
    package_data={'smsn': ['resources/']},
    install_requires=[
        'numpy',
        'joblib',
        'h5py',
        'pandas',
        'tqdm',
        'h5py',
        'numpy'
    ],
    entry_points={'console_scripts': [
        "smsn = smsn.launchers.smsn_launcher:main",
    ]},
)
