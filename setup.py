#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import os
import re
# from setuptools import setup
from numpy.distutils.core import setup, Extension
from os.path import dirname
from os.path import join
from setuptools import find_packages


def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()


setup(
    name='read_acq',
    version='0.2.0',
    license='MIT license',
    description='Read ACQ file types',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.md')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.md'))
    ),
    author='EDGES Collaboration',
    author_email='steven.g.murray@asu.edu',
    url='https://github.com/edges-collab/read_acq',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    install_requires=[
        'numpy',
        'tqdm',
        'scipy',
        'click'
    ],
    extras_require={
        'dev': [
            'bump2version',
        ],
        'all': [
            'h5py'
        ]
    },
    ext_modules=[
        Extension(
            'read_acq.decode',
            ['read_acq/decode.c'],
            extra_compile_args=['-Ofast', '-Wall'],
            libraries=['m', ],
        ),
    ],
    entry_points={
        'console_scripts': [
            'acq = read_acq.cli:main',
        ]
    },
)
