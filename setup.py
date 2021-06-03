#!/usr/bin/env python
"""Install the package."""

import re
from os.path import dirname, join

from numpy.distutils.core import Extension
from setuptools import find_packages, setup


def read(*names, **kwargs):
    """Read a file to a string."""
    return open(
        join(dirname(__file__), *names), encoding=kwargs.get("encoding", "utf8")
    ).read()


setup(
    name="read_acq",
    license="MIT license",
    description="Read ACQ file types",
    long_description="%s\n%s"
    % (
        re.compile("^.. start-badges.*^.. end-badges", re.M | re.S).sub(
            "", read("README.rst")
        ),
        re.sub(":[a-z]+:`~?(.*?)`", r"``\1``", read("CHANGELOG.rst")),
    ),
    author="EDGES Collaboration",
    author_email="steven.g.murray@asu.edu",
    url="https://github.com/edges-collab/read_acq",
    packages=find_packages("src"),
    package_dir={"": "src"},
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    setup_requires=["setuptools_scm"],
    use_scm_version=True,
    install_requires=["numpy", "tqdm", "scipy", "click"],
    extras_require={
        "dev": [
            "pre-commit",
            "pytest>=5<6",
            "pytest-cov",
            "edges-io",
        ],
        "h5": [
            "edges-io @ git+git://github.com/edges-collab/edges-io.git",
        ],
    },
    ext_modules=[
        Extension(
            "read_acq.decode",
            ["src/read_acq/decode.c"],
            extra_compile_args=["-Ofast", "-Wall"],
            libraries=[
                "m",
            ],
        ),
    ],
    entry_points={
        "console_scripts": [
            "acq = read_acq.cli:main",
        ]
    },
)
