#!/usr/bin/env python
"""Install the package."""

from distutils.extension import Extension

from setuptools import setup

setup(
    use_scm_version=True,
    ext_modules=[
        Extension(
            "read_acq.decode",
            ["src/read_acq/decode.c"],
            extra_compile_args=["-Ofast", "-Wall"],
            libraries=["m"],
        ),
    ],
)
