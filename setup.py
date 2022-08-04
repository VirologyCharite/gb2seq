#!/usr/bin/env python

from setuptools import setup
from pathlib import Path


# Modified from http://stackoverflow.com/questions/2058802/
# how-can-i-get-the-version-defined-in-setup-py-setuptools-in-my-package
def version():
    import re

    init = Path("gb2seq") / "__init__.py"
    with open(init) as fp:
        initData = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]+)['\"]", initData, re.M)
    if match:
        return match.group(1)
    else:
        raise RuntimeError(f"Unable to find version string in {init!r}.")


setup(
    name="gb2seq",
    version=version(),
    packages=["gb2seq"],
    url="https://github.com/virologycharite/gb2seq",
    download_url="https://github.com/virologycharite/gb2seq",
    author="Terry C. Jones",
    author_email="terence.jones@charite.de",
    keywords=["GenBank", "genetic sequences"],
    scripts=[
        "bin/describe-feature.py",
        "bin/describe-genome.py",
        "bin/describe-site.py",
    ],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    description=(
        "Python library and scripts for working with SARS-CoV-2 " "sequences."
    ),
    long_description=("See https://github.com/virologycharite/gb2seq for details."),
    license="MIT",
    install_requires=["dark-matter>=4.0.21"],
)
