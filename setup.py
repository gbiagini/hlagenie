#!/usr/bin/env python

from setuptools import setup

with open("README.md") as readme_file:
    readme = readme_file.read()

with open("requirements.txt") as requirements_file:
    requirements = requirements_file.read().split("\n")

with open("requirements-tests.txt") as requirements_file:
    test_requirements = requirements_file.read().split("\n")

setup(
    name="hlagenie",
    version="0.4.0",
    description="Sequence handing for HLA with Python",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Giovanni Biagini",
    author_email="dbiagini@tulane.edu",
    url="https://github.com/gbiagini/hlagenie",
    packages=[
        "hlagenie",
    ],
    provides=["hlagenie"],
    scripts=[
        "scripts/hlagenie",
        "scripts/hlagenie-match",
    ],
    install_requires=requirements,
    license="LGPL 3.0",
    zip_safe=False,
    keywords="hlagenie",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.11",
    ],
    test_suite="tests",
    tests_require=test_requirements,
    include_package_data=True,
)
