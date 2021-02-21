#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='methylartist',
    version='0.1',
    author='Adam Ewing',
    author_email='adam.ewing@gmail.com',
    description=("Tools for parsing and plotting nanopore methylation data"),
    license='MIT',
    url='https://github.com/adamewing/methylartist',
    scripts=['methylartist'],
    packages=find_packages(),
    install_requires = [
        'cython',
        'pysam',
        'scikit-bio',
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
        'seaborn',
        'ont-fast5-api'
    ]

)
