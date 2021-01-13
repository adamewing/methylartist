#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='tmnt',
    version='1.0',
    author='Adam Ewing',
    author_email='adam.ewing@gmail.com',
    description=("Tools for parsing and plotting nanopore methylation data aimed at TEs"),
    license='MIT',
    url='https://github.com/adamewing/tmnt',
    scripts=['tmnt'],
    packages=find_packages(),
    install_requires = [
        'cython',
        'pysam',
        'scikit-bio',
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
        'seaborn'
    ]

)
