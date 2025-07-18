#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='methylartist',
    version='1.5.2',
    author='Adam Ewing',
    author_email='adam.ewing@gmail.com',
    description=("Tools for parsing and plotting nanopore methylation data"),
    license='MIT',
    url='https://github.com/adamewing/methylartist',
    download_url='https://github.com/adamewing/methylartist/archive/refs/tags/1.5.2.tar.gz',
    scripts=['methylartist'],
    packages=find_packages(),
    install_requires = [
        'cython',
        'pysam',
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
        'seaborn',
        'bx-python',
        'scikit-learn',
        'tqdm',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],

)
