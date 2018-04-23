#!/usr/bin/env python

"""Setup script."""

import os
from setuptools import setup, find_packages

# Long description loaded from the README
README = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/README.rst'

# Get requirements
REQUIREMENTS = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/requirements.txt'

# Version of the soft
VERSION = "0.0"

# Package name
NAME = 'drptools'

# Packages (subdirectories in drptools/)
PACKAGES = find_packages()

CLASSIFIERS = ['Development Status :: 3 - Alpha',
               'Intended Audience :: Science/Research',
               'Topic :: Software Development :: Build Tools',
               'License :: OSI Approved :: MIT License',
               'Programming Language :: Python :: 3',
               'Topic :: Scientific/Engineering :: Astronomy']

setup(name=NAME,
      version=VERSION,
      description=("Some tools to explore an LSST Data Release Processing output directory"),
      license="MIT",
      classifiers=CLASSIFIERS,
      url="https://github.com/nicolaschotard/drptools",
      author="Nicolas Chotard",
      author_email="nchotard@in2p3.fr",
      packages=PACKAGES,
      long_description=open(README).read(),
     )
