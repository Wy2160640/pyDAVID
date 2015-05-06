#!/usr/bin/env python

import os
import sys

from distutils.core import setup
#from ez_setup import use_setuptools
#use_setuptools()
#from setuptools import setup

# convenience is king
opj = os.path.join

# setup and install
setup(name='pyDAVID',
      version='1.0',
      author='Adam Labadorf',
      author_email='alabadorf@gmail.com',
      py_modules=['pyDAVID'],
     )
