# Author: Ivan E. Cao-Berg (icaoberg@scs.cmu.edu)
#
# Copyright (C) 2011-2014 Murphy Lab
# Lane Center for Computational Biology
# School of Computer Science
# Carnegie Mellon University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# For additional information visit http://murphylab.web.cmu.edu or
# send email to murphy@cmu.edu

import os
from setuptools import setup

#load the current version number of the package
exec(compile(open('VERSION').read(),'VERSION', 'exec'))

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name = 'ricerca',
      version = __version__,
      description = ('Content-based search based using FALCON: '
      	'Feedback Adaptive Loop for Content-Based Retrieval'),
      long_description=read('README'),
      author = 'Ivan Cao-Berg',
      author_email = 'icaoberg@andrew.cmu.edu',
      install_requires = [
        'numpy>=1.4.1',
        'scipy>=0.7.2',
        'numexpr>=2.3',
        'cython>=0.20',
        'tables>=3.1.0',
        ],
      url = 'http://murphylab.web.cmu.edu/software/ricerca',
      classifiers=[
      	'Programming Language :: Python', 
      	'Intended Audience :: Science/Research',
      	'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering',
        'Topic :: Database',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: OS Independent',
        'Programming Language :: Python'],
      py_modules=['ricerca.content'])
