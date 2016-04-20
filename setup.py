"""
SolvationToolkit: Tools for setting up mixtures of small organic molecules for AMBER, GROMACS, OpenMM and others.

You can install SolvationTookit via

python setup.py install

"""

import sys
from setuptools import setup, find_packages

if sys.argv[-1] == 'setup.py':
    print("To install, run 'python setup.py install'")
    print()

if sys.version_info[:2] < (2, 7):
    print("SolvationToolkit requires Python 2.7 or later (%d.%d detected)." %
          sys.version_info[:2])
    sys.exit(-1)


##########################
VERSION = "0.2.1"
ISRELEASED = True
__version__ = VERSION
##########################


descr = """
SolvationToolkit is a simple toolkit for setting up input files for simulations of solutions of small organic molecules for various simulation packages.
"""

data = {'solvationtoolkit':['test/*.py'] }

setup(
    name                 = 'solvationtoolkit', 
    version              = __version__, 
    description          = 'Solvation Toolkit',
    long_description     = descr,
    url                  = 'https://github.com/mobleylab/solvationtoolkit',
    author               = 'Gaetano Calabro and David Mobley',
    author_email         = 'gcalabro -at- uci.edu',
    license              = 'LGPL',
    platforms            = ['Linux-64', 'Mac OSX-64', 'Unix-64'],
    packages             = find_packages(),  
    package_data         = data,
    include_package_data = True,
      
    #entry_points         = {'console_scripts':['start=scripts/start.py']},
    zip_safe             = False

)
