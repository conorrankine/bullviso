###############################################################################
############################### LIBRARY IMPORTS ###############################
###############################################################################

from setuptools import setup
from setuptools import find_packages
import pybull

###############################################################################
#################################### SETUP ####################################
###############################################################################

setup(
    name = 'pybull',
    version = pybull.__version__,
    author = 'Conor D. Rankine',
    author_email = 'conor.rankine@york.ac.uk',
    url = 'https://gitlab.com/conor.rankine/pybull',
    description = ('A package for generating configuration codes and '
        ' molecular structures of functionalised bullvalenes'),
    licence = 'GPL',
    packages = find_packages(),
    scripts = ['bin/pybull'],
    install_requires = [
        'numpy',
        'rdkit'
    ]
)
