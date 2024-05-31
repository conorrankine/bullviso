###############################################################################
############################### LIBRARY IMPORTS ###############################
###############################################################################

from setuptools import setup
from setuptools import find_packages
import bullviso

###############################################################################
#################################### SETUP ####################################
###############################################################################

setup(
    name = 'bullviso',
    version = bullviso.__version__,
    author = 'Conor D. Rankine',
    author_email = 'conor.rankine@york.ac.uk',
    url = 'https://gitlab.com/conorrankine/bullviso',
    description = ('A package for generating configuration codes and '
        ' molecular structures of functionalised bullvalenes'),
    licence = 'GPL',
    packages = find_packages(),
    entry_points = {
        'console_scripts': ['bullviso = bullviso.cli:main']
    },
    install_requires = [
        'rdkit',
        'tqdm'
    ]
)
