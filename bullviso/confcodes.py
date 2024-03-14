"""
BULLVISO
Copyright (C) 2022  Conor D. Rankine

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software 
Foundation, either Version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this program.  If not, see <https://www.gnu.org/licenses/>.
"""

###############################################################################
############################### LIBRARY IMPORTS ###############################
###############################################################################

from itertools import permutations
from bullviso.utils import rotate_tuple

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def gen_confcodes(n_func_sites: int, unique_only: bool = True) -> set:

    confcode_ = tuple(
        0 if i < 10 - n_func_sites else 1 for i in range(10)
    )

    confcodes = set(
        confcode for confcode in permutations(confcode_)
    )

    if unique_only:
        return filter_confcodes(confcodes)
    else:
        return confcodes

def gen_equiv_confcodes(confcode: tuple) -> set:

    equiv_confcodes = set(
        rotate_tuple(confcode[:-1], 3 * i) + confcode[-1:] for i in range(3)
    )

    return equiv_confcodes

def filter_confcodes(confcodes: set) -> set:

    unique_confcodes = set()

    for confcode in sorted(confcodes):
        equiv_confcodes = gen_equiv_confcodes(confcode)
        if unique_confcodes.isdisjoint(equiv_confcodes):
            unique_confcodes.add(confcode)

    return unique_confcodes

