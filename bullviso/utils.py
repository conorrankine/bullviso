"""
BULLVISO
Copyright (C) 2024  Conor D. Rankine

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
################################## FUNCTIONS ##################################
###############################################################################

def rotate_tuple(t: tuple, n: int) -> tuple:
    # rotates/rolls the tuple `t` forward by `n` positions, e.g. `rotate_tuple(
    # (1,2,3,4,5,6), 2)` returns the rotated/rolled tuple `(5,6,1,2,3,4)`
    
    return t[-n:] + t[:-n]

def count_list_elements(
    input_list: list
) -> int:
    """
    Counts the number of elements in a list, including the elements of any
    nested (sub)lists recursively.

    Args:
        input_list (list): Input list.

    Returns:
        int: Number of elements in the input list.
    """
    
    count = 0

    for i in input_list:
        if isinstance(i, list):
            count += count_list_elements(i)
        else:
            count += 1

    return count
