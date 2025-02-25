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
############################### LIBRARY IMPORTS ###############################
###############################################################################

from typing import Union, Generator, Any

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

def repeat_list_elements(
    input_list: list,
    repeats: Union[int, list[int]]
):
    """
    Repeats elements in a list.

    Args:
        input_list (list): Input list.
        repeats (Union[int, list[int]]): Integer or list of integers indicating
            the number of times that each (corresponding) element in the input
            list is to be repeated.

    Raises:
        ValueError: If `repeats` is a list and is not the same length as
            `input_list`, or if any of the elements of `repeats` are integers
            less than or equal to zero.

    Returns:
        _type_: List containing the elements of the input list repeated.
    """
    
    if isinstance(repeats, int):
        repeats = [repeats] * len(input_list)
    elif len(repeats) != len(input_list):
        raise ValueError(
            '`repeats` and `input_list` should have the same length'
        )

    if any(repeat <= 0 for repeat in repeats):
        raise ValueError(
            '`repeats` cannot contain integers less than or equal to zero'
        )
    
    input_list_repeated = []
    for i, repeat in zip(input_list, repeats):
        input_list_repeated.extend([i] * repeat)

    return input_list_repeated

def iterate_and_index(
    input_list: list[Union[Any, list[Any]]]
) -> Generator[tuple[int, Any], None, None]:
    """
    Iterates through a list that contains items and/or nested (sub)lists of
    items, yielding tuples comprising the index of each item in the outer list
    and either i) the item, or ii) each of the items in the nested (sub)list.

    Args:
        input_list (list[Union[Any, list[Any]]]): Input list.

    Yields:
        Generator[tuple[int, Any], None, None]: Generator that yields tuples
            comprising the index of each item in the outer list and either i)
            the item, or ii) each of the items in the nested (sub)list.
    """
    
    for idx, item in enumerate(input_list):
        if isinstance(item, list):
            for item_ in item:
                yield idx, item_
        else:
            yield idx, item
