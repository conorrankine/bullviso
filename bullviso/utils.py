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

def count_list_elements(
    input_list: list[Any]
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
    input_list: list[Any],
    repeats: Union[int, list[int]]
) -> list[Any]:
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

def unique_elements(
    input_list: list[Any]
) -> list[Any]:
    """
    Returns the unique elements in a list, preserving their order of occurance.

    Args:
        input_list (list[Any]): Input list.

    Returns:
        list[Any]: List containing the unique elements of the input list.
    """
    
    return [item for item in dict.fromkeys(input_list)]

def pad_list(
    input_list: list[Any],
    length: int,
    direction: str = 'right'
) -> list[Any]:
    """
    Zero-pads a list up to the specified length.

    Args:
        input_list (list[Any]): Input list.
        length (int): Integer specifying the length of the list post-padding;
            lists longer than this are returned as-is without padding.
        direction (str, optional): String specifying the direction of the
            padding; options are `left` or `right`. Defaults to 'right'.

    Raises:
        ValueError: If `direction` is not valid; i.e. if it is not one of
            either `left` or `right`.

    Returns:
        list[Any]: List containing the elements of the input list, padded with
            zeros up to the specified length.
    """
    
    if len(input_list) < length:
        padding = [0] * (length - len(input_list))
        if direction == 'left':
            return padding + input_list
        elif direction == 'right':
            return input_list + padding
        else:
            raise ValueError(
                '`direction` should be either `left` or `right`'
            )
    else:
        return input_list

def iterate_and_index(
    input_list: list[Union[Any, list[Any]]]
) -> Generator[tuple[int, Any], None, None]:
    """
    Iterates through a list that contains items and/or nested (sub)lists of
    items, yielding tuples comprising the index of each item in the outer list
    and either i) the item, or ii) each of the items in the nested (sub)list.

    Args:
        input_list (list[Union[Any, list[Any]]]): Input list.

    Raises:
        ValueError: If `input_list` has a maximum depth greater than 2, i.e.
            if the items in any nested (sub)list in `input_list` are
            themselves nested (sub)lists.

    Yields:
        Generator[tuple[int, Any], None, None]: Generator that yields tuples
            comprising the index of each item in the outer list and either i)
            the item, or ii) each of the items in the nested (sub)list.
    """
    
    if maxdepth(input_list) > 2:
        raise ValueError(
            '`input_list` cannot have a maximum depth greater than 2'
        )

    for idx, item in enumerate(input_list):
        if isinstance(item, list):
            for item_ in item:
                yield idx, item_
        else:
            yield idx, item

def maxdepth(
    input_list: Any
) -> int:
    """
    Returns the maximum depth of a list that contains nested (sub)lists.

    Args:
        input_list (list): Input list.

    Returns:
        int: Maximum depth of the input list.
    """
    
    if not isinstance(input_list, list):
        return 0
    else:
        return 1 + max(
            (maxdepth(item) for item in input_list), default = 0
        )
    
def roll(
    input_list: list[Any],
    n: int
) -> list[Any]:
    """
    Rolls/rotates a list forwards (if `n` is positive) or backwards (if `n` is
    negative) by `n` places.

    Args:
        input_list (list[Any]): Input list.
        n (int): Integer specifying the number of places to roll/rotate the
            input list forwards or backwards by.

    Returns:
        list[Any]: List containing the elements of the input list
            rolled/rotated forwards or backwards by `n` places.
    """
    
    return input_list[-n:] + input_list[:-n]
