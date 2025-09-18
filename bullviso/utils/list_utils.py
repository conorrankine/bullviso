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

# =============================================================================
#                               LIBRARY IMPORTS
# =============================================================================

from typing import Union, Generator, Any

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

def all_same_length(
    *input_lists: list[Any] 
) -> bool:
    """
    Validates that all of the input lists have the same length.

    Returns:
        bool: True if all of the input lists have the same length, else False.
    """
    
    if not input_lists:
        return True
    else:
        return len({len(input_list) for input_list in input_lists}) == 1

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

# =============================================================================
#                                     EOF
# =============================================================================