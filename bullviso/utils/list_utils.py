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

from typing import TypeVar

# =============================================================================
#                                   GLOBALS
# =============================================================================

T = TypeVar('T')

# =============================================================================
#                                  CONSTANTS
# =============================================================================

PADDING_LEFT = 'left'
PADDING_RIGHT = 'right'

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

def all_same_length(
    *input_lists: list
) -> bool:
    """
    Validates that all the input lists have the same length.

    Args:
        *input_lists (list): Variable number of input lists to validate.

    Returns:
        bool: `True` if all the input lists have the same length, else `False`.
    """
    
    if not input_lists:
        return True
    else:
        return len({len(input_list) for input_list in input_lists}) == 1

def pad_list(
    input_list: list[T],
    length: int,
    direction: str = PADDING_RIGHT
) -> list[T | int]:
    """
    Zero-pads a list up to the specified length.

    Args:
        input_list (list[T]): Input list.
        length (int): Length target for the padded list; lists longer than this
            are returned as-is, i.e., without padding.
        direction (str, optional): Direction of padding; use `PADDING_RIGHT`
            for right padding and `PADDING_LEFT` for left padding. Defaults to
            'PADDING_RIGHT'.

    Raises:
        ValueError: If `length` is negative (< 0).
        ValueError: If `direction` is not valid, i.e., if it is not one of
            {'left', 'right'}. 

    Returns:
        list[T | int]: List containing the elements of the input list, padded
            with zeros up to the specified length. If the input list is already
            at or above the length target, a copy is returned as-is.
    """
    
    if length < 0:
        raise ValueError(
            f'`length should be non-negative (>= 0): got {length}'
        )

    if len(input_list) < length:
        padding = [0] * (length - len(input_list))
        if direction == PADDING_LEFT:
            return padding + input_list
        elif direction == PADDING_RIGHT:
            return input_list + padding
        else:
            raise ValueError(
                f'`direction` should be one of {{{PADDING_LEFT!r}, '
                f'{PADDING_RIGHT!r}}}; got {direction!r}'
            )
    else:
        return input_list.copy()

def roll(
    input_list: list[T],
    n: int
) -> list[T]:
    """
    Rolls/rotates a list forwards (if `n` is positive) or backwards (if `n` is
    negative) by `n` places.

    Args:
        input_list (list[T]): Input list.
        n (int): Number of places to roll/rotate the input list by:
            - positive values = forward rotation (elements move right).
            - negative values = backward rotation (elements move left).

    Returns:
        list[T]: List containing the elements of the input list rolled/
            rotated forwards or backwards by `n` places.
    """
    
    if not input_list:
        return []
    
    n = n % len(input_list)

    return input_list[-n:] + input_list[:-n]

def maxdepth(
    input_data: object
) -> int:
    """
    Returns the maximum depth of any nested (sub)lists in the input data.

    This function recursively traverses the input data to find the deepest
    level of list nesting. Elements of the input data that are not lists, e.g.,
    other data types, have a depth of zero.

    Args:
        input_data (object): Input; can be list or any other data type.

    Returns:
        int: Maximum depth of any nested (sub)lists in the input.
    """
    
    if not isinstance(input_data, list):
        return 0
    else:
        return 1 + max(
            (maxdepth(item) for item in input_data), default = 0
        )

# =============================================================================
#                                     EOF
# =============================================================================