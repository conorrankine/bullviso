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

from itertools import permutations
from typing import Generator

###############################################################################
################################### CLASSES ###################################
###############################################################################

class BVBarcode:

    def __init__(self, barcode: tuple, grouped_barcode: tuple = None):
        """
        Instantiates a bullvalene isomer barcode.

        Args:
            barcode (tuple): A tuple of 10 integers in the range 0-9 inclusive
                representing the bullvalene isomer barcode.
            grouped_barcode (tuple): A tuple of 10 integers in the range 0-9
                inclusive representing the grouped bullvalene isomer barcode;
                e.g. for the barcode (0,0,0,0,0,0,1,2,3,4) and the grouped
                barcode (0,0,0,0,0,0,1,1,1,2), the implication is that
                substituents 1, 2, and 3 are equivalent, while substituent
                4 is unique (i.e. the bullvalene has two unique substituents).

        Raises:
            ValueError: If the barcode or grouped barcode i) are not of length
                10, ii) contain non-integer elements, or iii) contain integers
                less than 0 or greater than 9.
        """

        if len(barcode) == 10 and all(
            isinstance(x, int) and 0 <= x <= 9 for x in barcode
        ):        
            self._barcode = barcode
        else:
            raise ValueError(
                f'{barcode} is not a valid barcode'
            )
        
        if grouped_barcode is not None:
            if len(grouped_barcode) == 10 and all(
                isinstance(x, int) and 0 <= x <= 9 for x in grouped_barcode
            ):        
                self._grouped_barcode = grouped_barcode
            else:
                raise ValueError(
                    f'{grouped_barcode} is not a valid grouped barcode'
                )
        else:
            self._grouped_barcode = barcode

    def __hash__(self) -> int:
        """
        Returns a hash for the `BVBarcode` instance; two `BVBarcode` instances
        will have the same hash if they have the same grouped barcode
        representation (i.e. the same `self.grouped_barcode`).

        Returns:
            int: Hash for the `BVBarcode` instance.
        """
        
        return hash(self.grouped_barcode)
    
    def __eq__(self, barcode: object) -> bool:
        """
        Returns the result of an equality test between the `BVBarcode` instance
        and another `BVBarcode` instance (`barcode`); two `BVBarcode` instances
        are considered equal if they have the same grouped barcode
        representation (i.e. the same `self.grouped_barcode`).

        Args:
            barcode (BVBarcode): `BVBarcode` instance to test for equality.

        Returns:
            bool: True if both `BVBarcode` instances have the same grouped
                barcode representation, else False.
        """

        if isinstance(barcode, BVBarcode):
            return self.grouped_barcode == barcode.grouped_barcode
        else:
            return False

    @property
    def barcode(self) -> tuple:
        """
        Returns the bullvalene isomer barcode as a tuple.

        Returns:
            tuple: The bullvalene isomer barcode as a tuple.
        """

        return self._barcode
    
    @property
    def grouped_barcode(self) -> tuple:
        """
        Returns the grouped bullvalene isomer barcode as a tuple.

        Returns:
            tuple: The grouped bullvalene isomer barcode as a tuple.
        """

        return self._grouped_barcode
    
    def permutations(self) -> Generator['BVBarcode', None, None]:
        """
        Generates all permutations of the bullvalene isomer barcode, yielding a
        `BVBarcode` instance for each permutation; the `BVBarcode` instances
        that the method yields are instantiated with the permuted attributes 
        `barcode` and `grouped_barcode`.

        Yields:
            BVBarcode: `BVBarcode` instances for each bullvalene isomer barcode
            permutation.
        """

        for permutation_idx in permutations(
            tuple(i for i in range(len(self.barcode)))
        ):
            barcode = tuple(
                self.barcode[i] for i in permutation_idx
            )
            grouped_barcode = tuple(
                self.grouped_barcode[i] for i in permutation_idx
            )
            yield BVBarcode(
                barcode,
                grouped_barcode = grouped_barcode
            )
