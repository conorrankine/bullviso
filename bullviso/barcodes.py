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
from utils import rotate_tuple
from typing import Generator, Tuple

###############################################################################
################################### CLASSES ###################################
###############################################################################

class BVBarcode:

    def __init__(
        self,
        barcode: tuple[int, ...],
        grouped_barcode: tuple[int, ...] = None,
        canonicalize: bool = False
    ):
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
            canonicalize (bool): If `True`, the `barcode` and `grouped_barcode`
                attributes are converted to the canonical representation (i.e.
                the smallest grouped bullvalene isomer barcode
                lexicographically), else (if `False`) the `barcode` and
                `grouped_barcode` attributes are stored as-is.

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

        if canonicalize:
            self.canonicalize()
            self._canonical_barcode = self.barcode
            self._canonical_grouped_barcode = self.grouped_barcode
        else:
            canonical_equivalent = self.canonicalize(inplace = False)
            self._canonical_barcode = (
                canonical_equivalent.barcode
            )
            self._canonical_grouped_barcode = (
                canonical_equivalent.grouped_barcode
            )

    def __hash__(
        self
    ) -> int:
        """
        Returns a hash for the `BVBarcode` instance; two `BVBarcode` instances
        will have the same hash if their canonical representations (i.e. the
        smallest grouped bullvalene isomer barcodes lexicographically) are the
        same.

        Returns:
            int: Hash for the `BVBarcode` instance.
        """
        
        return hash(
            self.canonical_grouped_barcode
        )
    
    def __eq__(
        self,
        barcode: 'BVBarcode'
    ) -> bool:
        """
        Returns the result of an equality test between the `BVBarcode` instance
        and another `BVBarcode` instance (`barcode`); two `BVBarcode` instances
        are considered equal if their canonical representations (i.e. the
        smallest grouped bullvalene isomer barcodes lexicographically) are the
        same..

        Args:
            barcode (BVBarcode): `BVBarcode` instance to test for equality.

        Returns:
            bool: True if both `BVBarcode` instances have the same grouped
                barcode representation, else False.
        """

        if isinstance(barcode, BVBarcode):
            return (
                self.canonical_grouped_barcode
                == barcode.canonical_grouped_barcode
            )
        else:
            return False

    @property
    def barcode(
        self
    ) -> tuple[int, ...]:
        """
        Returns the bullvalene isomer barcode as a tuple.

        Returns:
            tuple: The bullvalene isomer barcode as a tuple.
        """

        return self._barcode
    
    @property
    def grouped_barcode(
        self
    ) -> tuple[int, ...]:
        """
        Returns the grouped bullvalene isomer barcode as a tuple.

        Returns:
            tuple: The grouped bullvalene isomer barcode as a tuple.
        """

        return self._grouped_barcode
    
    @property
    def canonical_barcode(
        self
    ) -> tuple[int, ...]:
        """
        Returns the canonical bullvalene isomer barcode as a tuple.

        Returns:
            tuple: The canonical bullvalene isomer barcode as a tuple.
        """

        return self._canonical_barcode

    @property
    def canonical_grouped_barcode(
        self
    ) -> tuple[int, ...]:
        """
        Returns the canonical grouped isomer barcode as a tuple.

        Returns:
            tuple: The canonical grouped bullvalene isomer barcode as a tuple.
        """

        return self._canonical_grouped_barcode
    
    def permutations(
        self
    ) -> Generator['BVBarcode', None, None]:
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

    def equivalents(
        self
    ) -> Tuple['BVBarcode']:
        """
        Returns a tuple of three `BVBarcode` instances that are equivalent
        by rotation around the threefold symmetry axis of the bullvalene; the
        `BVBarcode` instances that the method returns are instantiated with
        the threefold-rotated attributes `barcode` and `grouped_barcode`.

        Returns:
            Tuple[BVBarcode]: Tuple of three `BVBarcode` instances that are
            equivalent by rotation around the threefold symmetry axis of the
            bullvalene.
        """

        equivalents = []

        for i in range(3):
            barcode = (
                rotate_tuple(
                    self.barcode[:-1], 3 * i
                ) + self.barcode[-1:]
            )
            grouped_barcode = (
                rotate_tuple(
                    self.grouped_barcode[:-1], 3 * i
                ) + self.grouped_barcode[-1:]
            )
            equivalents.append(
                BVBarcode(
                    barcode,
                    grouped_barcode = grouped_barcode
                )
            )

        return tuple(equivalents)
    
    def canonicalize(
        self,
        inplace: bool = True
    ) -> 'BVBarcode':
        """
        Canonicalizes the bullvalene isomer barcode; updates the `barcode` and
        `grouped_barcode` attributes of the `BVBarcode` instance to the
        canonical representation (i.e. the smallest grouped bullvalene isomer
        barcode lexicographically) if `inplace = True`, else returns a new
        `BVBarcode` instance with the `barcode` and `grouped_barcode`
        attributes corresponding to the canonical representation.

        Args:
            inplace (bool): If `True`, the `barcode` and `grouped_barcode`
                attributes are updated in place, else (if `False`) a new
                `BVBarcode` instance is returned with the `barcode` and
                `grouped_barcode` attributes corresponding to the canonical
                representation.

        Returns:
            BVBarcode: Canonicalized `BVBarcode` instance if `inplace = False`, 
                else `None` if `inplace = True`.
        """

        canonical_equivalent = min(
            self.equivalents(),
            key = lambda equivalent: equivalent.grouped_barcode
        )

        if inplace:
            self._barcode = canonical_equivalent.barcode
            self._grouped_barcode = canonical_equivalent.grouped_barcode
        else:
            return BVBarcode(
                canonical_equivalent.barcode,
                grouped_barcode = canonical_equivalent.grouped_barcode
            )

    def _get_equivalent_barcodes(
        self
    ) -> tuple[tuple[tuple[int, ...], tuple[int, ...]]]:
        """
        Returns a tuple of three (`barcode`, `grouped_barcode`) pairs, each
        represented as a tuple of two tuples; the three (`barcode`,
        `grouped_barcode`) pairs are equivalent by rotation around the
        threefold symmetry axis of the bullvalene, i.e. they correspond to
        the same configurational isomer of the bullvalene.

        Returns:
            tuple[tuple[tuple[int, ...], tuple[int, ...]]]: Tuple of three
            (`barcode`, `grouped_barcode`) pairs that are equivalent by
            rotation around the threefold symmetry axis of the bullvalene.
        """
        
        return tuple(
            (rotate_tuple(
                self._barcode[:-1], i * 3
            ) + self._barcode[-1:], 
            rotate_tuple(
                self._grouped_barcode[:-1], i * 3
            ) + self._grouped_barcode[-1:])
            for i in range(3)
        )
