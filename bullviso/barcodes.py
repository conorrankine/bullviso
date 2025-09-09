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

from . import utils
from itertools import permutations
from typing import Generator, TypeVar

# =============================================================================
#                                   GLOBALS
# =============================================================================

T = TypeVar('T', bound = 'BVBarcode')

# =============================================================================
#                                   CLASSES
# =============================================================================

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

        self._canonicalized_barcode, self._canonicalized_grouped_barcode = (
            min(self._get_equivalent_barcodes(), key = lambda t: t[1])
        )

        if canonicalize:
            self.canonicalize()

    @classmethod
    def from_substituents(
        cls: type[T],
        sub_smiles: list[str],
        sub_attach_idx: list[int | list[int]],
        canonicalize: bool = True
    ) -> T:
        """
        Creates a new bullvalene barcode from a list of substituent SMILEs
        strings and a list of attachment ooints.

        Args:
            sub_smiles (list[str]): List of SMILEs strings specifying the
                substituents attached to the bullvalene.
            sub_attach_idx (list[Union[int, list[int]]]): List of substituent
                attachment points, supplied (for each substituent) as either:
                    - int: defining a single attachment point;
                    - list[int]: defining multiple attachment points.
            canonicalize (bool, optional): If True, the bullvalene barcode is
                canonicalised. Defaults to True.

        Raises:
            ValueError: If `sub_smiles` and `sub_attach_idx` are not of equal
                length.
            ValueError: If the number of attachment points is greater than 9.

        Returns:
            BVBarcode: Bullvalene barcode.
        """
        
        if not utils.all_same_length(sub_smiles, sub_attach_idx):
            raise ValueError(
                f'`sub_smiles` and `sub_attach_idx` should have equal '
                f'length; got lists with length {len(sub_smiles)} and '
                f'{len(sub_attach_idx)}'
            )

        groups = [
            f'{sub_smiles[n]}_{idx}' for (n, idx)
                in utils.iterate_and_index(sub_attach_idx)
        ]

        n_attachment_points = len(groups)
        if n_attachment_points > 9:
            raise ValueError(
                f'too many attachment points defined: got '
                f'{n_attachment_points} (maximum allowed = 9)'
            )

        equivalent_group_map = {
            group: i for i, group in enumerate(
                utils.unique_elements(groups), start = 1
            )
        }

        barcode = utils.pad_list(
            [i for i in range(1, len(groups) + 1)],
            length = 10,
            direction = 'left'
        )

        grouped_barcode = utils.pad_list(
            [equivalent_group_map[group] for group in groups],
            length = 10,
            direction = 'left'
        )

        return cls(
            barcode,
            grouped_barcode = grouped_barcode,
            canonicalize = canonicalize
        )

    def __str__(
        self
    ) -> str:
        """
        Returns a string representation for the `_grouped_barcode` attribute
        of the `BVBarcode` instance.

        Returns:
            str: String representation for the `_grouped_barcode` attribute.
        """
        
        return ''.join(map(str, self._grouped_barcode))

    def __hash__(
        self
    ) -> int:
        """
        Returns a hash for the `BVBarcode` instance; two `BVBarcode` instances
        will have the same hash if their canonical representations (i.e. their
        `_canonicalized_grouped_barcode` attributes) are the same.

        Returns:
            int: Hash for the `BVBarcode` instance.
        """
        
        return hash(
            self._canonicalized_grouped_barcode
        )
    
    def __eq__(
        self,
        barcode: 'BVBarcode'
    ) -> bool:
        """
        Returns the result of an equality test between the `BVBarcode` instance
        and another `BVBarcode` instance (`barcode`); two `BVBarcode` instances
        are considered equal if their canonical representations (i.e. their
        `_canonicalized_grouped_barcode` attributes) are the same.

        Args:
            barcode (BVBarcode): `BVBarcode` instance to test for equality.

        Returns:
            bool: True if both `BVBarcode` instances have the same
            canonicalized grouped barcode representation, else False.
        """

        if isinstance(barcode, type(self)):
            return (
                self._canonicalized_grouped_barcode
                == barcode._canonicalized_grouped_barcode
            )
        else:
            return False

    @property
    def barcode(
        self
    ) -> tuple[int, ...]:
        """
        Returns the isomer barcode as a tuple.

        Returns:
            tuple: The isomer barcode as a tuple.
        """

        return self._barcode
    
    @property
    def grouped_barcode(
        self
    ) -> tuple[int, ...]:
        """
        Returns the grouped isomer barcode as a tuple.

        Returns:
            tuple: The grouped isomer barcode as a tuple.
        """

        return self._grouped_barcode
    
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

        for permutation_idx in permutations(range(len(self._barcode))):
            barcode = tuple(
                self._barcode[i] for i in permutation_idx
            )
            grouped_barcode = tuple(
                self._grouped_barcode[i] for i in permutation_idx
            )
            yield type(self)(
                barcode,
                grouped_barcode = grouped_barcode
            )

    def equivalents(
        self
    ) -> tuple['BVBarcode']:
        """
        Returns a tuple of three `BVBarcode` instances that are equivalent by
        rotation around the threefold symmetry axis of the bullvalene, i.e.
        they correspond to the same configurational isomer.

        Returns:
            Tuple[BVBarcode]: Tuple of three `BVBarcode` instances that are
            equivalent by rotation around the threefold symmetry axis of the
            bullvalene.
        """

        return tuple(
            type(self)(barcode, grouped_barcode = grouped_barcode)
            for barcode, grouped_barcode in self._get_equivalent_barcodes()
        )
    
    def connections(
        self
    ) -> tuple['BVBarcode']:
        
        raise NotImplementedError(
            '`connections()` is currently a placeholder method'
        )
    
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

        if inplace:
            self._barcode = self._canonicalized_barcode
            self._grouped_barcode = self._canonicalized_grouped_barcode
        else:
            return type(self)(
                self._canonicalized_barcode,
                grouped_barcode = self._canonicalized_grouped_barcode
            )

    def is_canonicalized(
        self
    ) -> bool:
        """
        Checks if the bullvalene isomer barcode is canonicalized, i.e. if
        the `_grouped_barcode` and `_canonicalized_grouped_barcode` attributes
        of the `BVBarcode` instance are the same; if canonicalized, the
        bullvalene isomer barcode corresponds to the smallest of the set of
        equivalent barcodes lexicographically.

        Returns:
            bool: `True` if the bullvalene isomer barcode is canonicalized,
                else `False`.
        """

        return (
            self._grouped_barcode == self._canonicalized_grouped_barcode
        )
    
    def is_chiral(
        self
    ) -> bool:
        """
        Checks if the bullvalene isomer barcode corresponds to a chiral
        bullvalene, i.e. (a1,a2,a3) != (b1,b2,b3) != (c1,c2,c3) for the
        bullvalene isomer barcode (a1,a2,a3,b1,b2,b3,c1,c2,c3,d1).

        Returns:
            bool: `True` if the bullvalene isomer barcode corresponds to a
                chiral bullvalene, else `False`.
        """

        return (
            len(set(
                self._grouped_barcode[(i * 3):(i * 3) + 3] for i in range(3)
            )) == 3
        )
    
    def is_connected(
        self
    ) -> bool:
        
        raise NotImplementedError(
            '`is_connected()` is currently a placeholder method'
        )

    def _get_equivalent_barcodes(
        self
    ) -> tuple[tuple[tuple[int, ...], tuple[int, ...]]]:
        """
        Returns a tuple of three (`barcode`, `grouped_barcode`) pairs, each
        represented as a tuple of two tuples; the three (`barcode`,
        `grouped_barcode`) pairs are equivalent by rotation around the
        threefold symmetry axis of the bullvalene, i.e. they correspond to
        the same configurational isomer.

        Returns:
            tuple[tuple[tuple[int, ...], tuple[int, ...]]]: Tuple of three
            (`barcode`, `grouped_barcode`) pairs that are equivalent by
            rotation around the threefold symmetry axis of the bullvalene.
        """
        
        return tuple(
            (utils.roll(
                self._barcode[:-1], i * 3
            ) + self._barcode[-1:], 
            utils.roll(
                self._grouped_barcode[:-1], i * 3
            ) + self._grouped_barcode[-1:])
            for i in range(3)
        )

    def _get_connected_barcodes(
        self
    ) ->  tuple[tuple[tuple[int, ...], tuple[int, ...]]]:
        
        raise NotImplementedError(
            '`_get_connected_barcodes()` is currently a placeholder method'
        )

# =============================================================================
#                                     EOF
# =============================================================================