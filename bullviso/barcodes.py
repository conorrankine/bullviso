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
from itertools import permutations, count
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
        canonicalize: bool = False
    ):
        """
        Initialises a new bullvalene barcode (`BVBarcode`) instance.

        Args:
            barcode (tuple): A tuple of 10 integers in the range 0-9
                inclusive representing the bullvalene barcode.
            canonicalize (bool, optional): If `True`, the bullvalene barcode
                is canonicalized. Defaults to `False`.

        Raises:
            ValueError: If the bullvalene barcode i) is not of length 10, ii)
                contains non-integer elements, or iii) contains integers less
                than 0 or greater than 9.
        """
       
        if len(barcode) == 10 and all(
            isinstance(x, int) and 0 <= x <= 9 for x in barcode
        ):        
            self._barcode = barcode
        else:
            raise ValueError(
                f'{barcode} is not a valid barcode'
            )

        self._canonical_barcode = min(self._get_equivalent_barcodes())

        self._barcode_labels = self._generate_barcode_labels()

        if canonicalize:
            self.canonicalize()

    @classmethod
    def from_substituents(
        cls: type[T],
        sub_smiles: list[str],
        sub_attach_idx: list[int | list[int]],
        canonicalize: bool = False
    ) -> T:
        """
        Creates a new bullvalene barcode from a list of substituent SMILEs
        strings and a list of attachment points.

        Args:
            sub_smiles (list[str]): List of SMILEs strings specifying the
                substituents attached to the bullvalene.
            sub_attach_idx (list[Union[int, list[int]]]): List of substituent
                attachment points, supplied (for each substituent) as either:
                    - int: defining a single attachment point;
                    - list[int]: defining multiple attachment points.
            canonicalize (bool, optional): If `True`, the bullvalene barcode
                is canonicalised. Defaults to `False`.

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
            [equivalent_group_map[group] for group in groups],
            length = 10,
            direction = 'left'
        )

        return cls(barcode, canonicalize = canonicalize)

    def __str__(
        self
    ) -> str:
        """
        Returns a string representation of the bullvalene barcode.

        Returns:
            str: String representation of the bullvalene barcode.
        """
        
        return ''.join(map(str, self._barcode))

    def __hash__(
        self
    ) -> int:
        """
        Returns a hash for the bullvalene barcode; two bullvalene barcodes
        have the same hash if their canonical representations are the same.

        Returns:
            int: Hash for the bullvalene barcode.
        """
        
        return hash(self._canonical_barcode)
    
    def __eq__(
        self,
        barcode: 'BVBarcode'
    ) -> bool:
        """
        Returns the result of an equality test between the bullvalene barcode
        and another bullvalene barcode; the equality test passes if the
        canonical representations of both bullvalene barcodes are the same.

        Args:
            barcode (BVBarcode): Bullvalene barcode.

        Returns:
            bool: `True` if both bullvalene barcodes have the same canonical
                representation, else `False`.
        """

        if isinstance(barcode, type(self)):
            return (self._canonical_barcode == barcode._canonical_barcode)
        else:
            return False

    @property
    def barcode(
        self
    ) -> tuple[int, ...]:
        """
        Returns the bullvalene barcode as a tuple.

        Returns:
            tuple: The bullvalene barcode as a tuple.
        """

        return self._barcode
    
    @property
    def barcode_labels(
        self
    ) -> tuple[int, ...]:
        """
        Returns the bullvalene barcode labels as a tuple.

        Returns:
            tuple: The bullvalene barcode labels as a tuple.
        """

        return self._barcode_labels
    
    def permutations(
        self
    ) -> Generator['BVBarcode', None, None]:
        """
        Yields all possible permutations of the bullvalene barcode.

        Yields:
            BVBarcode: `BVBarcode` instances for each possible permutation of
                the bullvalene barcode.
        """

        for permutation_idx in permutations(range(10)):
            barcode = tuple(self._barcode[i] for i in permutation_idx)
            yield type(self)(barcode)

    def equivalents(
        self
    ) -> tuple['BVBarcode']:
        """
        Returns the three bullvalene barcodes equivalent by rotation around the
        threefold symmetry axis of the bullvalene.

        Returns:
            Tuple[BVBarcode]: Tuple of three `BVBarcode` instances equivalent
            by rotation around the threefold symmetry axis of the bullvalene.
        """

        return tuple(
            type(self)(barcode) for barcode in self._get_equivalent_barcodes()
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
        Canonicalizes the bullvalene barcode.
        
        A bullvalene barcode is canonical if it is the lexicographically
        smallest option of the three bullvalene barcodes equivalent by rotation
        around the threefold symmetry axis of the bullvalene.

        Note: the `_barcode_labels` attribute is refreshed on canonicalization.

        Args:
            inplace (bool): If `True`, the `_barcode` attribute is updated in
                in place, else a new `BVBarcode` instance is returned with the
                `_barcode` attribute canonicalized.

        Returns:
            BVBarcode: Canonicalized `BVBarcode` instance if `inplace = False`
                else `None` if `inplace = True`.
        """

        if inplace:
            self._barcode = self._canonical_barcode
            self._barcode_labels = self._generate_barcode_labels()
        else:
            return type(self)(
                self._canonical_barcode
            )

    def is_canonical(
        self
    ) -> bool:
        """
        Checks if the bullvalene barcode is canonical.
        
        A bullvalene barcode is canonical if it is the lexicographically
        smallest option of the three bullvalene barcodes equivalent by rotation
        around the threefold symmetry axis of the bullvalene.

        Returns:
            bool: `True` if the bullvalene barcode is canonical, else `False`.
        """

        return (
            self._barcode == self._canonical_barcode
        )
    
    def is_chiral(
        self
    ) -> bool:
        """
        Checks if the bullvalene barcode is chiral.

        A bullvalene barcode (a1,a2,a3,b1,b2,b3,c1,c2,c3,d1) is chiral if
        (a1,a2,a3) != (b1,b2,b3) != (c1,c2,c3).

        Returns:
            bool: `True` if the bullvalene barcode is chiral, else `False`.
        """

        return (
            len(set(
                self._barcode[(i * 3):(i * 3) + 3] for i in range(3)
            )) == 3
        )
    
    def is_connected(
        self
    ) -> bool:
        
        raise NotImplementedError(
            '`is_connected()` is currently a placeholder method'
        )

    def _generate_barcode_labels(
        self
    ) -> tuple[int, ...]:
        """
        Generates a tuple of unique integer labels for the non-zero elements
        of a bullvalene barcode.

        Each non-zero element in the bullvalene barcode is assigned a unique
        integer label according to its ascending rank order; zeros are
        preserved. Rank ordering is stable, i.e., earlier occurences of equal-
        valued elements receive smaller unique integer labels.

        Returns:
            tuple[int, ...]: Tuple of equal length to the bullvalene barcode
                with non-zero elements replaced by unique integer labels and
                zeros preserved.
        """
        
        nonzero_bits = [
            (i, bit) for i, bit in enumerate(self._barcode) if bit != 0
        ]

        sorted_nonzero_bits = sorted(nonzero_bits, key = lambda t: t[1])

        counter = count(1)
        label_map = {idx: next(counter) for idx, _ in sorted_nonzero_bits}

        return [
            label_map.get(i, 0) for i in range(len(self._barcode))
        ]

    def _get_equivalent_barcodes(
        self
    ) -> tuple[tuple[int, ...]]:
        """
        Generates a tuple of three bullvalene barcodes (in integer tuple
        format) equivalent by rotation around the threefold symmetry axis of
        the bullvalene.

        Returns:
            tuple[tuple[int, ...]]: Tuple of three bullvalene barcodes (in
            integer tuple format) equivalent by rotation around the threefold
            symmetry axis of the bullvalene.
        """
        
        return tuple(
            (utils.roll(self._barcode[:-1], i * 3) + self._barcode[-1:])
            for i in range(3)
        )

    def _get_connected_barcodes(
        self
    ) ->  tuple[tuple[int, ...]]:
        
        raise NotImplementedError(
            '`_get_connected_barcodes()` is currently a placeholder method'
        )

# =============================================================================
#                                     EOF
# =============================================================================