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

from __future__ import annotations
from importlib import resources
from itertools import count
from typing import Iterable, Iterator, Sequence
from rdkit import Chem
from .utils.list_utils import (
    all_same_length,
    pad_list
)
from .utils.rdkit_utils import (
    remove_explicit_hydrogens,
    set_atom_stereochemistry,
    set_bond_stereochemistry
)
from .barcodes import BVBarcode, BVTSBarcode

__all__ = [
    "Substituent",
    "Substituents",
    "Bullvalene"
]

# =============================================================================
#                                  CONSTANTS
# =============================================================================

BULLVALENE_STEREO_MAP = {
    0: 'CCW', 3: 'CW', 6: 'CCW', 9: 'CW'
}
BULLVALENE_TS_STEREO_MAP = {
    0: 'CW', 2: 'CCW', 3: 'CCW', 5: 'CW', 6: 'CCW', 9: 'CW'
}

BULLVALENE_STRUCTURE_PACKAGE = 'bullviso.structures'
BULLVALENE_MOL_FILE = 'bv.sdf'
BULLVALENE_TS_MOL_FILE = 'bv_ts.sdf'

# =============================================================================
#                                   CLASSES
# =============================================================================

class Substituent():

    __slots__ = (
        '_smiles',
        '_mol',
        '_attach_idx',
        '_attach_sym_cls',
        '_is_multidentate'
    )

    def __init__(
        self,
        smiles: str,
        attach_idx: int | Iterable[int]
    ) -> None:
        """
        Initialises a Substituent instance from a SMILES string and one or
        more atom indices defining the attachment point(s).

        Args:
            smiles (str): SMILES string.
            attach_idx (int | Iterable[int]): Atom index, or iterable of atom
                indices, specifying the attachment point(s) on the substituent.

        Raises:
            ValueError: If `smiles` does not define a valid substituent.
            ValueError: If no attachment indices are supplied.
            ValueError: If more than nine attachment indices are supplied.
            ValueError: If any attachment index is out of bounds for the
                substituent, i.e., less than zero, or greater than / equal to
                the maximum number of atoms in the substituent.
        """

        self._smiles = smiles
        
        mol = Chem.MolFromSmiles(self._smiles)
        if mol is None:
            raise ValueError(
                f'error generating a valid substituent from SMILES string '
                f'\"{self._smiles}\"'
            )

        if isinstance(attach_idx, int):
            attach_idx = [attach_idx]
        attach_idx = tuple(attach_idx)
        if not attach_idx:
            raise ValueError(
                f'empty `attach_idx` for \"{self._smiles}\"; no attachment '
                f'indices are defined for this substituent'
            )
        if len(attach_idx) > 9:
            raise ValueError(
                f'>9 elements in `attach_idx` for \"{self._smiles}\"; too many '
                f'attachment indices are defined for this substituent'
            )
        if min(attach_idx) < 0 or max(attach_idx) >= mol.GetNumAtoms():
            raise ValueError(
                f'invalid attachment index defined for \"{self._smiles}\" in '
                f'{attach_idx}; valid attachment indices for this substituent '
                f'are in the range 0 <= `idx` < {mol.GetNumAtoms()}'
            )
        self._mol = mol
        self._attach_idx = attach_idx

        canonical_ranks = Chem.CanonicalRankAtoms(self._mol, breakTies = False)
        self._attach_sym_cls = tuple(
            canonical_ranks[idx] for idx in self._attach_idx
        )

        self._is_multidentate = len(self._attach_idx) > 1

    @property
    def smiles(
        self
    ) -> str:
        """
        Returns:
            str: SMILES string representing the substituent.
        """

        return self._smiles

    @property
    def mol(
        self
    ) -> Chem.Mol:
        """
        Returns:
            Chem.Mol: Copy of the RDKit Chem.Mol instance representing the
                substituent.
        """

        return Chem.Mol(self._mol)

    @property
    def attach_idx(
        self
    ) -> tuple[int, ...]:
        """
        Returns:
            tuple[int, ...]: Attachment indices for the substituent.
        """

        return self._attach_idx

    @property
    def attach_sym_cls(
        self
    ) -> tuple[int, ...]:
        """
        Returns:
            tuple[int, ...]: Symmetry classes of the atoms indexed by the
                attachment indices for the substituent.
        """

        return self._attach_sym_cls

    @property
    def is_multidentate(
        self
    ) -> bool:
        """
        Returns:
            bool: `True` if the substituent is multidentate, i.e., if it has
                multiple attachment points, else `False`.
        """

        return self._is_multidentate

class Substituents():

    __slots__ = (
        '_substituents',
    )

    def __init__(
        self,
        substituents: Iterable[Substituent] | None = None
    ) -> None:
        """
        Initialises a Substituents instance from an iterable of Substituent
        instances.

        Args:
            substituents (Iterable[Substituent] | None, optional): Iterable of
                Substituent instances to add to the Substituents container. If
                `None`, an empty Substituents container is initialised.
                Defaults to `None`.
        """

        self._substituents: list[Substituent] = []
        if substituents is not None:
            for substituent in substituents:
                self.add(substituent)

    def add(
        self,
        substituent: Substituent
    ) -> None:
        """
        Adds a substituent to the Substituents container.

        Args:
            substituent (Substituent): Substituent to add to the Substituents
                container.

        Raises:
            TypeError: If `substituent` is not a Substituent instance.
            ValueError: If adding `substituent` to the Substituents container
                would result in more than nine total attachment points.
        """

        if not isinstance(substituent, Substituent):
            raise TypeError(
                f'expected a `Substituent` instance; got an object of type '
                f'{type(substituent).__name__}'
            )
        n_attach_idx_proj = self.n_attach_idx + len(substituent.attach_idx)
        if n_attach_idx_proj > 9:
            raise ValueError(
                f'support is only available for up to 9 attachment points; '
                f'adding this substituent gives a projected number of '
                f'attachment points of {n_attach_idx_proj}'
            )
        self._substituents.append(substituent)

    def to_barcode(
        self,
        canonicalize: bool = True,
        transition_state: bool = False
    ) -> BVBarcode | BVTSBarcode:
        """
        Returns the bullvalene barcode for the Substituents container.

        Args:
            canonicalize (bool, optional): If `True`, the bullvalene barcode is
                canonicalised. Defaults to `True`.
            transition_state (bool, optional): If `True`, a transition-state
                bullvalene barcode is returned. Defaults to `False`.

        Returns:
            BVBarcode | BVTSBarcode: Bullvalene barcode.
        """

        hash_strs = self._get_hash_strs()

        unique_hash_strs = list(dict.fromkeys(hash_strs))

        hash_str_to_group_map = {
            hash_str: i + 1 for i, hash_str in enumerate(unique_hash_strs)
        }

        barcode = pad_list(
            [hash_str_to_group_map[hash_str] for hash_str in hash_strs],
            length = 10,
            direction = 'left'
        )

        if transition_state:
            return BVTSBarcode(barcode, canonicalize = canonicalize)
        else:
            return BVBarcode(barcode, canonicalize = canonicalize)

    @property
    def n_substituents(
        self
    ) -> int:
        """
        Returns:
            int: Number of substituents in the Substituents container.
        """

        return len(self._substituents)

    @property
    def n_attach_idx(
        self
    ) -> int:
        """
        Returns:
            int: Total number of attachment points across all substituents
                in the Substituents container.
        """

        return sum(len(s.attach_idx) for s in self._substituents)

    def _get_hash_strs(
        self
    ) -> list[str]:
        """
        Generates hash strings used to map substituent attachment points to
        bullvalene barcode bit labels.

        Each hash string encodes the substituent SMILES string together with
        the symmetry class of the corresponding attachment atom. For
        multidentate substituents, the substituent index is also included so
        that attachment points belonging to different instances of the same
        multidentate substituent are assigned distinct barcode bit labels.

        Returns:
            list[str]: Hash strings for each substituent attachment point.
        """

        hash_strs: list[str] = []
        for substituent_idx, substituent in enumerate(self._substituents):
            for attach_sym_cls in substituent.attach_sym_cls:
                if substituent.is_multidentate:
                    hash_str = '{}_{}_{}'.format(
                        substituent.smiles,
                        attach_sym_cls,
                        substituent_idx
                    )
                else:
                    hash_str = '{}_{}'.format(
                        substituent.smiles,
                        attach_sym_cls
                    )
                hash_strs.append(hash_str)

        return hash_strs

    def __iter__(
        self
    ) -> Iterator[Substituent]:
        """
        Returns:
            Iterator[Substituent]: Iterator over the substituents in the
                Substituents container.
        """

        return iter(self._substituents)

    def __len__(
        self
    ) -> int:
        """
        Returns:
            int: Number of substituents in the Substituents container.
        """

        return len(self._substituents)

    def __getitem__(
        self,
        idx: int
    ) -> Substituent:
        """
        Args:
            idx (int): Index.

        Returns:
            Substituent: Substituent at the specified index in the Substituents
                container.
        """

        return self._substituents[idx]

class Bullvalene():

    def __init__(
        self,
        substituents: Substituents,
        transition_state: bool = False
    ) -> None:
        """
        Initialises a Bullvalene instance from a Substituents instance
        specifying the attached substituents and their configuration; templates
        either a regular or Cope rearrangement transition-state bullvalene.

        Args:
            substituents (Substituents): Substituents.
            transition_state (bool, optional): If `True`, a Cope rearrangement
                transition-state template is initialised. Defaults to `False`.

        Raises:
            ValueError: If `substituents` is not a `Substituents` instance.
        """

        if not isinstance(substituents, Substituents):
            raise ValueError(
                f'expected a `Substituents` instance; got an object of type '
                f'{type(substituents).__name__}'
            )

        self.substituents = substituents
        self.transition_state = transition_state

        self._barcode = substituents.to_barcode(
            canonicalize = True,
            transition_state = transition_state
        )
        
        self._template = self._load_template(
            transition_state = transition_state
        )
        
        self._barcode_bit_to_sub_atom_offset = (
            self._get_barcode_bit_to_sub_atom_offset()
        )

    def build(
        self,
        barcode: BVBarcode | BVTSBarcode,
        sanitize: bool = True,
        set_properties: bool = True,
        set_stereochemistry: bool = True
    ) -> Chem.Mol:
        """
        Builds a substituted bullvalene structure corresponding to the supplied
        bullvalene barcode by connecting the substituents to the bullvalene
        template.

        Args:
            barcode (BVBarcode | BVTSBarcode): Bullvalene barcode.
            sanitize (bool, optional): If `True`, the substituted bullvalene is
                sanitised before being returned. Defaults to `True`.
            set_properties (bool, optional): If `True`, the `barcode` and
                `transition_state` properties are set for the substituted
                bullvalene. Defaults to `True`.
            set_stereochemistry (bool, optional): If `True`, the substituted
                bullvalene is assigned stereochemistry. Defaults to `True`.

        Raises:
            ValueError: If `barcode` is incompatible with the stored
                substituents.

        Returns:
            Chem.Mol: Substituted bullvalene.
        """

        self._validate_barcode_compatibility(barcode)

        substituted_bullvalene = Chem.RWMol(Chem.Mol(self._template))
        for substituent in self.substituents:
            substituted_bullvalene.InsertMol(Chem.Mol(substituent.mol))

        for attach_idx_bullvalene, barcode_bit in enumerate(barcode.barcode_labels):
            if barcode_bit == 0:
                continue
            attach_idx_substituent = (
                self._template.GetNumAtoms()
                + self._barcode_bit_to_sub_atom_offset[barcode_bit]
            )
            remove_explicit_hydrogens(
                substituted_bullvalene, attach_idx_bullvalene
            )
            substituted_bullvalene.AddBond(
                attach_idx_bullvalene,
                attach_idx_substituent,
                Chem.rdchem.BondType.SINGLE
            )
        substituted_bullvalene = substituted_bullvalene.GetMol()

        if sanitize:
            Chem.SanitizeMol(substituted_bullvalene)

        if set_stereochemistry:
            stereo_map = (
                BULLVALENE_TS_STEREO_MAP if self.transition_state
                else BULLVALENE_STEREO_MAP
            )
            set_atom_stereochemistry(substituted_bullvalene, stereo_map)

        if set_properties:
            substituted_bullvalene.SetProp('barcode', str(barcode))
            substituted_bullvalene.SetBoolProp(
                'transition_state', self.transition_state
            )

        return substituted_bullvalene

    @property
    def barcode(
        self
    ) -> BVBarcode | BVTSBarcode:
        """
        Returns:
            BVBarcode | BVTSBarcode: Bullvalene barcode.
        """

        return self._barcode

    @property
    def template_mol(
        self
    ) -> Chem.Mol:
        """
        Returns:
            Chem.Mol: Copy of the bullvalene template.
        """

        return Chem.Mol(self._template)

    def _get_barcode_bit_to_sub_atom_offset(
        self
    ) -> dict[int, int]:
        """
        Precomputes the mapping of barcode bit to substituent atom offset.

        Returns:
            dict[int, int]: Mapping of barcode bit value to substituent atom
                index offset within the combined inserted substituent block.
        """

        barcode_bit_to_sub_atom_offset: dict[int, int] = {}

        barcode_bits = iter(
            bit for bit in self._barcode.barcode_labels if bit != 0
        )

        atom_offset = 0
        for substituent in self.substituents:
            for attach_idx in substituent.attach_idx:
                barcode_bit_to_sub_atom_offset[next(barcode_bits)] = (
                    atom_offset + attach_idx
                )
            atom_offset += substituent.mol.GetNumAtoms()

        return barcode_bit_to_sub_atom_offset

    def _validate_barcode_compatibility(
        self,
        barcode: BVBarcode | BVTSBarcode
    ) -> None:
        """
        Validates that a barcode is compatible with this bullvalene instance.

        Args:
            barcode (BVBarcode | BVTSBarcode): Barcode to validate.

        Raises:
            ValueError: If the barcode is incompatible with this bullvalene
                instance's substituent set.
        """

        if isinstance(barcode, BVTSBarcode) != self.transition_state:
            raise ValueError(
                f'barcode ({barcode}) is incompatible with this bullvalene '
                f'instance'
            )

        expected_barcode = self._barcode
        expected_nonzero_bits = {
            bit for bit in expected_barcode.barcode_labels if bit != 0
        }
        actual_nonzero_bits = {
            bit for bit in barcode.barcode_labels if bit != 0
        }
        if actual_nonzero_bits != expected_nonzero_bits:
            raise ValueError(
                f'barcode ({barcode}) is incompatible with this bullvalene '
                f'instance'
            )

    @staticmethod
    def _load_template(
        transition_state: bool = False
    ) -> Chem.Mol:
        """
        Loads a bullvalene or, optionally, Cope rearrangement transition state,
        template geometry as an RDKit `Chem.Mol` instance from the
        `bullviso.structures` library with preset 3D/Cartesian coordinates.

        Args:
            transition_state (bool, optional): If `True`, a Cope rearrangement
                transition state template geometry is returned. Defaults to
                `False`.

        Raises:
            ValueError: If the RDKit `Chem.Mol` instance fails to load from the
                relevant mol (.sdf) file in the `bullviso.structures` library.

        Returns:
            Chem.Mol: Bullvalene template.
        """

        template_file_name = (
            BULLVALENE_TS_MOL_FILE if transition_state else BULLVALENE_MOL_FILE
        )
        
        with resources.open_text(
            BULLVALENE_STRUCTURE_PACKAGE, template_file_name
        ) as template_file:
            template = Chem.MolFromMolBlock(template_file.read())
        if template is None:
            raise ValueError(
                f'error loading template geometry from `{template_file_name}` '
                f'in `{BULLVALENE_STRUCTURE_PACKAGE}`: check that the file '
                f'exists and points to a valid mol (.sdf) file'
            )

        return template

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

def substituents_from_specifications(
    substituent_smiles: Sequence[str],
    substituent_counts: Sequence[int],
    attach_idx: Sequence[Sequence[int]]
) -> Substituents:
    """
    Generates a Substituents instance from a user-friendly shorthand multi-
    substituent specification.

    Each element of `substituent_smiles`, `substituent_counts`, and
    `attach_idx` defines one substituent specification: the substituent SMILES
    string, the substituent count, and the atom index/indices defining the
    substituent attachment point(s), respectively.

    Args:
        substituent_smiles (Sequence[str]): Sequence of SMILES strings (one
            per substituent) specifying the substituent structure.
        substituent_counts (Sequence[int]): Sequence of integers (one per
            substituent) specifying the substituent count. 
        attach_idx (Sequence[Sequence[int]]): Sequence of integer sequences
            (one per substituent) defining the atom index/indices of the
            substituent attachment point(s) for the corresponding substituent.

    Raises:
        ValueError: If `substituent_smiles`, `substituent_counts`, and
            `attach_idx` are not all of equal length.
        ValueError: If any element of `substituent_counts` is <= 0.

    Returns:
        Substituents: Substituents instance generated from the supplied
            shorthand multi-substituent specification(s).
    """

    if not all_same_length(
        substituent_smiles, substituent_counts, attach_idx
    ):
        raise ValueError(
            f'`substituent_smiles` ({len(substituent_smiles)} elements), '
            f'`substituent_counts` ({len(substituent_counts)} elements), and '
            f'`attach_idx` ({len(attach_idx)} elements) should all be of '
            f'equal length'
        )

    for smiles_, substituent_count_ in zip(
        substituent_smiles, substituent_counts
    ):
        if substituent_count_ <= 0:
            raise ValueError(
                f'`substituent_count` should be a positive (>0) integer: got '
                f'{substituent_count_} for \"{smiles_}\"'
            )
    
    substituents = Substituents()

    for smiles_, substituent_count_, attach_idx_ in zip(
        substituent_smiles, substituent_counts, attach_idx
    ):
        for _ in range(substituent_count_):
            substituent = Substituent(
                smiles = smiles_,
                attach_idx = attach_idx_
            )
            substituents.add(substituent)

    return substituents

# =============================================================================
#                                     EOF
# =============================================================================
