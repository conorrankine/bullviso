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
from enum import Enum
from typing import Iterable, Iterator, Sequence
from rdkit import Chem
from .utils.list_utils import all_same_length, pad_list
from .barcodes import BVBarcode, BVTSBarcode

__all__ = [
    "Substituent",
    "Substituents",
    "Bullvalene"
]

# =============================================================================
#                                  CONSTANTS
# =============================================================================

ATOM_STEREO_FLAGS = Enum(
    'ATOM_STEREO_FLAGS', {
        'CW': Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
        'CCW': Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
    }
)

BOND_STEREO_FLAGS = Enum(
    'BOND_STEREO_FLAGS', {
        'CIS': Chem.rdchem.BondStereo.STEREOCIS,
        'TRANS': Chem.rdchem.BondStereo.STEREOTRANS
    }
)

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

    def __init__(
        self,
        smiles: str,
        attach_idx: int | Iterable[int]
    ) -> None:

        self.smiles = smiles
        
        self.mol = Chem.MolFromSmiles(self.smiles)
        if self.mol is None:
            raise ValueError(
                f'error generating a valid substituent from SMILES string '
                f'\"{self.smiles}\"'
            )

        if isinstance(attach_idx, int):
            attach_idx = [attach_idx]
        attach_idx = tuple(attach_idx)
        if not attach_idx:
            raise ValueError(
                f'empty `attach_idx` for \"{self.smiles}\"; no attachment '
                f'indices are defined for this substituent'
            )
        if len(attach_idx) > 9:
            raise ValueError(
                f'>9 elements in `attach_idx` for \"{self.smiles}\"; too many '
                f'attachment indices are defined for this substituent'
            )
        if min(attach_idx) < 0 or max(attach_idx) >= self.mol.GetNumAtoms():
            raise ValueError(
                f'invalid attachment index defined for \"{self.smiles}\" in '
                f'{attach_idx}; valid attachment indices for this substituent '
                f'are in the range 0 <= `idx` < {self.mol.GetNumAtoms()}'
            )
        self.attach_idx = attach_idx

        canonical_ranks = Chem.CanonicalRankAtoms(self.mol, breakTies = False)
        self.attach_sym_cls = tuple(
            canonical_ranks[idx] for idx in self.attach_idx
        )

        self.is_multidentate = len(self.attach_idx) > 1

class Substituents():

    def __init__(
        self,
        substituents: Iterable[Substituent] | None = None
    ) -> None:

        self._substituents: list[Substituent] = []
        self._n_substituents: int = 0
        self._n_attach_idx: int = 0
        if substituents is not None:
            for substituent in substituents:
                self.add(substituent)

    def add(
        self,
        substituent: Substituent
    ) -> None:

        if not isinstance(substituent, Substituent):
            raise ValueError(
                f'expected a `Substituent` instance; got an object of type '
                f'{type(substituent).__name__}'
            )
        n_attach_idx_proj = self._n_attach_idx + len(substituent.attach_idx)
        if n_attach_idx_proj > 9:
            raise ValueError(
                f'support is only available for up to 9 attachment points; '
                f'adding this substituent gives a projected number of '
                f'attachment points of {n_attach_idx_proj}'
            )
        self._substituents.append(substituent)

        self._n_substituents += 1
        self._n_attach_idx += len(substituent.attach_idx)

    def to_barcode(
        self,
        canonicalize: bool = True,
        transition_state: bool = False
    ) -> BVBarcode | BVTSBarcode:

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

    def _get_hash_strs(
        self
    ) -> list[str]:

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

        return iter(self._substituents)

    def __len__(
        self
    ) -> int:

        return len(self._substituents)

    def __getitem__(
        self,
        idx: int
    ) -> Substituent:

        return self._substituents[idx]

class Bullvalene():

    def __init__(
        self,
        substituents: Substituents,
        transition_state: bool = False
    ) -> None:

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
            _remove_explicit_hydrogens(
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
            _set_atom_stereochemistry(substituted_bullvalene, stereo_map)

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
            BVBarcode | BVTSBarcode: Copy of the bullvalene barcode.
        """

        return self._barcode

    @property
    def template_mol(
        self
    ) -> Chem.Mol:
        """
        Returns:
            Chem.Mol: Copy of the bullvalene template geometry.
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
            Chem.Mol: Bullvalene template geometry.
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
    Generates a `Substituents` instance from a user-friendly shorthand multi-
    substituent specification; ensures that a unique bullvalene barcode bit
    label is assigned to each attachment index.

    Args:
        substituent_smiles (Sequence[str]): Sequence of SMILES strings (one
            per substituent).
        substituent_counts (Sequence[int]): Sequence of integers (one per
            substituent) defining the substituent count. 
        attach_idx (Sequence[Sequence[int]]): Sequence of integer sequences
            (one per substituent) defining the atom indices of the attachment
            point(s) for the corresponding substituent.

    Raises:
        ValueError: If `substituent_smiles`, `substituent_count`, and
            `attach_idx` are not all of equal length.
        ValueError: If any element of `substituent_count` is <= 0.

    Returns:
        Substituents: `Substituents` instance with a unique bullvalene barcode
            bit label assigned to each attachment index.
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

def _remove_explicit_hydrogens(
    mol: Chem.Mol,
    atom_idx: int
) -> None:
    """
    Removes all explicit hydrogen atoms from an atom in a molecule. Removal is
    an in-place operation, i.e., the molecule is directly modified.

    Args:
        mol (Chem.Mol): Molecule.
        atom_idx (int): Atom index.

    Raises:
        IndexError: If the atom index is out of bounds, i.e., less than zero
            or greater than the number of atoms in the molecule.
    """

    if not (0 <= atom_idx < mol.GetNumAtoms()):
        raise IndexError(
            f'atom index {atom_idx} is out of bounds for a molecule with '
            f'{mol.GetNumAtoms()} atoms'
        )

    atom = mol.GetAtomWithIdx(atom_idx)
    
    atom.SetNumExplicitHs(0)
    atom.UpdatePropertyCache()

def _set_atom_stereochemistry(
    mol: Chem.Mol,
    stereo_map: dict[int, str]
) -> None:
    """
    Sets stereochemical tags of atoms in a molecule according to a dictionary
    mapping of atom indices to stereochemical (chirality) flags.

    Stereochemical flags:
        - 'cw' = clockwise (tetrahedral)
        - 'ccw' = counterclockwise (tetrahedral)

    Args:
        mol (Chem.Mol): Molecule.
        stereo_map (dict): Dictionary mapping of atomic indices to
            stereochemical (chirality) flags.

    Raises:
        IndexError: If any atom index in `stereo_map` is out of bounds, i.e.,
            less than zero or greater than the number of atoms in the molecule.
        ValueError: If any stereochemical (chirality) flag is not recognised,
            i.e., if it is not one of {'cw', 'ccw'}.
    """
    
    for atom_idx, stereo_flag in stereo_map.items():
        if not (0 <= atom_idx < mol.GetNumAtoms()):
            raise IndexError(
                f'atom index {atom_idx} is out of bounds for a molecule with '
                f'{mol.GetNumAtoms()} atoms'
            )
        atom = mol.GetAtomWithIdx(atom_idx)
        if stereo_flag.upper() not in ATOM_STEREO_FLAGS.__members__:
            raise ValueError(
                f'invalid stereochemical flag: got {stereo_flag}'
            )
        atom.SetChiralTag(ATOM_STEREO_FLAGS[stereo_flag].value)

def _set_bond_stereochemistry(
    mol: Chem.Mol,
    stereo_map: dict[tuple[int, int], str]
) -> None:
    """
    Sets stereochemical tags of bonds in a molecule according to a dictionary
    mapping of bonds (pairs of atom indices) to stereochemical flags.

    Args:
        mol ((Chem.Mol)): Molecule.
        stereo_map (dict[tuple[int, int], str]): Dictionary mapping of bonds
            (pairs of atom indices) to stereochemical flags.


    Raises:
        IndexError: If any atom index defining part of a bond in `stereo_map`
            is out of bounds, i.e., less than zero or greater than the number
            of atoms in the molecule.
        ValueError: If any stereochemical flag is not recognised, i.e., if it
            is not one of {'cis', 'trans'}.
    """
    
    for atom_idxs, stereo_flag in stereo_map.items():
        for atom_idx in atom_idxs:
            if not (0 <= atom_idx < mol.GetNumAtoms()):
                raise IndexError(
                    f'atom index {atom_idx} is out of bounds for a molecule '
                    f'with {mol.GetNumAtoms()} atoms'
                )
        bond = mol.GetBondBetweenAtoms(*atom_idxs)
        if stereo_flag.upper() not in BOND_STEREO_FLAGS.__members__:
            raise ValueError(
                f'invalid stereochemical flag: got {stereo_flag}'
            )
        bond.SetStereo(BOND_STEREO_FLAGS[stereo_flag].value)

# =============================================================================
#                                     EOF
# =============================================================================
