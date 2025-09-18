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
from dataclasses import dataclass
from importlib import resources
from functools import lru_cache
from itertools import count
from enum import Enum
from rdkit import Chem
from .utils import all_same_length
from .barcodes import BVBarcode, BVTSBarcode

__all__ = [
    "Substituent",
    "build_bullvalene_from_barcode",
    "load_bullvalene_from_library"
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

# =============================================================================
#                                   CLASSES
# =============================================================================

@dataclass
class Substituent:
    """
    Defines a substituent that can be attached to a bullvalene.

    Attributes: 
        smiles (str): SMILES string for the substituent.
        attachment_idxs (list[int]): [#TODO: DESCRIPTION]
        barcode_bits (list[int]): [#TODO: DESCRIPTION]   
    """

    smiles: str
    attachment_idxs: list[int]
    barcode_bits: list[int]

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

def substituents_from_specifications(
    substituent_smiles: list[str],
    substituent_counts: list[int],
    attachment_idxs: list[list[int]]
) -> tuple[Substituent, ...]:
    """
    Generates a tuple of `Substituent` instances from a user-friendly shorthand
    substituent specification.

    This function automates the boilerplate associated with creating multiple
    `Substituent` instances manually and ensures that a unique bullvalene
    barcode bit value is assigned to each attachment point.

    Args:
        substituent_smiles (list[str]): Sequence of SMILES strings (one per
            substituent).
        substituent_counts (list[int]): Sequence of integers (one per
            substituent) defining the substituent count, i.e., the number of
            times that the corresponding substituent should appear. 
        attachment_idxs (list[list[int]]): Sequence of integer sequences (one
            per substituent) defining the atom indices of the attachment
            point(s) for the corresponding substituent.

    Raises:
        ValueError: If `substituent_smiles`, `substituent_count`, and
            `attachment_idxs` are not all of equal length.
        ValueError: If any element of `substituent_count` is < 0.
        ValueError: If any element of `attachment_idxs` is empty.
        ValueError: If more than 9 attachment points are specified (summed over
            all substituents).

    Returns:
        tuple[Substituent, ...]: Tuple of `Substituent` instances with unique
            bullvalene barcode bit values assigned to each attachment point.
    """

    if not all_same_length(
        substituent_smiles, substituent_counts, attachment_idxs
    ):
        raise ValueError(
            f'`substituent_smiles`, `substituent_count`, and `attachment_idxs` '
            f'should all be of equal length; got sequences with lengths '
            f'{len(substituent_smiles)}, {len(substituent_counts)}, and '
            f'{len(attachment_idxs)} (respectively)'
        )
    
    substituents: list[Substituent] = []

    barcode_bit_value = count(start = 1)

    for smiles_, substituent_count_, attachment_idxs_ in zip(
        substituent_smiles, substituent_counts, attachment_idxs
    ):
        if substituent_count_ <= 0:
            raise ValueError(
                f'`substituent_count` should be positive: got '
                f'{substituent_count_} for SMILES = \'{smiles_}\''
            )
        if not attachment_idxs_:
            raise ValueError(
                f'`attachment_idxs` should not be empty: got '
                f'{attachment_idxs_} for SMILES = \'{smiles_}\''
            )
        for _ in range(substituent_count_):
            substituents.append(
                Substituent(
                    smiles = smiles_,
                    attachment_idxs = attachment_idxs_,
                    barcode_bits = [
                        next(barcode_bit_value) for _ in attachment_idxs_
                    ]
                )
            )

    num_attachment_idxs = next(barcode_bit_value) - 1
    if num_attachment_idxs > 9:
        raise ValueError(
            f'support is only available for up to 9 attachment points: got '
            f'{num_attachment_idxs} attachment points across all substituents'
        )

    return tuple(substituents)

def build_bullvalene_from_barcode(
    barcode: BVBarcode | BVTSBarcode,
    substituents: tuple[Substituent, ...],
    sanitize: bool = True,
    set_stereochemistry: bool = True
) -> Chem.Mol:
    """
    Builds a substituted bullvalene with the supplied substituents attached.

    If the bullvalene barcode supplied is a `BVTSBarcode` instance, a
    bullvalene transition state is built instead of an equilibrium/'stable-
    state' bullvalene.

    Args:
        barcode (BVBarcode | BVTSBarcode): `BVBarcode` or `BVTSBarcode`
            instance defining the substitution configuration of the bullvalene.
        substituents (tuple[Substituent, ...]): Tuple of `Substituent`
            instances defining the substituents to attach to the bullvalene
            and the barcode bit <-> attachment point mapping(s).
        sanitize (bool, optional): If `True`, the substituted bullvalene is
            sanitized by RDKit using `Chem.SanitizeMol()` before return.
            Defaults to `True`.
        set_stereochemistry (bool, optional): If `True`, the stereochemistry of
            the substituted bullvalene is set for compatibility with BULLVISO.
            Defaults to `True`.

    Returns:
        Chem.Mol: Substituted bullvalene.
    """

    _validate_barcode_substituent_compatibility(barcode, substituents)

    transition_state = isinstance(barcode, BVTSBarcode)
    bullvalene = load_bullvalene_from_library(
        transition_state = transition_state
    )

    substituted_bullvalene = _add_substituents(
        bullvalene,
        substituents,
        barcode
    )

    substituted_bullvalene = substituted_bullvalene.GetMol()

    if sanitize:
        Chem.SanitizeMol(substituted_bullvalene)

    if set_stereochemistry:
        stereo_map = (
            BULLVALENE_STEREO_MAP if not transition_state
            else BULLVALENE_TS_STEREO_MAP
        )
        _set_atom_stereochemistry(substituted_bullvalene, stereo_map)

    return substituted_bullvalene

@lru_cache(maxsize = 2)
def load_bullvalene_from_library(
    transition_state: bool = False
) -> Chem.Mol:
    """
    Loads a bullvalene (or, optionally, Cope rearrangement transition state)
    RDKit `Chem.Mol` instance from the `bullviso.structures` library with
    preset 3D/Cartesian coordinates.

    Args:
        transition_state (bool, optional): If `True`, a Cope rearrangement
            transition state bullvalene is returned rather than a stable-state
            bullvalene. Defaults to `False`.

    Raises:
        ValueError: If the RDKit `Chem.Mol` instance fails to load from the
            relevant mol (.sdf) file in the `bullviso.structures` library.

    Returns:
        Chem.Mol: Bullvalene.
    """
    
    mol_file_name = 'bv.sdf' if not transition_state else 'bv_ts.sdf'
    with resources.open_text('bullviso.structures', mol_file_name) as mol_file:
        bullvalene = Chem.MolFromMolBlock(mol_file.read())
    
    if bullvalene is None:
        raise ValueError(
            f'failed to load structure from \'{mol_file_name}\' in '
            f'`bullviso.structures`: check that the file exists and contains a '
            f'valid mol (.sdf) file'
        )
    
    return bullvalene

def _validate_barcode_substituent_compatibility(
    barcode: BVBarcode | BVTSBarcode,
    substituents: tuple[Substituent, ...]
) -> None:
    """
    Validates compatibility between the bullvalene barcode and substituents.

    This function checks that every barcode bit is associated with a
    corresponding substituent attachment point (i.e., there are no extra bits),
    and that no substituents have attachment points that are not associated
    with a corresponding barcode bit (i.e., there are no missing bits).

    Args:
        barcode (BVBarcode | BVTSBarcode): `BVBarcode` or `BVTSBarcode`
            instance defining the substitution configuration of the bullvalene.
        substituents (tuple[Substituent, ...]): Tuple of `Substituent`
            instances defining the substituents to attach to the bullvalene
            and the barcode bit <-> attachment point mapping(s).

    Raises:
        ValueError: If there are barcode bits that are not associated with a
            substituent attachment point (i.e., there are extra bits).
        ValueError: If there are substituent attachment points that are not
            associated with barcode bits (i.e., there are missing bits).
    """
    
    nonzero_barcode_bits = set(
        bit for bit in barcode.barcode_labels if bit != 0
    )

    used_barcode_bits: set[int] = set()
    for substituent in substituents:
        used_barcode_bits.update(substituent.barcode_bits)

    extra_barcode_bits = nonzero_barcode_bits - used_barcode_bits
    if extra_barcode_bits:
        raise ValueError(
            f'bullvalene barcode ({barcode}) contains extra bits '
            f'(extra bits = {extra_barcode_bits}) not associated with any '
            f'of the supplied substituent'
        )
    
    missing_barcode_bits = used_barcode_bits - nonzero_barcode_bits
    if missing_barcode_bits:
        raise ValueError(
            f'bullvalene barcode ({barcode}) is missing bits '
            f'(missing bits = {missing_barcode_bits}) associated with one or '
            f'more of the supplied substituent(s)'
        )

def _add_substituents(
    bullvalene: Chem.Mol,
    substituents: tuple[Substituent, ...],
    barcode: BVBarcode | BVTSBarcode
) -> Chem.RWMol:
    """
    Attaches multiple substituents to a bullvalene.

    Args:
        bullvalene (Chem.Mol): Bullvalene.
        substituents (tuple[Substituent, ...]): Tuple of `Substituent`
            instances defining the substituents to attach to the bullvalene
            and the barcode bit <-> attachment point mapping(s).
        barcode (BVBarcode | BVTSBarcode): `BVBarcode` or `BVTSBarcode`
            instance defining the substitution configuration of the bullvalene.

    Returns:
        Chem.RWMol: Bullvalene (mutable) with the substituents attached.
    """

    substituted_bullvalene = Chem.RWMol(bullvalene)
    for substituent in substituents:
        substituted_bullvalene = _add_substituent(
            substituted_bullvalene,
            substituent,
            barcode
        )

    return substituted_bullvalene

def _add_substituent(
    bullvalene: Chem.RWMol,
    substituent: Substituent,
    barcode: BVBarcode | BVTSBarcode
) -> Chem.RWMol:
    """
    Attaches a single substituent to a bullvalene.

    Args:
        bullvalene (Chem.RWMol): Bullvalene (mutable).
        substituent (Substituent): `Substituent` instance defining the
            substituent to attach to the bullvalene and the barcode bit <->
            attachment point mapping(s).
        barcode (BVBarcode | BVTSBarcode): `BVBarcode` or `BVTSBarcode`
            instance defining the substitution configuration of the bullvalene.

    Raises:
        ValueError: If the substituent SMILES string is invalid.

    Returns:
        Chem.RWMol: Bullvalene (mutable) with the substituent attached.
    """

    substituent_mol = Chem.MolFromSmiles(substituent.smiles)
    if substituent_mol is None:
        raise ValueError(
            f'invalid substituent SMILES string: got {substituent.smiles}'
        )

    bullvalene.InsertMol(substituent_mol)
    
    offset = (
        bullvalene.GetNumAtoms() - substituent_mol.GetNumAtoms()
    )
    
    for barcode_bit, attachment_idx in zip(
        substituent.barcode_bits, substituent.attachment_idxs
    ):
        attachment_idx_bullvalene = barcode.barcode_labels.index(barcode_bit)
        attachment_idx_substituent = attachment_idx + offset
        _remove_explicit_hydrogens(
            bullvalene, attachment_idx_bullvalene
        )
        bullvalene.AddBond(
            attachment_idx_bullvalene,
            attachment_idx_substituent,
            Chem.rdchem.BondType.SINGLE
        )

    return bullvalene

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