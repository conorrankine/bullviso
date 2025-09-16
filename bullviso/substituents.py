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
from itertools import count
from rdkit import Chem
from .utils import all_same_length
from .barcodes import BVBarcode, BVTSBarcode

__all__ = [
    "Substituent",
    "build_bullvalene_from_barcode",
    "load_bullvalene_from_library"
]

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
    num_substituents: list[int],
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
        num_substituents (list[int]): Sequence of integers (one per
            substituent) defining the substituent count, i.e., the number of
            times that the corresponding substituent should appear. 
        attachment_idxs (list[list[int]]): Sequence of integer sequences (one
            per substituent) defining the atom indices of the attachment
            point(s) for the corresponding substituent.

    Raises:
        ValueError: If `substituent_smiles`, `num_substituents`, and
            `attachment_idxs` are not all of equal length.
        ValueError: If any element of `num_substituents` is < 0.
        ValueError: If any element of `attachment_idxs` is empty.
        ValueError: If more than 9 attachment points are specified (summed over
            all substituents).

    Returns:
        tuple[Substituent, ...]: Tuple of `Substituent` instances with unique
            bullvalene barcode bit values assigned to each attachment point.
    """

    if not all_same_length(
        substituent_smiles, num_substituents, attachment_idxs
    ):
        raise ValueError(
            f'`substituent_smiles`, `num_substituents`, and `attachment_idxs` '
            f'should all be of equal length; got sequences with lengths '
            f'{len(substituent_smiles)}, {len(num_substituents)}, and '
            f'{len(attachment_idxs)} (respectively)'
        )
    
    substituents: list[Substituent] = []

    barcode_bit_value = count(start = 1)

    for smiles_, num_substituents_, attachment_idxs_ in zip(
        substituent_smiles, num_substituents, attachment_idxs
    ):
        if num_substituents_ <= 0:
            raise ValueError(
                f'`num_substituents` should be positive: got '
                f'{num_substituents_} for SMILES = \'{smiles_}\''
            )
        if not attachment_idxs_:
            raise ValueError(
                f'`attachment_idxs` should not be empty: got '
                f'{attachment_idxs_} for SMILES = \'{smiles_}\''
            )
        for _ in range(num_substituents_):
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
            f'{num_attachment_idxs}'
        )

    return tuple(substituents)

def build_bullvalene_from_barcode(
    barcode: BVBarcode | BVTSBarcode,
    substituents: tuple[Substituent, ...],
    sanitize: bool = True
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

    Returns:
        Chem.Mol: Substituted bullvalene.
    """

    bullvalene = load_bullvalene_from_library(
        transition_state = isinstance(barcode, BVTSBarcode)
    )

    substituted_bullvalene = _add_substituents(
        bullvalene,
        substituents,
        barcode
    )

    substituted_bullvalene = substituted_bullvalene.GetMol()

    if sanitize:
        Chem.SanitizeMol(substituted_bullvalene)

    return substituted_bullvalene

def load_bullvalene_from_library(
    transition_state: bool = False
) -> Chem.Mol:
    
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
    """

    atom = mol.GetAtomWithIdx(atom_idx)
    atom.SetNumExplicitHs(0)
    atom.UpdatePropertyCache()

# =============================================================================
#                                     EOF
# =============================================================================