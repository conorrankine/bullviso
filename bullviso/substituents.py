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
from rdkit import Chem
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
        attach_idxs (list[int]): [#TODO: DESCRIPTION]
        barcode_bits (list[int]): [#TODO: DESCRIPTION]   
    """

    smiles: str
    attach_idxs: list[int]
    barcode_bits: list[int]

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

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
    
    for barcode_bit, attach_idx in zip(
        substituent.barcode_bits, substituent.attach_idxs
    ):
        attach_idx_bullvalene = barcode.barcode_labels.index(barcode_bit)
        attach_idx_substituent = attach_idx + offset
        _remove_explicit_hydrogens(
            bullvalene, attach_idx_bullvalene
        )
        bullvalene.AddBond(
            attach_idx_bullvalene,
            attach_idx_substituent,
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