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

from enum import Enum
from rdkit import Chem

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

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

def get_conf_props(
    mol: Chem.Mol,
    prop: str,
    prop_type: str = 'str'
) -> list[tuple[int, str | float | int]]:
    """
    Returns a list of ([CONFORMER_ID], [PROPERTY_VALUE]) tuples for each
    conformer of the specified molecule.

    Args:
        mol (Chem.Mol): Molecule.
        prop (str): Property name to retrieve values for.
        prop_type (str, optional): Property type; one of {'str', 'double',
            'int'}. Defaults to 'str'.

    Raises:
        ValueError: If `prop_type` is not one of {'str', 'double', 'int'}.

    Returns:
        list[tuple[int, str | float | int]]: List of ([CONFORMER_ID],
            [PROPERTY_VALUE]) tuples.
    """

    getter = {
        'str': Chem.rdchem.Conformer.GetProp,
        'int': Chem.rdchem.Conformer.GetIntProp,
        'double': Chem.rdchem.Conformer.GetDoubleProp
    }.get(prop_type)

    if getter is None:
        raise ValueError(
            f'`prop_type` should be one of {{"str", "double", "int"}}; got '
            f'{prop_type!r}'
        )

    return [
        (conf.GetId(), getter(conf, prop)) for conf in mol.GetConformers()
    ]

def reorder_confs(
    mol: Chem.Mol,
    conf_ids: list[int]
) -> None:
    """
    Reorders conformers of the specified molecule to match the supplied
    sequence of conformer IDs.

    Args:
        mol (Chem.Mol): Molecule.
        conf_ids (list[int]): Ordered list of conformer IDs to retain.

    Raises:
        ValueError: If the supplied sequence of conformer IDs contains
            duplicated entries.
        ValueError: If any conformer IDs are missing from the molecule.
    """

    if len(conf_ids) != len(set(conf_ids)):
        raise ValueError(
            '`conf_ids` contains duplicated conformer IDs'
        )

    _validate_conf_ids(mol, conf_ids)
    
    conf_by_id = {
        conf.GetId(): Chem.Conformer(conf) for conf in mol.GetConformers()
    }

    mol.RemoveAllConformers()
    for conf_id in conf_ids:
        mol.AddConformer(conf_by_id[conf_id], assignId = True)

def remove_confs(
    mol: Chem.Mol,
    conf_ids: list[int]
) -> None:
    """
    Removes conformers of the specified molecule by conformer ID.

    Args:
        mol (Chem.Mol): Molecule.
        conf_ids (list[int]): List of conformer IDs to remove.

    Raises:
        ValueError: If any conformer IDs are missing from the molecule.
    """

    _validate_conf_ids(mol, conf_ids)

    conf_ids_to_remove = set(conf_ids)
    for conf in list(mol.GetConformers()):
        if conf.GetId() in conf_ids_to_remove:
            mol.RemoveConformer(conf.GetId())

def _validate_conf_ids(
    mol: Chem.Mol,
    conf_ids: list[int]
) -> None:
    """
    Validates that all supplied conformer IDs exist in the molecule.

    Args:
        mol (Chem.Mol): Molecule.
        conf_ids (list[int]): List of conformer IDs to validate.

    Raises:
        ValueError: If any conformer IDs are missing from the molecule.
    """

    conf_by_id = set(conf.GetId() for conf in mol.GetConformers())
    missing_conf_ids = sorted(set(conf_ids) - conf_by_id)
    if missing_conf_ids:
        raise ValueError(
            f'`conf_ids` contains conformer IDs that are missing from the '
            f'specified molecule: {{{", ".join(map(str, missing_conf_ids))}}}'
        )

def remove_explicit_hydrogens(
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

def set_atom_stereochemistry(
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

def set_bond_stereochemistry(
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
