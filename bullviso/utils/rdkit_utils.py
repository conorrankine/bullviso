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

from rdkit import Chem

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

# =============================================================================
#                                     EOF
# =============================================================================
