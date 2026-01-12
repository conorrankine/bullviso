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

# =============================================================================
#                                     EOF
# =============================================================================
