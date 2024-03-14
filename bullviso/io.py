"""
BULLVISO
Copyright (C) 2022  Conor D. Rankine

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

###############################################################################
############################### LIBRARY IMPORTS ###############################
###############################################################################

from pathlib import Path
from rdkit import Chem

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def mol_to_xyz_f(
    filepath: Path,
    mol: Chem.Mol,
    conf_idx: int = 0
) -> None:
    
    with open(filepath.with_suffix('.xyz'), 'w') as f:
        # write number of atoms
        f.write(
            f'{mol.GetNumAtoms()}\n\n'
        )
        # write Cartesian coordinate lines
        for i, atom in enumerate(mol.GetAtoms()):
            coord = mol.GetConformer(conf_idx).GetAtomPosition(i)
            coord_line_fmt = '{:<2}{:>14.8f}{:>14.8f}{:>14.8f}\n'
            f.write(
                coord_line_fmt.format(
                    atom.GetSymbol(), coord.x, coord.y, coord.z
                )
            )
    
    return None
