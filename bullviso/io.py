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

def mol_to_out_f(
    filepath: Path,
    filetype: str,
    mol: Chem.Mol,
    conf_idx: int = 0
) -> None:
    
    mol_to_out_f_ = {
        'xyz': mol_to_xyz_f,
        'gaussian': mol_to_gaussian_input_f,
        'orca': mol_to_orca_input_f
    }

    mol_to_out_f_[filetype](
        filepath = filepath,
        mol = mol,
        conf_idx = conf_idx
    )

    return None

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

def mol_to_gaussian_input_f(
    filepath: Path,
    mol: Chem.Mol,
    conf_idx: int = 0,
    method: str = 'PBE1PBE',
    basis: str = 'DEF2SVPP',
    charge: int = 0,
    multiplicity: int = 1,
    n_proc: int = 1,
    memory: int = 4000
) -> None:
    
    with open(filepath.with_suffix('.gjf'), 'w') as f:
        # write .chk filepath and resource specifications
        f.write(
            f'%CHK={filepath.stem}.chk\n%MEM={memory}MB\n%NPROC={n_proc}\n'
        )
        # write route line
        f.write(
            f'#N {method}/{basis} OPT FREQ\n\n'
        )
        # write title line
        f.write(
            f'CONFIGURATIONAL ISOMER: {filepath.stem}\n\n'
        )
        # write charge, multiplicity, and Cartesian coordinate lines
        f.write(
            f'{charge} {multiplicity}\n'
        )
        for i, atom in enumerate(mol.GetAtoms()):
            coord = mol.GetConformer(conf_idx).GetAtomPosition(i)
            coord_line_fmt = '{:<2}{:>14.8f}{:>14.8f}{:>14.8f}\n'
            f.write(
                coord_line_fmt.format(
                    atom.GetSymbol(), coord.x, coord.y, coord.z
                )
            )
        # write terminating blank line
        f.write(
            '\n'
        )

    return None

def mol_to_orca_input_f(
    filepath: Path,
    mol: Chem.Mol,
    conf_idx: int = 0,
    method: str = 'PBE0',
    basis: str = 'DEF2-SV(P)',
    charge: int = 0,
    multiplicity: int = 1,
    n_proc: int = 1,
    memory: int = 4000
) -> None:
    
    with open(filepath.with_suffix('.in'), 'w') as f:
        # write route line
        f.write(
            f'! {method} {basis} OPT FREQ\n\n'
        )
        # write charge, multiplicity, and Cartesian coordinate lines
        f.write(
            f'* XYZ {charge} {multiplicity}\n'
        )
        for i, atom in enumerate(mol.GetAtoms()):
            coord = mol.GetConformer(conf_idx).GetAtomPosition(i)
            coord_line_fmt = '{:<2}{:>14.8f}{:>14.8f}{:>14.8f}\n'
            f.write(
                coord_line_fmt.format(
                    atom.GetSymbol(), coord.x, coord.y, coord.z
                )
            )
        f.write(
            '*'
        )
