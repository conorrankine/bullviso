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

###############################################################################
############################### LIBRARY IMPORTS ###############################
###############################################################################

import subprocess
import tempfile
from pathlib import Path
from rdkit import Chem
from rdkit.Geometry import rdGeometry

###############################################################################
################################### CLASSES ###################################
###############################################################################

class XTBOptimiser:
    """_
    XTB wrapper that implements and exposes the same interface as RDKit
    forcefield (rdForceField.ForceField) objects.
    """

    def __init__(
        self,
        mol: Chem.Mol,
        conf_id: int = -1,
        method: str = 'GFN2-xTB',
        xtb_path: str = 'xtb',
        n_proc: int = 1
    ):
        
        self.mol = mol
        self.conf_id = conf_id
        self.conf = mol.GetConformer(conf_id)

        self.charge = sum(
            atom.GetFormalCharge() for atom in mol.GetAtoms()
        )
        self.uhf = sum(
            atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()
        )

        self.method = method

        self.xtb_path = xtb_path

        self.n_proc = n_proc

    def CalcEnergy(
        self
    ) -> float:
        
        try:
            energy, _ = self._run_xtb_calculation()
            return energy
        except Exception as e:
            print(f'XTB energy calculation failed: {e}')
            return float('inf')
        
    def CalcGradient(
        self
    ) -> None:
        
        raise NotImplementedError(
            'support for XTB gradients coming in a future version of BULLVISO'
        )
        
    def Minimize(
        self,
        maxIts: int = 600
    ) -> int:
        
        try:
            energy, coords = self._run_xtb_calculation(
                minimize = True,
                max_iter = maxIts
            )
            for i, atom_coords in enumerate(coords):
                self.conf.SetAtomPosition(i, atom_coords)
            return 0
        except Exception as e:
            print(f'XTB geometry optimisation failed: {e}')
            return 1

    def _run_xtb_calculation(
        self,
        minimize: bool = False,
        max_iter: int = 600
    ) -> float:
        
        with tempfile.TemporaryDirectory() as tmpdir:
            Chem.MolToXYZFile(
                self.mol,
                f'{tmpdir}/mol.xyz',
                confId = self.conf_id
            )

            cmd = [self.xtb_path, f'{tmpdir}/mol.xyz']

            if minimize:
                cmd.extend(['--opt'])
                cmd.extend(['--cycles', str(max_iter)])

            if self.method == 'GFN1-xTB':
                cmd.extend(['--gfn', '1'])

            if self.charge != 0:
                cmd.extend(['--chrg', str(self.charge)])
            if self.uhf != 0:
                cmd.extend(['--uhf', str(self.uhf)])

            cmd.extend(['--parallel', str(self.n_proc)])

            result = subprocess.run(
                cmd,
                cwd = Path(tmpdir),
                capture_output = True,
                text = True
            )

            if result.returncode != 0:
                raise RuntimeError(
                    f'XTB did not finish successfully - spooling stderr: '
                    f'{result.stderr}'
                )
            
            energy = _get_xtb_energy_au(result.stdout)

            if minimize:
                coords = _get_coords_from_xyz_file(f'{tmpdir}/xtbopt.xyz')
            else:
                coords = None

            return energy, coords

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def _get_xtb_energy_au(
    xtb_output: str
) -> float:
    
    if not xtb_output.strip():
        raise RuntimeError(
            'XTB output is empty'
        )
    for line in reversed(xtb_output.split('\n')):
        if 'TOTAL ENERGY' in line:
            for item in line.split():
                try:
                    return float(item)
                except ValueError:
                    continue
    energy_lines = [
        line.strip() for line in xtb_output.split('\n')
        if 'TOTAL ENERGY' in line
    ]
    if energy_lines:
        raise RuntimeError(
            'couldn\'t read energy from XTB output: \'TOTAL ENERGY\' is in '
            'the output, but the energy value couldn\'t be parsed'
        )
    else:
        raise RuntimeError(
            'couldn\'t read energy from XTB output: \'TOTAL ENERGY\' is not '
            'in the output'
        )

def _get_gradient_from_grad_file(
    grad_file: str
) -> None:

    raise NotImplementedError(
        'support for XTB gradients coming in a future version of BULLVISO'
    )

def _get_coords_from_xyz_file(
    xyz_file: str
) -> tuple[rdGeometry.Point3D]:
    
    try:
        mol = Chem.MolFromXYZFile(xyz_file)
        if mol is None:
            raise RuntimeError(
                f'couldn\'t read a molecule from {xyz_file} '
                f'[Chem.MolFromXYZFile({xyz_file}) returned `None`]'
            )
        conf = mol.GetConformer()
        return [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    except Exception as e:
        raise RuntimeError(
            f'couldn\'t read coordinates from {xyz_file}: {e}'
        )
