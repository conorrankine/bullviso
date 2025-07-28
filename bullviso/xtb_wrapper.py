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
            energy = self._run_xtb_calculation()
            return energy
        except Exception as e:
            print(f'XTB energy calculation failed: {e}')
            return float('inf')

    def _run_xtb_calculation(
        self
    ) -> float:
        
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_file = Path(f'{tmpdir}/mol.xyz')
            Chem.MolToXYZFile(self.mol, xyz_file, confId = self.conf_id)

            cmd = [self.xtb_path, xyz_file]

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

            return energy

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def _get_xtb_energy_au(
    xtb_output: str
) -> float:
        
    for line in reversed(xtb_output.split('\n')):
        if 'TOTAL ENERGY' in line:
            for item in line.split():
                try:
                    return float(item)
                except ValueError:
                    continue

    raise RuntimeError(
        'couldn\'t read the energy from the XTB output'
    )
