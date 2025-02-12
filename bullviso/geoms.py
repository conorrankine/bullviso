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

import copy
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem import rdMolDescriptors

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def generate_confs(
        mol: Chem.Mol,
        forcefield: str = 'uff',
        prune_rms_thresh: float = 0.5,
        num_threads: int = 1
) -> Chem.Mol:
    """
    Embeds `k` conformers of a Chem.Mol molecule `mol` and optimises the
    conformers using a forcefield `forcefield`. For molecules with less
    than eight rotatable bonds, `k` = 30; for molecules with more than eight
    rotatable bonds, `k` = 120.

    Args:
        mol (Chem.Mol): A molecule to embed conformers for.
        forcefield (str, optional): The forcefield to use for optimising the
            generated conformers; options are the Universal Forcefield ('uff')
            or the Merck Molecular Forcefield ('mmff'). Defaults to 'uff'.
        prune_rms_thresh (float, optional): The RMSD threshold for pruning the
            generated conformers; conformers below the RMSD threshold are
            considered the same and are pruned. Defaults to 0.5 Angstrom.
        num_threads(int, optional): The number of threads to use for generating
            conformers; `num_threads` threads are used for each of the
            embedding and optimisation procedures. Defaults to 1.

    Returns:
        Chem.Mol: A molecule with embedded conformers.
    """
    
    params = getattr(Chem.rdDistGeom, "ETKDGv2")()
    params.pruneRmsThresh = prune_rms_thresh
    params.numThreads = num_threads

    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_confs = 30 if n_rotatable_bonds < 8 else 120

    if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
        mol = Chem.AddHs(mol)

    rdDistGeom.EmbedMultipleConfs(
        mol, numConfs = n_confs, params = params
    )

    _ff_optimise = {
        'uff': rdForceFieldHelpers.UFFOptimizeMoleculeConfs,
        'mmff': rdForceFieldHelpers.MMFFOptimizeMoleculeConfs
    }
    ff_optimise = _ff_optimise.get(forcefield)

    ff_optimise(
        mol, numThreads = num_threads
    )

    _ff = {
        'uff': rdForceFieldHelpers.UFFGetMoleculeForceField,
        'mmff': rdForceFieldHelpers.MMFFGetMoleculeForceField
    }
    ff = _ff.get(forcefield)

    if forcefield == 'mmff':
        mol_props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)

    for conf in mol.GetConformers():
        if forcefield == 'mmff':
            conf_ff = ff(
                mol, mol_props, confId = conf.GetId()
            )
        else:
            conf_ff = ff(
                mol, confId = conf.GetId()
            )            
        conf.SetDoubleProp(
            'energy', conf_ff.CalcEnergy()
        )

    mol = order_confs_by_energy(mol)

    return mol

def order_confs_by_energy(
        mol: Chem.Mol,
) -> Chem.Mol:
    """
    Orders the conformers of a Chem.Mol molecule `mol` from lowest to highest
    energy; requires the 'energy' property to be set for each conformer, i.e.
    accessible via `conf.GetProp('energy')` for each instance of `conf`
    returned via `mol.GetConformers()`.

    Args:
        mol (Chem.Mol): A molecule with embedded conformers.

    Returns:
        Chem.Mol: A molecule with embedded conformers ordered from lowest to
            highest energy.
    """
       
    conf_energies = [
        conf.GetProp('energy') for conf in mol.GetConformers()
    ]

    mol_ = copy.deepcopy(mol)
    
    ordered_confs = [
        conf for _, conf in sorted(
            zip(conf_energies, mol_.GetConformers()),
            key = lambda x: x[0]
        )
    ]
    
    mol.RemoveAllConformers()
    for conf in ordered_confs:
        mol.AddConformer(conf, assignId = True)

    return mol
