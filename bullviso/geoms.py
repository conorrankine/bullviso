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
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem import rdMolDescriptors

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def functionalise(
    mol: Chem.Mol,
    func_group: Chem.Mol,
    confcode: tuple,
    func_group_attach_idx: int = 1
) -> Chem.Mol:
    """
    Functionalises a Chem.Mol molecule `mol` with the Chem.Mol functional group
    `func_group` according to a configuration code tuple `confcode`; N 
    (0 <= N <= `mol.GetNumAtoms()`) instances of `func_group` are attached to
    `mol` via single bonds at the position(s) indicated in `confcode`.

    Args:
        mol (Chem.Mol): A molecule to functionalise.
        func_group (Chem.Mol): A functional group.
        confcode (tuple): A configuration code encoding the functionalisation
            of `mol` as a tuple of binary (eg. 1/0) elements with length
            `mol.GetNumAtoms()`; each binary element indicates whether
            or not that (functionalisable) site in `mol` is functionalised.
        func_group_attach_idx (int): An attachment index indicating the
            position on the functional group `func_group` to be attached to
            `mol` via single bonds at the position(s) indicated in `confcode`.
            The atomic numbering scheme starts at 1. Defaults to 1.

    Returns:
        Chem.Mol: A functionalised molecule.
    """

    if len(confcode) != mol.GetNumAtoms():
        raise ValueError('the configuration code `confcode` is incompatible '\
            'with the molecule `mol`')
    
    if func_group_attach_idx > func_group.GetNumAtoms():
        raise ValueError('the attachment index `func_group_attach_idx` is '\
            'incompatible with the functional group `func_group`'
        )

    func_sites_i = tuple(
        i for i in range(len(confcode)) if confcode[i]
    )

    func_sites_j = tuple(
        mol.GetNumAtoms() + (i * func_group.GetNumAtoms()) + (
            func_group_attach_idx - 1) for i in range(len(func_sites_i))
    )

    func_mol = Chem.Mol(mol)
    
    for func_site_i, func_site_j in zip(func_sites_i, func_sites_j):
        func_mol = Chem.EditableMol(
            Chem.CombineMols(func_mol, func_group)
        )
        func_mol.AddBond(
            func_site_i, func_site_j, order = Chem.rdchem.BondType.SINGLE
        )
        func_mol = func_mol.GetMol()

    Chem.SanitizeMol(func_mol)

    func_mol = Chem.AddHs(func_mol)

    return func_mol

def generate_confs(
        mol: Chem.Mol,
        prune_rms_thresh: float = 0.5,
        num_threads: int = 1
) -> tuple:
    """
    Embeds `k` conformers of a Chem.Mol molecule `mol` and optimises the
    conformers using the Universal Forcefield (UFF). For molecules with less
    than eight rotatable bonds, `k` = 30; for molecules with more than eight
    rotatable bonds, `k` = 120.

    Args:
        mol (Chem.Mol): A molecule to embed conformers for.
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

    rdDistGeom.EmbedMultipleConfs(
        mol, numConfs = n_confs, params = params
    )

    rdForceFieldHelpers.UFFOptimizeMoleculeConfs(
        mol, maxIters = 120
    )

    for conf in mol.GetConformers():
        ff = rdForceFieldHelpers.UFFGetMoleculeForceField(
            mol, confId = conf.GetId()
        )
        conf.SetDoubleProp(
            'energy', ff.CalcEnergy()
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
