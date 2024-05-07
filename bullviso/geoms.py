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

def generate_conformations(
        mol: Chem.Mol,
        m_confs: int = 10,
        rms_threshold: float = 0.5
) -> tuple:
    """
    Embeds `n_confs` conformations of a Chem.Mol molecule `mol` with an RMS
    threshold of `rms_threshold` and optimises the conformations using the
    Universal Forcefield (UFF); if a parallel installation of `rdkit` is
    available, these routines can run parallelised over `n_threads` threads.

    Args:
        mol (Chem.Mol): A molecule to generate conformations for.
        n_threads (int, optional): The number of threads to use if a parallel
            installation of `rdkit` is available. Defaults to 1.
        m_confs (int, optional): The number of conformations of `mol` to
            return. Defaults to 10.
        rms_threshold (float, optional): The RMS threshold for embedding the
            trial conformations of `mol`. Defaults to 0.5.

    Returns:
        unconverged (list): A list of binary elements (e.g. 1/0; True/False)
            indicating UFF convergence for each conformation.
        uff_energies (list): A list of UFF energies for each conformation.
    """
    
    params = getattr(Chem.rdDistGeom, "ETKDGv2")()
    params.pruneRmsThresh = rms_threshold

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
    
    conf_energies = [
        conf.GetProp('energy') for conf in mol.GetConformers()
    ]

    mol_copy = copy.deepcopy(mol)
    ordered_confs = [
        conf for _, conf in sorted(
            zip(conf_energies, mol_copy.GetConformers()), key = lambda x: x[0]
        )
    ]
    mol.RemoveAllConformers()
    for conf in ordered_confs[:m_confs]:
        mol.AddConformer(conf, assignId = True)

    return mol
