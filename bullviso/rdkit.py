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
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs, EmbedParameters
from rdkit.Geometry import rdGeometry

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def generate_confs(
        mol: Chem.Mol,
        forcefield: str = 'uff',
        prune_rms_thresh: float = 0.5,
        coord_map: dict[int, rdGeometry.Point3D] = None,
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
        coord_map(dict[int, rdGeometry.Point3D], optional): Coordinate map
            dictionary mapping atom indices to their 3D coordinates
            (represented as rdGeometry.Point3D instances); these atoms are
            fixed/frozen during the embedding procedure. Defaults to None.
        num_threads(int, optional): The number of threads to use for generating
            conformers; `num_threads` threads are used for each of the
            embedding and optimisation procedures. Defaults to 1.

    Returns:
        Chem.Mol: A molecule with embedded conformers.
    """
    
    params = getattr(Chem.rdDistGeom, "ETKDGv2")()
    params.pruneRmsThresh = prune_rms_thresh
    params.coordMap = {} if coord_map is None else coord_map
    params.numThreads = num_threads

    mol = embed_confs(
        mol, params = params
    )

    mol = optimise_confs(
        mol, forcefield = forcefield
    )

    mol = order_confs_by_energy(mol)

    return mol

def embed_confs(
    mol: Chem.Mol,
    n_confs: int = None,
    params: EmbedParameters = None
) -> Chem.Mol:
    """
    Embeds `n_confs` conformers of a molecule `mol`; optional parameters to
    control conformer embedding can be passed as an RDKit
    rdDistGeom.EmbedParameter object.

    Args:
        mol (Chem.Mol): Molecule.
        n_confs (int, optional): Number of conformers to embed; if `None`,
            the default is to embed 30 conformers if the number of rotatable
            bonds is less than 8, else 120 conformers if the number of
            rotatable bonds is 8 or greater. Defaults to `None`.
        params (EmbedParameters, optional): RDKit rdDistGeom.EmbedParameters
            object that holds optional parameters to control conformer
            embedding as attributes. Defaults to `None`.

    Returns:
        Chem.Mol: Molecule with embedded conformers.
    """
    
    if n_confs is None:
        n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        n_confs = 30 if n_rotatable_bonds < 8 else 120

    if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
        mol = Chem.AddHs(mol)

    EmbedMultipleConfs(
        mol, numConfs = n_confs, params = params
    )

    return mol

def optimise_confs(
    mol: Chem.Mol,
    forcefield: str = 'uff'
) -> Chem.Mol:
    """
    Optimises conformers of a molecule `mol` using a molecular mechanics / 
    forcefield method; the Universal Forcefield (UFF) and Merck Molecular
    Forcefield (MMFF) are available via RDKit.

    Args:
        mol (Chem.Mol): Molecule.
        forcefield (str, optional): Forcefield for conformer optimisation;
            choices are 'uff' and 'mmff'. Defaults to 'uff'.

    Returns:
        Chem.Mol: Molecule with forcefield-optimised conformers.
    """
    
    if mol.GetNumConformers() > 0:
        for conf in mol.GetConformers():
            if forcefield == 'mmff':
                energy = AllChem.MMFFOptimizeMolecule(
                    mol, confId = conf.GetId()
                )
            elif forcefield == 'uff':
                energy = AllChem.UFFOptimizeMolecule(
                    mol, confId = conf.GetId()
                )
            else:
                raise ValueError(
                    f'{forcefield} is not a recognised forcefield'
                )
            conf.SetDoubleProp('energy', energy)

    return mol

def order_confs_by_energy(
        mol: Chem.Mol,
) -> Chem.Mol:
    """
    Orders the forcefield-optimised conformers of a molecule `mol` in
    ascending order by energy; requires the `energy` property to be set for
    each conformer, i.e., accessible via `conf.GetProp('energy')` for each
    instance of `conf` returned via, e.g., `mol.GetConformers()`.

    Args:
        mol (Chem.Mol): Molecule.

    Returns:
        Chem.Mol: Molecule with forcefield-optimised conformers in ascending
            order by energy.
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

def get_coord_map(
    mol: Chem.Mol,
    conf_idx: int = -1,
    atom_idx: list[int] = None
) -> dict[int, rdGeometry.Point3D]:
    """
    Returns a coordinate map dictionary for a molecule `mol` that maps the
    atom indices to their 3D coordinates (represented as rdGeometry.Point3D
    instances).

    Args:
        mol (Chem.Mol): Molecule.
        conf_idx (int, optional): Index of the conformer to return the
            coordinate map for. Defaults to -1.
        atom_idx (list[int], optional): List of indices defining the atoms to
            include in the coordinate map; if `None`, all atoms are included
            in the coordinate map. Defaults to `None`.

    Returns:
        dict[int, rdGeometry.Point3D]: Coordinate map dictionary mapping atom
            indices to their 3D coordinates (represented as rdGeometry.Point3D
            instances).
    """
    
    if atom_idx is None:
        atom_idx = [i for i in range(mol.GetNumAtoms())]

    conf = mol.GetConformer(conf_idx)
    
    return {i: conf.GetAtomPosition(i) for i in atom_idx}
