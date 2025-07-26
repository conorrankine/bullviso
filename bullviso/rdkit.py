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
from rdkit.ForceField import rdForceField

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def generate_confs(
        mol: Chem.Mol,
        prune_rms_thresh: float = 0.5,
        coord_map: dict[int, rdGeometry.Point3D] = None,
        ff_type: str = 'mmff',
        constrained_ff_opt: bool = True,
        random_seed: int = -1,
        num_threads: int = 1
) -> Chem.Mol:
    """
    Embeds (using the ETKDGv2 distance geometry approach) and optimises (using
    a molecular mechanics / forcefield method) conformers of a molecule `mol`;
    the Merck Molecular Forcefield (MMFF) and Universal Forcefield (UFF) are
    available via RDKit.

    Args:
        mol (Chem.Mol): Molecule.
        prune_rms_thresh (float, optional): RMSD threshold for pruning
            conformers; conformers with RMSDs below the RMSD threshold are
            considered equivalent and are pruned. Defaults to 0.5 (Angstroem).
        coord_map (dict[int, rdGeometry.Point3D], optional): Coordinate map
            dictionary mapping atom indices to their 3D coordinates
            (represented as rdGeometry.Point3D instances); these atoms are
            fixed/frozen during the embedding procedure. Defaults to `None`.
        ff_type (str, optional): Forcefield type; choices are 'mmff' and 
            'uff'. Defaults to 'mmff'.
        constrained_ff_opt (bool, optional): Toggles constrained conformer
            optimisation; if `True`, and if `coord_map` is not `None`, the
            atoms that are fixed/frozen during the embedding procedure are
            also fixed/frozen during conformer optimisation.
        random_seed (int, optional): Number to use as the random seed for the
            embedding procedure; if -1, the random seed is obtained via
            random number generation. Defaults to -1.
        num_threads (int, optional): Number of threads to use in parallel
            processing operations. Defaults to 1.

    Returns:
        Chem.Mol: A molecule with embedded conformers.
    """
    
    params = getattr(Chem.rdDistGeom, "ETKDGv2")()
    
    params.SetCoordMap({} if coord_map is None else coord_map)
    
    params.pruneRmsThresh = prune_rms_thresh
    params.randomSeed = random_seed
    params.numThreads = num_threads

    mol = embed_confs(
        mol, params = params
    )

    if constrained_ff_opt and coord_map:
        fixed_atom_idx = [i for i in coord_map.keys()]
    else:
        fixed_atom_idx = None

    mol = optimise_confs(
        mol, ff_type = ff_type, fixed_atom_idx = fixed_atom_idx
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
    ff_type: str = 'mmff',
    fixed_atom_idx: list[int] = None
) -> Chem.Mol:
    """
    Optimises conformers of a molecule `mol` using a molecular mechanics / 
    forcefield method; the Merck Molecular Forcefield (MMFF) and Universal
    Forcefield (UFF) are available via RDKit.

    Args:
        mol (Chem.Mol): Molecule.
        ff_type (str, optional): Forcefield type; choices are 'mmff' and 
            'uff'. Defaults to 'mmff'.
        fixed_atom_idx (list[int], optional): List of atom indices for atoms
            to fix/freeze during conformer optimisation. Defaults to `None`.

    Returns:
        Chem.Mol: Molecule with forcefield-optimised conformers.
    """
    
    if mol.GetNumConformers() > 0:

        keep_conf_idxs = []

        for conf in mol.GetConformers():            
            ff = _get_forcefield(ff_type, mol, conf_idx = conf.GetId())
            if fixed_atom_idx:
                ff = _add_atomic_position_constraints(
                    ff, ff_type, fixed_atom_idx
                )                      
            opt_result = ff.Minimize()
            if opt_result == 0:
                conf.SetDoubleProp('energy', ff.CalcEnergy())
                keep_conf_idxs.append(conf.GetId())

        for conf_idx in reversed(range(mol.GetNumConformers())):
            if conf_idx not in keep_conf_idxs:
                mol.RemoveConformer(conf_idx)

    return mol

def _get_forcefield(
    ff_type: str,
    mol: Chem.Mol,
    conf_idx: int = -1
) -> rdForceField.ForceField:
    """
    Returns the specified type of forcefield for a molecule `mol` as an
    rdForceField.ForceField instance.

    Args:
        ff_type (str): Forcefield type; choices are 'mmff' and 'uff'.
        mol (Chem.Mol): Molecule.
        conf_idx (int, optional): Index of the conformer to return the
            specified type of forcefield for. Defaults to -1.

    Raises:
        ValueError: If `ff_type` is not either 'mmff' or 'uff'.

    Returns:
        rdForceField.ForceField: Forcefield.
    """
    
    forcefield_getter_functions = {
        'mmff': _get_mmff_forcefield,
        'uff': _get_uff_forcefield
    }

    try:
        forcefield_getter_function = (
            forcefield_getter_functions[ff_type]
        )
    except KeyError:
        raise ValueError(
            f'{ff_type} is not a recognised forcefield'
        ) from None
    
    return forcefield_getter_function(
        mol, conf_idx = conf_idx
    )

def _get_mmff_forcefield(
    mol: Chem.Mol,
    conf_idx: int = -1
) -> rdForceField.ForceField:
    """
    Returns a Merck Molecular Forcefield (MMFF) forcefield for a molecule
    `mol` as an rdForceField.ForceField instance.

    Args:
        mol (Chem.Mol): Molecule.
        conf_idx (int, optional): Index of the conformer to return the
            MMFF forcefield for. Defaults to -1.

    Returns:
        rdForceField.ForceField: MMFF forcefield.
    """
    
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
    return AllChem.MMFFGetMoleculeForceField(
        mol, mmff_props, confId = conf_idx
    )
                
def _get_uff_forcefield(
    mol: Chem.Mol,
    conf_idx: int = -1
) -> rdForceField.ForceField:
    """
    Returns a Universal Forcefield (UFF) forcefield for a molecule `mol` as
    an rdForceField.ForceField instance.

    Args:
        mol (Chem.Mol): Molecule.
        conf_idx (int, optional): Index of the conformer to return the
            UFF forcefield for. Defaults to -1.

    Returns:
        rdForceField.ForceField: UFF forcefield.
    """
    
    return AllChem.UFFGetMoleculeForceField(
        mol, confId = conf_idx
    )

def _add_atomic_position_constraints(
    ff: rdForceField.ForceField,
    ff_type: str,
    fixed_atom_idx: list[int]
) -> rdForceField.ForceField:
    """
    Adds atomic position constraints to a forcefield, e.g., to keep atoms
    fixed/frozen during conformer optimisation; the Merck Molecular Forcefield
    (MMFF) and Universal Forcefield (UFF) are supported.

    Args:
        ff (rdForceField.ForceField): Forcefield.
        ff_type (str): Forcefield type; choices are 'mmff' and 'uff'.
        fixed_atom_idx (list[int]): List of atom indices for atoms to apply
            atomic position constraints to.

    Raises:
        ValueError: If `ff_type` is not either 'mmff' or 'uff'.

    Returns:
        rdForceField.ForceField: Forcefield with atomic position constraints.
    """
    
    position_constraint_functions = {
        'mmff': ff.MMFFAddPositionConstraint,
        'uff': ff.UFFAddPositionConstraint
    }

    try:
        position_constraint_function = (
            position_constraint_functions[ff_type]
        )
    except KeyError:
        raise ValueError(
            f'{ff_type} is not a recognised forcefield'
        ) from None    

    for atom_idx in fixed_atom_idx:
        position_constraint_function(atom_idx, 0.0, 1.0E5)

    return ff

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
