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
from rdkit.ML.Cluster import Butina

###############################################################################
################################## CONSTANTS ##################################
###############################################################################

EMBED_METHOD = 'ETKDGv3'
EMBED_MAXATTEMPTS = 1000
EMBED_PRUNERMSTHRESH = 0.5

ROTATABLE_BOND_THRESHOLD = 8
MIN_CONFS = 60
MAX_CONFS = 300

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def generate_confs(
    mol: Chem.Mol,
    ff_type: str = 'mmff',
    max_iter: int = 600,
    coord_map: dict[int, rdGeometry.Point3D] = None,
    energy_threshold: float = 10.0,
    rmsd_threshold: float = 0.5,
    rmsd_atom_idxs: list[int] = None,
    seed: int = -1,
    n_proc: int = 1
) -> Chem.Mol:
    """
    Orchestrates a multi-step conformer generation, optimisation, and selection
    workflow to produce a diverse set of low-energy molecular conformations.

    The conformers of the input molecule `mol` are:

        1. embedded using the ETKDGv3 distance geometry algorithm;

        2. optimised using either the Merck Molecular Forcefield (MMFF) or
           Universal Forcefield (UFF);

        3. filtered to remove high-energy instances that exceed a cutoff energy
           threshold relative to the lowest-energy conformer;

        4. clustered using the Butina algorithm to group structurally similar
           instances based on pairwise RMSD;

        5. selected from the Butina clusters, retaining only the lowest-energy
           instance belonging to each;

        6. sorted in ascending order by energy.

    Args:
        mol (Chem.Mol): Molecule.
        ff_type (str, optional): Forcefield type for conformer optimisation;
            choices are 'mmff' and 'uff'. Defaults to 'mmff'.
        max_iter (int, optional): Maximum number of iterations for conformer
            optimisation. Defaults to 600.
        coord_map (dict[int, rdGeometry.Point3D], optional): Coordinate map
            dictionary mapping atom indices to their 3D coordinates
            (represented as rdGeometry.Point3D instances); these atoms are
            fixed/frozen during conformer optimisation. Defaults to `None`.
        energy_threshold (float, optional): Maximum allowed energy difference
            (in kcal/mol) relative to the lowest-energy conformation.
            Defaults to 10.0 (kcal/mol).
        rmsd_threshold (float, optional): RMSD threshold for Butina clustering;
            conformers with RMSDs below the RMSD threshold are considered to
            belong to the same Butina cluster. Defaults to 0.5 (Angstroem).
        rmsd_atom_idxs (list[int], optional): List of atom indices defining the
            atoms to align pre-calculation of the pairwise RMSDs; if `None`,
            all atoms are used to align. Defaults to `None`.
        seed (int, optional): Seed for conformer embedding; if -1, the seed is
            obtained via pseudo-random number generation. Defaults to -1.
        n_proc (int, optional): Number of parallel processes for conformer
            embedding and optimisation; if 1, conformer embedding and
            optimisation is serial. Defaults to 1.

    Returns:
        Chem.Mol: Molecule with a diverse set of low-energy molecular
            conformations sorted in ascending order by energy.
    """
    
    params = getattr(Chem.rdDistGeom, EMBED_METHOD)()
    params.maxAttempts = EMBED_MAXATTEMPTS
    params.pruneRmsThresh = EMBED_PRUNERMSTHRESH
    params.coordMap = {} if coord_map is None else coord_map
    params.randomSeed = seed
    params.numThreads = n_proc

    mol = embed_confs(
        mol, params = params
    )

    if coord_map:
        fixed_atom_idx = [i for i in coord_map.keys()]
    else:
        fixed_atom_idx = None

    mol = optimise_confs(
        mol,
        ff_type = ff_type,
        fixed_atom_idx = fixed_atom_idx,
        max_iter = max_iter
    )

    mol = filter_low_energy_confs(
        mol,
        energy_threshold = energy_threshold
    )

    mol = select_cluster_representatives(
        *cluster_confs(
            mol,
            rmsd_threshold = rmsd_threshold,
            rmsd_atom_idxs = rmsd_atom_idxs
        )
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
        n_confs = (
            MIN_CONFS if n_rotatable_bonds < ROTATABLE_BOND_THRESHOLD
            else MAX_CONFS
        )

    if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
        mol = Chem.AddHs(mol)

    EmbedMultipleConfs(
        mol, numConfs = n_confs, params = params
    )

    return mol

def optimise_confs(
    mol: Chem.Mol,
    ff_type: str = 'mmff',
    fixed_atom_idx: list[int] = None,
    max_iter: int = 600
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
        max_iter (int, optional): Maximum number of iterations for conformer
            optimisation. Defaults to 600.

    Returns:
        Chem.Mol: Molecule with forcefield-optimised conformers.
    """
    
    if mol.GetNumConformers() > 0:

        keep_conf_idxs = []

        for conf_idx, conf in enumerate(mol.GetConformers()):            
            ff = _get_forcefield(ff_type, mol, conf_id = conf.GetId())
            if fixed_atom_idx:
                ff = _add_atomic_position_constraints(
                    ff, ff_type, fixed_atom_idx
                )                      
            opt_result = ff.Minimize(maxIts = max_iter)
            if opt_result == 0:
                conf.SetDoubleProp('energy', ff.CalcEnergy())
                keep_conf_idxs.append(conf_idx)

        _prune_confs(mol, keep_conf_idxs)

    return mol

def filter_low_energy_confs(
    mol: Chem.Mol,
    energy_threshold: float = 10.0
) -> Chem.Mol:
    """
    Removes conformers of a molecule `mol` with energies greater than the
    specified threshold `energy_threshold` (in kcal/mol) relative to the
    lowest-energy conformation.

    Absolute conformer energies (in kcal/mol) are expected to be stored as
    properties under the key 'energy' and, consequently, accessible via
    `conf.GetDoubleProp('energy')` for each conformer `conf`.

    Args:
        mol (Chem.Mol): Molecule.
        energy_threshold (float, optional): Maximum allowed energy difference
            (in kcal/mol) relative to the lowest-energy conformation.
            Defaults to 10.0 (kcal/mol).

    Returns:
        Chem.Mol: Molecule with only low-energy conformers retained.
    """
    
    keep_conf_idxs = []

    min_energy = min(
        conf.GetDoubleProp('energy') for conf in mol.GetConformers()
    )

    for conf_idx, conf in enumerate(mol.GetConformers()):
        if (conf.GetDoubleProp('energy') - min_energy) <= energy_threshold:
            keep_conf_idxs.append(conf_idx)

    _prune_confs(mol, keep_conf_idxs)

    return mol

def cluster_confs(
    mol: Chem.Mol,
    rmsd_threshold: float = 0.5,
    rmsd_atom_idxs: list[int] = None
) -> tuple[Chem.Mol, tuple[tuple[int]]]:
    """
    Clusters conformers of a molecule `mol` based on their pairwise RMSDs
    using the Butina algorithm; as a side effect, the conformers of the
    molecule are left in an aligned state.

    Args:
        mol (Chem.Mol): Molecule.
        rmsd_threshold (float, optional): RMSD threshold for Butina clustering;
            conformers with RMSDs below the RMSD threshold are considered to
            belong to the same Butina cluster. Defaults to 0.5 (Angstroem).
        rmsd_atom_idxs (list[int], optional): List of atom indices defining the
            atoms to align pre-calculation of the pairwise RMSDs; if `None`,
            all atoms are used to align. Defaults to `None`.

    Returns:
        tuple[Chem.Mol, tuple[tuple[int]]]: Molecule with aligned conformers,
            and tuple of Butina clusters where each Butina cluster is a tuple
            of conformer indices.
    """    
    
    rms_matrix = AllChem.GetConformerRMSMatrix(
        mol,
        atomIds = rmsd_atom_idxs
    )

    clusters = Butina.ClusterData(
        rms_matrix,
        mol.GetNumConformers(),
        rmsd_threshold,
        isDistData = True
    )

    return mol, clusters

def select_cluster_representatives(
    mol: Chem.Mol,
    clusters: tuple[tuple[int]]
) -> Chem.Mol:
    """
    Removes conformers of a molecule `mol` such that only the lowest-energy
    conformation belonging to each Butina cluster in `clusters` is retained.

    Absolute conformer energies (in kcal/mol) are expected to be stored as
    properties under the key 'energy' and, consequently, accessible via
    `conf.GetDoubleProp('energy')` for each conformer `conf`.

    Args:
        mol (Chem.Mol): Molecule.
        clusters (tuple[tuple[int]]): Tuple of Butina clusters where each
            Butina cluster is a tuple of conformer indices.

    Returns:
        Chem.Mol: Molecule with only the lowest-energy conformation belonging
            to each Butina cluster retained.
    """
    
    keep_conf_idxs = []

    energies = [
        conf.GetDoubleProp('energy') for conf in mol.GetConformers()
    ]

    for cluster in clusters:
        keep_conf_idxs.append(
            min(cluster, key = lambda i: energies[i])
        )
        
    _prune_confs(mol, keep_conf_idxs)

    return mol

def order_confs_by_energy(
        mol: Chem.Mol,
) -> Chem.Mol:
    """
    (Re)orders conformers of a molecule `mol` in ascending order by energy.

    Absolute conformer energies (in kcal/mol) are expected to be stored as
    properties under the key 'energy' and, consequently, accessible via
    `conf.GetDoubleProp('energy')` for each conformer `conf`.

    Args:
        mol (Chem.Mol): Molecule.

    Returns:
        Chem.Mol: Molecule with conformers in ascending order by energy.
    """
       
    energies = [
        conf.GetDoubleProp('energy') for conf in mol.GetConformers()
    ]

    mol_ = copy.deepcopy(mol)
    
    ordered_confs = [
        conf for _, conf in sorted(
            zip(energies, mol_.GetConformers()),
            key = lambda x: x[0]
        )
    ]
    
    mol.RemoveAllConformers()
    for conf in ordered_confs:
        mol.AddConformer(conf, assignId = True)

    return mol

def get_coord_map(
    mol: Chem.Mol,
    conf_id: int = -1,
    atom_idx: list[int] = None
) -> dict[int, rdGeometry.Point3D]:
    """
    Returns a coordinate map dictionary for a molecule `mol` that maps the
    atom indices to their 3D coordinates (represented as rdGeometry.Point3D
    instances).

    Args:
        mol (Chem.Mol): Molecule.
        conf_id (int, optional): Conformer ID to return the coordinate map
            for. Defaults to -1.
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

    conf = mol.GetConformer(conf_id)
    
    return {i: conf.GetAtomPosition(i) for i in atom_idx}

def _prune_confs(
    mol: Chem.Mol,
    keep_conf_idxs: list[int]
) -> Chem.Mol:
    """
    Removes conformers of a molecule `mol` by index such that only conformers
    corresponding to indices in `keep_conf_idxs` are retained.

    Args:
        mol (Chem.Mol): Molecule.
        keep_conf_idxs (list[int]): List of conformer indices defining the set
            of conformers to retain.

    Returns:
        Chem.Mol: Molecule with only the conformers corresponding to indices in
            `keep_conf_idxs` retained.
    """
    
    for conf_idx in reversed(range(mol.GetNumConformers())):
        if conf_idx not in keep_conf_idxs:
            mol.RemoveConformer(conf_idx)

    return mol

def _get_forcefield(
    ff_type: str,
    mol: Chem.Mol,
    conf_id: int = -1
) -> rdForceField.ForceField:
    """
    Returns the specified type of forcefield for a molecule `mol` as an
    rdForceField.ForceField instance.

    Args:
        ff_type (str): Forcefield type; choices are 'mmff' and 'uff'.
        mol (Chem.Mol): Molecule.
        conf_id (int, optional): Conformer ID to return the specified type of
            forcefield for. Defaults to -1.

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
        mol, conf_id = conf_id
    )

def _get_mmff_forcefield(
    mol: Chem.Mol,
    conf_id: int = -1
) -> rdForceField.ForceField:
    """
    Returns a Merck Molecular Forcefield (MMFF) forcefield for a molecule
    `mol` as an rdForceField.ForceField instance.

    Args:
        mol (Chem.Mol): Molecule.
        conf_id (int, optional): Conformer ID to return the MMFF forcefield
            for. Defaults to -1.

    Returns:
        rdForceField.ForceField: MMFF forcefield.
    """
    
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
    return AllChem.MMFFGetMoleculeForceField(
        mol, mmff_props, confId = conf_id
    )
                
def _get_uff_forcefield(
    mol: Chem.Mol,
    conf_id: int = -1
) -> rdForceField.ForceField:
    """
    Returns a Universal Forcefield (UFF) forcefield for a molecule `mol` as
    an rdForceField.ForceField instance.

    Args:
        mol (Chem.Mol): Molecule.
        conf_id (int, optional): Conformer ID to return the UFF forcefield
            for. Defaults to -1.

    Returns:
        rdForceField.ForceField: UFF forcefield.
    """
    
    return AllChem.UFFGetMoleculeForceField(
        mol, confId = conf_id
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
