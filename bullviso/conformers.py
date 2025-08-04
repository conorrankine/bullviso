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
from typing import Union
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs, EmbedParameters
from rdkit.Chem.rdMolAlign import AlignMolConformers
from rdkit.Geometry import rdGeometry
from rdkit.ForceField import rdForceField
from rdkit.ML.Cluster import Butina
from .xtb_wrapper import XTBOptimiser

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
        fixed_atom_idxs = set([i for i in coord_map.keys()])
    else:
        fixed_atom_idxs = None

    optimise_confs(
        mol,
        ff_type = ff_type,
        fixed_atom_idxs = fixed_atom_idxs,
        max_iter = max_iter
    )

    filter_low_energy_confs(
        mol,
        energy_threshold = energy_threshold
    )

    select_cluster_representatives(
        mol,
        clusters = cluster_confs(
            mol,
            rmsd_threshold = rmsd_threshold,
            rmsd_atom_idxs = rmsd_atom_idxs
        )
    )

    order_confs_by_energy(
        mol
    )

    align_confs(
        mol,
        rmsd_atom_idxs = rmsd_atom_idxs
    )

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

    mol = copy.deepcopy(mol)
    if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
        mol = Chem.AddHs(mol)

    EmbedMultipleConfs(
        mol, numConfs = n_confs, params = params
    )

    return mol

def optimise_confs(
    mol: Chem.Mol,
    ff_type: str = 'mmff',
    fixed_atom_idxs: set[int] = None,
    max_iter: int = 600
) -> None:
    """
    Optimises conformers of a molecule `mol` using a molecular mechanics / 
    forcefield method; the Merck Molecular Forcefield (MMFF) and Universal
    Forcefield (UFF) are available via RDKit.

    Args:
        mol (Chem.Mol): Molecule.
        ff_type (str, optional): Forcefield type; choices are 'mmff' and 
            'uff'. Defaults to 'mmff'.
        fixed_atom_idxs (set[int], optional): Set of atomic indices defining
            the fixed atoms. Defaults to `None`.
        max_iter (int, optional): Maximum number of iterations for conformer
            optimisation. Defaults to 600.
    """
    
    if mol.GetNumConformers() > 0:

        keep_conf_ids = []

        for conf in mol.GetConformers():            
            ff = _get_optimiser(ff_type, mol, conf_id = conf.GetId())
            if fixed_atom_idxs:
                _fix_atoms(
                    ff, fixed_atom_idxs
                )                      
            opt_result = ff.Minimize(maxIts = max_iter)
            if opt_result == 0:
                conf.SetDoubleProp('energy', ff.CalcEnergy())
                keep_conf_ids.append(conf.GetId())

        _prune_confs(mol, keep_conf_ids)

def filter_low_energy_confs(
    mol: Chem.Mol,
    energy_threshold: float = 10.0
) -> None:
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
    """
    
    if mol.GetNumConformers() > 0:

        keep_conf_ids = []

        min_energy = min(
            conf.GetDoubleProp('energy') for conf in mol.GetConformers()
        )

        for conf in mol.GetConformers():
            if (conf.GetDoubleProp('energy') - min_energy) < energy_threshold:
                keep_conf_ids.append(conf.GetId())

        _prune_confs(mol, keep_conf_ids)

def align_confs(
    mol: Chem.Mol,
    rmsd_atom_idxs: list[int] = None
) -> None:
    """
    Aligns conformers of a molecule `mol` using the first conformer (i.e., the
    conformer at index zero) as the reference conformation.

    Args:
        mol (Chem.Mol): Molecule.
        rmsd_atom_idxs (list[int], optional): List of atom indices defining the
            atoms to use as alignment points; if `None`, all atoms are used as
            alignment points. Defaults to `None`.
    """
    
    AlignMolConformers(
        mol,
        atomIds = rmsd_atom_idxs,
        reflect = True
    )

def cluster_confs(
    mol: Chem.Mol,
    rmsd_threshold: float = 0.5,
    rmsd_atom_idxs: list[int] = None
) -> tuple[tuple[int]]:
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
        tuple[tuple[int]]: Tuple of Butina clusters where each Butina cluster
            is represented as a tuple of conformer IDs.
    """    

    rms_matrix = AllChem.GetConformerRMSMatrix(
        mol,
        atomIds = rmsd_atom_idxs,
    )

    clusters = Butina.ClusterData(
        rms_matrix,
        mol.GetNumConformers(),
        rmsd_threshold,
        isDistData = True
    )

    conf_ids = [conf.GetId() for conf in mol.GetConformers()]

    conf_id_clusters = (
        tuple(tuple(conf_ids[i] for i in cluster) for cluster in clusters)
    )

    return conf_id_clusters

def select_cluster_representatives(
    mol: Chem.Mol,
    clusters: tuple[tuple[int]]
) -> None:
    """
    Removes conformers of a molecule `mol` such that only the lowest-energy
    conformation belonging to each Butina cluster in `clusters` is retained.

    Absolute conformer energies (in kcal/mol) are expected to be stored as
    properties under the key 'energy' and, consequently, accessible via
    `conf.GetDoubleProp('energy')` for each conformer `conf`.

    Args:
        mol (Chem.Mol): Molecule.
        clusters (tuple[tuple[int]]): Tuple of Butina clusters where each
            Butina cluster is a tuple of conformer IDs.
    """
    
    keep_conf_ids = []

    energies = {
        conf.GetId(): conf.GetDoubleProp('energy')
        for conf in mol.GetConformers()
    }

    for cluster in clusters:
        keep_conf_ids.append(
            min(cluster, key = lambda conf_id: energies[conf_id])
        )
        
    _prune_confs(mol, keep_conf_ids)

def order_confs_by_energy(
        mol: Chem.Mol,
) -> None:
    """
    (Re)orders conformers of a molecule `mol` in ascending order by energy.

    Absolute conformer energies (in kcal/mol) are expected to be stored as
    properties under the key 'energy' and, consequently, accessible via
    `conf.GetDoubleProp('energy')` for each conformer `conf`.

    Args:
        mol (Chem.Mol): Molecule.
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
    keep_conf_ids: list[int]
) -> None:
    """
    Removes conformers of a molecule `mol` by ID such that only conformers with
    IDs listed in `keep_conf_ids` are retained.

    Args:
        mol (Chem.Mol): Molecule.
        keep_conf_ids (list[int]): List of conformer IDs defining the set of
            conformers to retain.
    """
    
    for conf in list(mol.GetConformers()):
        if conf.GetId() not in keep_conf_ids:
            mol.RemoveConformer(conf.GetId())

def _get_optimiser(
    optimiser_type: str,
    mol: Chem.Mol,
    conf_id: int = -1
) -> Union[rdForceField.ForceField, XTBOptimiser]:
    """
    Returns the specified type of optimiser for a molecule `mol` as either an
    rdForceField.ForceField or XTBOptimiser instance.

    Args:
        optimiser_type (str): Optimiser type; supported options are 'mmff',
            'uff', and 'xtb'.
        mol (Chem.Mol): Molecule.
        conf_id (int, optional): Conformer ID to return the specified type of
            optimiser for. Defaults to -1.

    Raises:
        ValueError: If `optimiser_type` is not one of 'mmff', 'uff', or 'xtb'.

    Returns:
        Union[rdForceField.ForceField, XTBOptimiser]: Optimiser.
    """
    
    optimiser_getters = {
        'mmff': _get_mmff_forcefield,
        'uff': _get_uff_forcefield,
        'xtb': _get_xtb_optimiser
    }

    try:
        optimiser_getter = (
            optimiser_getters[optimiser_type]
        )
    except KeyError:
        raise ValueError(
            f'{optimiser_type} is not a supported optimiser: supported '
            f'optimisers include {{{", ".join(optimiser_getters)}}}'
        ) from None
    
    return optimiser_getter(
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

def _get_xtb_optimiser(
    mol: Chem.Mol,
    conf_id: int = -1,
    **kwargs
) -> XTBOptimiser:
    """
    Returns an XTB optimiser for a molecule `mol` as an XTBOptimiser instance.

    Args:
        mol (Chem.Mol): Molecule.
        conf_id (int, optional): Conformer ID to return the XTB optimiser for.
            Defaults to -1.

    Returns:
        XTBOptimiser: XTB optimiser.
    """
    
    return XTBOptimiser(
        mol, conf_id = conf_id, **kwargs
    )

def _fix_atoms(
    optimiser: Union[rdForceField.ForceField, XTBOptimiser],
    fixed_atom_idxs: set[int]
) -> None:
    """
    Fixes atomic positions for an optimiser by atomic index; atoms that are
    fixed do not have their Cartesian coordinates modified during a subsequent
    geometry optimisation using the optimiser.

    Args:
        optimiser (Union[rdForceField.ForceField, XTBOptimiser]): Optimiser.
        fixed_atom_idxs (set[int]): Set of atomic indices defining the fixed
            atoms.
    """

    for atom_idx in fixed_atom_idxs:
        optimiser.AddFixedPoint(atom_idx)
