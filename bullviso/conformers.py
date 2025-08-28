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
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
from rdkit.Chem.rdMolAlign import AlignMolConformers
from rdkit.Geometry import rdGeometry
from rdkit.ForceField import rdForceField
from rdkit.ML.Cluster import Butina
from .xtb_wrapper import XTBCalculator

###############################################################################
################################## CONSTANTS ##################################
###############################################################################

EMBED_METHOD_DEFAULT = 'ETKDGv3'

EMBED_N_CONFS_DEFAULT = {
     8 :  64,
    10 : 128,
    12 : 256,
}

EMBED_N_CONFS_DEFAULT_MAX = 512

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def generate_confs(
    mol: Chem.Mol,
    embed_n_confs: int = None,
    embed_rmsd_threshold: float = 0.0,
    embed_timeout: float = None,
    calculator_type: str = 'mmff',
    max_iter: int = 600,
    energy_threshold: float = 10.0,
    rmsd_threshold: float = 0.5,
    rmsd_atom_idxs: list[int] = None,
    coord_map: dict[int, rdGeometry.Point3D] = None,
    n_proc: int = 1,
    seed: int = -1,
) -> Chem.Mol:
    """
    Orchestrates a multi-step conformer generation, optimisation, and selection
    workflow to produce a diverse set of low-energy molecular conformations.

    The conformers of the input molecule `mol` are:

        1. embedded using the ETKDGv3 distance geometry algorithm;

        2. optimised using either the Merck Molecular Forcefield (MMFF),
           Universal Forcefield (UFF), or Extended Tight Binding (XTB)
           semiempirical quantum mechanical framework;

        3. filtered to remove high-energy instances that exceed a cutoff energy
           threshold relative to the lowest-energy conformer;

        4. clustered using the Butina algorithm to group structurally similar
           instances based on pairwise RMSD;

        5. selected from the Butina clusters, retaining only the lowest-energy
           instance belonging to each;

        6. sorted in ascending order by energy.

    Args:
        mol (Chem.Mol): Molecule.
        embed_n_confs (int, optional): Number of conformers to embed; if
            `None`, a default is determined based on the number of rotatable
            bonds in the molecule. Defaults to `None`.
        embed_rmsd_threshold (float, optional): RMSD threshold for deduplicating
            embeddings; conformers with RMSDs below the RMSD threshold are
            considered to be duplicates. Defaults to 0.0 (Angstroem).
        embed_timeout (float, optional): Timeout (in seconds) for conformer
            embedding. Defaults to `None`.
        calculator_type (str, optional): Calculator type; supported options are
            'mmff', 'uff', and 'xtb'. Defaults to 'mmff'.
        max_iter (int, optional): Maximum number of iterations for conformer
            optimisation. Defaults to 600.
        energy_threshold (float, optional): Maximum allowed energy difference
            (in kcal/mol) relative to the lowest-energy conformation.
            Defaults to 10.0 (kcal/mol).
        rmsd_threshold (float, optional): RMSD threshold for Butina clustering;
            conformers with RMSDs below the RMSD threshold are considered to
            belong to the same Butina cluster. Defaults to 0.5 (Angstroem).
        rmsd_atom_idxs (list[int], optional): List of atom indices defining the
            atoms to align pre-calculation of the pairwise RMSDs; if `None`,
            all atoms are used to align. Defaults to `None`.
        coord_map (dict[int, rdGeometry.Point3D], optional): Coordinate map
            dictionary mapping atom indices to their 3D coordinates
            (represented as rdGeometry.Point3D instances); these atoms are fixed
            during conformer embedding and optimisation. Defaults to `None`.
        n_proc (int, optional): Number of (parallel) processes for conformer
            embedding and optimisation. Defaults to 1.
        seed (int, optional): Seed for conformer embedding; if -1, the seed is
            obtained via pseudo-random number generation. Defaults to -1.

    Returns:
        Chem.Mol: Molecule with a diverse set of low-energy molecular
            conformations sorted in ascending order by energy.
    """

    mol = embed_confs(
        mol,
        n_confs = embed_n_confs,
        rmsd_threshold = embed_rmsd_threshold,
        timeout = embed_timeout,
        coord_map = coord_map,
        n_proc = n_proc,
        seed = seed
    )

    if coord_map:
        fixed_atom_idxs = set([i for i in coord_map.keys()])
    else:
        fixed_atom_idxs = None

    optimise_confs(
        mol,
        calculator_type = calculator_type,
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
            rmsd_threshold = rmsd_threshold
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
    rmsd_threshold: float = 0.0,
    timeout: float = None,
    coord_map: dict[int, rdGeometry.Point3D] = None,
    n_proc: int = 1,
    seed: int = -1
) -> Chem.Mol:
    """
    Embeds `n_confs` 3D conformers of a molecule `mol`.

    Notes:
        - The embedding algorithm uses the distance geometry method defined
          by the `EMBED_METHOD_DEFAULT` module-level variable.
        - The molecule is deep-copied before embedding; the original molecule
          is not modified.
        - If the molecule does not have explicit hydrogen atoms, these are
          added before embedding using RDKit's `Chem.AddHs()` function.        

    Args:
        mol (Chem.Mol): Molecule.
        n_confs (int, optional): Number of conformers to embed; if `None`, a
            default is determined based on the number of rotatable bonds in the
            molecule. Defaults to `None`.
        rmsd_threshold (float, optional): RMSD threshold for deduplicating
            embeddings; conformers with RMSDs below the RMSD threshold are
            considered to be duplicates. Defaults to 0.0 (Angstroem).
        timeout (float, optional): Timeout (in seconds) for conformer
            embedding. Defaults to `None`.
        coord_map (dict[int, rdGeometry.Point3D], optional): Coordinate map
            dictionary mapping atom indices to their 3D coordinates
            (represented as rdGeometry.Point3D instances); these atoms are
            fixed during conformer embedding. Defaults to `None`.
        n_proc (int, optional): Number of (parallel) processes for conformer
            embedding. Defaults to 1.
        seed (int, optional): Seed for conformer embedding; if -1, the seed is
            obtained via pseudo-random number generation. Defaults to -1.

    Returns:
        Chem.Mol: Molecule with embedded conformers.
    """

    mol = copy.deepcopy(mol)
    if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
        mol = Chem.AddHs(mol)

    if not n_confs:
        n_confs = _determine_n_confs_to_embed(mol)

    params = getattr(Chem.rdDistGeom, EMBED_METHOD_DEFAULT)()
    params.pruneRmsThresh = rmsd_threshold
    params.numThreads = n_proc
    params.randomSeed = seed
    if coord_map:
        params.SetCoordMap(coord_map)
    if timeout:
        params.timeout = timeout

    EmbedMultipleConfs(
        mol, numConfs = n_confs, params = params
    )

    return mol

def optimise_confs(
    mol: Chem.Mol,
    calculator_type: str = 'mmff',
    fixed_atom_idxs: set[int] = None,
    max_iter: int = 600
) -> None:
    """
    Optimises conformers of a molecule `mol`; the Merck Molecular Forcefield
    (MMFF), Universal Forcefield (UFF), and Extended Tight Binding (XTB)
    semiempirical quantum mechanical framework are available.

    Note:
        MMFF and UFF are available through RDKit, while XTB requires a working
        installation of XTB accessible via the system PATH.

    Args:
        mol (Chem.Mol): Molecule.
        calculator_type (str, optional): Calculator type; supported options are
            'mmff', 'uff', and 'xtb'. Defaults to 'mmff'.
        fixed_atom_idxs (set[int], optional): Set of atomic indices defining
            the fixed atoms. Defaults to `None`.
        max_iter (int, optional): Maximum number of iterations for conformer
            optimisation. Defaults to 600.
    """
    
    if mol.GetNumConformers() > 0:

        keep_conf_ids = []

        for conf in mol.GetConformers():            
            calculator = _get_calculator(
                calculator_type, mol, conf_id = conf.GetId()
            )
            if fixed_atom_idxs:
                _fix_atoms(
                    calculator, fixed_atom_idxs
                )                      
            opt_result = calculator.Minimize(maxIts = max_iter)
            if opt_result == 0:
                conf.SetDoubleProp('energy', calculator.CalcEnergy())
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
        atomIds = rmsd_atom_idxs
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
        atomIds = rmsd_atom_idxs
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

def _determine_n_confs_to_embed(
    mol: Chem.Mol,
) -> int:
    """
    Returns a default number of conformers to embed for a molecule `mol` based
    on the number of rotatable bonds in the molecule.
    
    The default number of conformers to embed is determined using predefined
    rotatable bond count thresholds in `EMBED_N_CONFS_DEFAULT`; if the
    rotatable bond count of the molecule exceeds all thresholds, a maximum
    value (`EMBED_N_CONFS_DEFAULT_MAX`) is returned.

    Args:
        mol (Chem.Mol): Molecule.

    Returns:
        int: Number of conformers to embed.
    """
    
    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    for threshold, n_confs in sorted(EMBED_N_CONFS_DEFAULT.items()):
        if n_rotatable_bonds <= threshold:
            return n_confs
    return EMBED_N_CONFS_DEFAULT_MAX

def _get_calculator(
    calculator_type: str,
    mol: Chem.Mol,
    conf_id: int = -1,
    **kwargs
) -> Union[rdForceField.ForceField, XTBCalculator]:
    """
    Returns the specified type of calculator for a molecule `mol` as either an
    rdForceField.ForceField or XTBCalculator instance.

    Args:
        calculator_type (str): Calculator type; supported options are 'mmff',
            'uff', and 'xtb'.
        mol (Chem.Mol): Molecule.
        conf_id (int, optional): Conformer ID to return the specified type of
            calculator for. Defaults to -1.
        **kwargs: Additional keyword arguments (kwargs) passed to the
            calculator class constructor.

    Raises:
        ValueError: If `calculator_type` is not one of 'mmff', 'uff', or 'xtb'.

    Returns:
        Union[rdForceField.ForceField, XTBCalculator]: Calculator.
    """
    
    calculators = {
        'mmff': {'function': _get_mmff_forcefield, 'accepts_kwargs': False},
        'uff':  {'function': _get_uff_forcefield,  'accepts_kwargs': False},
        'xtb':  {'function': _get_xtb_calculator,  'accepts_kwargs': True}
    }

    try:
        calculator_info = calculators[calculator_type]
    except KeyError:
        raise ValueError(
            f'{calculator_type} is not a supported calculator: supported '
            f'calculators include {{{", ".join(calculators)}}}'
        ) from None
    
    if calculator_info['accepts_kwargs']:
        return calculator_info['function'](
            mol, conf_id = conf_id, **kwargs
        )
    elif kwargs:
        raise TypeError(
            f'{calculator_type} does not accept additional keyword arguments: '
            f'{list(kwargs.keys())}'
        )
    else:
        return calculator_info['function'](
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

def _get_xtb_calculator(
    mol: Chem.Mol,
    conf_id: int = -1,
    **kwargs
) -> XTBCalculator:
    """
    Returns an XTB calculator for a molecule `mol` as a XTBCalculator instance.

    Args:
        mol (Chem.Mol): Molecule.
        conf_id (int, optional): Conformer ID to return the XTB calculator for.
            Defaults to -1.
        **kwargs: Additional keyword arguments (kwargs) passed to the
            calculator class constructor.

    Returns:
        XTBCalculator: XTB calculator.
    """
    
    return XTBCalculator(
        mol, conf_id = conf_id, **kwargs
    )

def _fix_atoms(
    calculator: Union[rdForceField.ForceField, XTBCalculator],
    fixed_atom_idxs: set[int]
) -> None:
    """
    Fixes atomic positions for a calculator by atomic index; atoms that are
    fixed do not have their Cartesian coordinates modified during a subsequent
    geometry optimisation using the calculator.

    Args:
        calculator (Union[rdForceField.ForceField, XTBCalculator]): Calculator.
        fixed_atom_idxs (set[int]): Set of atomic indices defining the fixed
            atoms.
    """

    for atom_idx in fixed_atom_idxs:
        calculator.AddFixedPoint(atom_idx)
