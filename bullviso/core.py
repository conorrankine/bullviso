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

# =============================================================================
#                               LIBRARY IMPORTS
# =============================================================================

import bullviso as bv
import datetime
from . import utils
from tqdm import tqdm
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from rdkit import Chem

# =============================================================================
#                                   CLASSES
# =============================================================================

@dataclass(frozen = True)
class BullvisoParams:
    sub_smiles: list[str]
    n_subs: list[int]
    sub_attach_idx: list[list[int]]
    transition_state: bool
    m_confs: int
    embed_n_confs: int
    embed_rmsd_threshold: float
    embed_timeout: int
    embed_seed: int
    calculator_type: str
    max_iter: int
    energy_threshold: float
    rmsd_threshold: float
    n_proc: int
    output_dir: Path
    output_filetype: str

    def __post_init__(self):
        
        if not utils.all_same_length(
            self.sub_smiles,
            self.n_subs,
            self.sub_attach_idx
        ):
            raise ValueError(
                f'`sub_smiles`, `n_subs`, and `sub_attach_idx` should have '
                f'equal length; got lists with length {len(self.sub_smiles)}, '
                f'{len(self.n_subs)}, and {len(self.sub_attach_idx)}'
        )

        # TODO: add additional parameter validation for future version

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

def _bullviso(
    params: BullvisoParams
) -> None:

    datetime_ = datetime.datetime.now()
    print(f'launched @ {datetime_.strftime("%H:%M:%S (%Y-%m-%d)")}\n')

    with resources.open_text('bullviso.assets.banners','banner.txt') as f:
        for line in f.readlines():
            print(line.rstrip())
    print('\n')

    print(' ' * 3 + f'{"SMILEs":<30} {"number":>10} {"attached @":>15}')
    print('-' * 60)
    for i, (sub_smile, n_sub, attach_idx) in enumerate(
        zip(params.sub_smiles, params.n_subs, params.sub_attach_idx), start = 1
    ):
        print(f'{i}. {sub_smile:<30} {str(n_sub):>10} {str(attach_idx):>15}')
    print('-' * 60 + '\n')

    transition_states = (False, True) if params.transition_state else (False,)

    for transition_state in transition_states:

        bullvalene = bv.substituents.load_bullvalene_from_library(
            transition_state = transition_state
        )

        substituents = bv.substituents.substituents_from_specifications(
            params.sub_smiles,
            params.n_subs,
            params.sub_attach_idx
        )

        barcode_type = bv.BVTSBarcode if transition_state else bv.BVBarcode
        canonical_barcode = barcode_type.from_substituents(substituents)
        print(f'canonical barcode: {canonical_barcode}\n')

        print('identifying unique bullvalene barcodes:')
        barcodes = set(
            barcode for barcode in tqdm(
                canonical_barcode.permutations(),
                ncols = 60,
                total = 3628800,
                bar_format='{l_bar}{bar}| [{elapsed}<{remaining}]'
            )
        )
        print('')

        coord_map = bv.conformers.get_coord_map(bullvalene)

        print('building bullvalene isomers:')
        for barcode in tqdm(
            barcodes,
            ncols = 60,
            bar_format='{l_bar}{bar}| [{elapsed}<{remaining}]'
        ):
            
            mol = bv.substituents.build_bullvalene_from_barcode(
                barcode,
                substituents
            )
            mol = bv.conformers.generate_confs(
                mol,
                embed_n_confs = params.embed_n_confs,
                embed_rmsd_threshold = params.embed_rmsd_threshold,
                embed_timeout = params.embed_timeout,
                embed_seed = params.embed_seed,
                calculator_type = params.calculator_type,
                max_iter = params.max_iter,
                coord_map = coord_map,
                energy_threshold = params.energy_threshold,
                rmsd_threshold = params.rmsd_threshold,
                align_atom_idxs = [i for i in range(10)],
                n_proc = params.n_proc
            )
            conf_ids = _get_conf_ids(mol, max_conf_ids = params.m_confs)
            for conf_id in conf_ids:
                _write_conf_to_file(
                    mol,
                    conf_id,
                    params.output_dir,
                    params.output_filetype
                )
        print('')

    datetime_ = datetime.datetime.now()
    print(f'finished @ {datetime_.strftime("%H:%M:%S (%Y-%m-%d)")}')

def _get_conf_ids(
    mol: Chem.Mol,
    max_conf_ids: int = None
) -> list[int]:
    """
    Returns a list of conformer IDs for an RDKit `Chem.Mol` molecule.

    Args:
        mol (Chem.Mol): Molecule.
        max_conf_ids (int, optional): Maximum number of conformer IDs to
            return; if None, or if `max_conf_ids` exceeds the number of stored
            conformers, all conformer IDs are returned. Defaults to None.

    Returns:
        list[int]: List of conformer IDs.
    """
    
    try:
        confs = [conf for conf in mol.GetConformers()]
        if max_conf_ids is None:
            return [conf.GetId() for conf in confs]
        return [conf.GetId() for conf in confs[:max_conf_ids]]
    except Exception as e:
        raise ValueError(
            f'failed to retrieve conformer IDs: {e}'
        )

def _write_conf_to_file(
    mol: Chem.Mol,
    conf_id: int,
    output_dir: Path,
    output_filetype: str
) -> None:
    """
    Writes a single conformer of the supplied RDKit `Chem.Mol` molecule to a
    file in an organised directory structure.
    
    The schema for the organised directory structure is:

    `output_dir` / [BARCODE] / [BARCODE]_[CONF_ID] / 
        [BARCODE]_[CONF_ID].[OUTPUT_FILETYPE]

    Note: `conf_id` is 0-based internally but output as 1-based in filenames
    for user-friendliness, e.g., `conf_id = 0` -> [CONF_ID] = '001'.
    
    Args:
        mol (Chem.Mol): Molecule.
        conf_id (int): Conformer ID to write to file.
        output_dir (Path): Output directory; the output directory is created
            if it does not already exist.
        output_filetype (str): Output filetype.

    Raises:
        ValueError: If the molecule does not have a `barcode` property.
        ValueError: If the molecule does not have a `transition_state` property.
    """

    if not mol.HasProp('barcode'):
        raise ValueError(
            'molecule does not have a \'barcode\' property; set one using '
            '`mol.SetProp(\'barcode\', [BARCODE])`'
        )
    barcode = mol.GetProp('barcode')

    if not mol.HasProp('transition_state'):
        raise ValueError(
            'molecule does not have a \'transition_state\' property; set one '
            'using `mol.SetBoolProp(\'transition_state\', [TRUE/FALSE])`'
        )
    transition_state = mol.GetBoolProp('transition_state')

    try:
        conf_dir = (
            output_dir
            / ('minima' if not transition_state else 'transition_states')
            / barcode
            / f'{barcode}_{(conf_id+1):03d}'
        )
        conf_dir.mkdir(parents = True, exist_ok = True)
        conf_file = conf_dir / f'{barcode}_{(conf_id+1):03d}'
        bv.io.mol_to_file(
            conf_file,
            mol,
            conf_id = conf_id,
            filetype = output_filetype,
        )
    except Exception as e:
        raise ValueError(
            f'failed to write conformer(s) to file: {e}'
        )
    
# =============================================================================
#                                     EOF
# =============================================================================
