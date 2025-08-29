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
import networkx as nx
import tqdm
import datetime
from . import utils
from argparse import Namespace
from pathlib import Path

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

def run_bullviso(
    args: Namespace
) -> None:

    datetime_ = datetime.datetime.now()
    print(f'launched @ {datetime_.strftime("%H:%M:%S (%Y-%m-%d)")}\n')

    header_f = Path(__file__).parent / 'assets' / 'banners' / 'banner.txt'
    with open(header_f, 'r') as f:
        for line in f.readlines():
            print(line.rstrip())
    print('\n')

    print(' ' * 3 + f'{"SMILEs":<30} {"number":>10} {"attached @":>15}')
    print('-' * 60)
    for i, (sub_smile, n_sub, attach_idx) in enumerate(
        zip(args.sub_smiles, args.n_subs, args.sub_attach_idx), start = 1
    ):
        print(f'{i}. {sub_smile:<30} {str(n_sub):>10} {str(attach_idx):>15}')
    print('-' * 60 + '\n')

    sub_smiles = utils.repeat_list_elements(
        args.sub_smiles, args.n_subs
    )

    sub_attach_idx = utils.repeat_list_elements(
        args.sub_attach_idx, args.n_subs
    )

    canonical_barcode = bv.barcodes.create_barcode(
        sub_smiles,
        sub_attach_idx
    )
    print(f'canonical barcode: {canonical_barcode}\n')

    connectivity_map = {
        i: f'sub{sub_n+1}_{sub_attach_idx_}'
            for i, (sub_n, sub_attach_idx_) in enumerate(
                utils.iterate_and_index(sub_attach_idx), start = 1
            )
    }

    print('identifying inequivalent permutations...')
    barcodes = set(
        barcode for barcode in tqdm.tqdm(
            canonical_barcode.permutations(), ncols = 60, total = 3628800
        )
    )
    print('...done!\n')

    model_bullvalene = bv.io.sdf_to_mol(
        Path(__file__).parent / 'structures' / 'bv.sdf'
    )

    bullvalene_G = bv.graphs.mol_to_graph(
        model_bullvalene, node_label_prefix = 'bullvalene_'
    )

    bv.graphs.set_atom_stereochemistry(
        bullvalene_G,
        atom_stereo_map = {
            'bullvalene_1'  : 'ccw',
            'bullvalene_4'  : 'cw',
            'bullvalene_7'  : 'ccw',
            'bullvalene_10' : 'cw'
        }
    )

    sub_G = [
        bv.graphs.smiles_to_graph(
            sub_smile, node_label_prefix = f'sub{i}_'
        ) for i, sub_smile in enumerate(sub_smiles, start = 1)
    ]
        
    super_G = nx.compose_all([bullvalene_G, *sub_G])

    coord_map = bv.conformers.get_coord_map(
        model_bullvalene
    )

    print('generating geometries and writing output...')
    for barcode in tqdm.tqdm(barcodes, ncols = 60):
        super_G_ = super_G.copy()
        for i, bit in enumerate(barcode.barcode, start = 1):
            if bit != 0:
                super_G_.add_edge(
                    f'bullvalene_{i}', connectivity_map[bit]
                )
        mol = bv.graphs.graph_to_mol(
            super_G_
        )
        mol = bv.conformers.generate_confs(
            mol,
            embed_n_confs = args.embed_n_confs,
            embed_rmsd_threshold = args.embed_rmsd_threshold,
            embed_timeout = args.embed_timeout,
            embed_seed = args.embed_seed,
            calculator_type = args.calculator_type,
            max_iter = args.max_iter,
            coord_map = coord_map,
            energy_threshold = args.energy_threshold,
            rmsd_threshold = args.rmsd_threshold,
            rmsd_atom_idxs = [i for i in range(10)],
            n_proc = args.n_proc
        )
        confs = list(mol.GetConformers())
        for conf_idx in range(min(args.m_confs, mol.GetNumConformers())):
            conf = confs[conf_idx]
            output_dir = args.output_dir / (
                f'./{barcode}/{barcode}_{conf_idx+1:03d}/'
            )
            if not output_dir.is_dir():
                output_dir.mkdir(parents = True)
            bv.io.mol_to_file(
                output_dir / f'./{barcode}_{conf_idx+1:03d}',
                mol,
                filetype = args.output_filetype,
                conf_id = conf.GetId()
            )
    print('...done!\n')

    datetime_ = datetime.datetime.now()
    print(f'finished @ {datetime_.strftime("%H:%M:%S (%Y-%m-%d)")}')

# =============================================================================
#                                     EOF
# =============================================================================