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

from argparse import ArgumentParser, Namespace
from bullviso.graphs import compose_bullvalene_supergraph_from_smiles

###############################################################################
############################## ARGUMENT PARSING ###############################
###############################################################################

def parse_args() -> Namespace:
    """
    Parses command line arguments for `bullviso:cli.py`.

    Returns:
        argparse.Namespace: Parsed command line arguments as an
        argparse.Namespace object that holds the arguments as attributes.
    """

    p = ArgumentParser()

    p.add_argument(
        'sub_smiles', type = str, nargs = '+',
        help = 'SMILES string representation for each unique substituent'
    )
    p.add_argument(
        '--n_sub', '-n', type = int, nargs = '+', default = 1,
        help = ('number of each unique substituent to add')
    )
    p.add_argument(
        '--sub_attachment_idx', '-a', type = int, nargs = '+', default = 1,
        help =('atomic index of the substituent-bullvalene attachment point '
            'for each unique substituent')
    )
    # p.add_argument('--m_confs', '-m', type = int, default = 1,
    #     help = ('number of conformational isomers to generate')
    # )
    # p.add_argument('--forcefield', '-ff', type = str, default = 'uff',
    #     choices = ('uff', 'mmff'),
    #     help = ('forcefield for optimising conformational isomers')
    # )
    # p.add_argument('--prune_rms_thresh', '-rmsd', type = float, default = 0.5,
    #     help = ('RMSD threshold for pruning conformational isomers')
    # )
    # p.add_argument('--num_threads', '-nt', type = int, default = 1,
    #     help = ('number of threads for generating conformational isomers')
    # )
    # p.add_argument('--out_f_type', '-o', type = str, default = 'xyz',
    #     choices = ('xyz', 'gaussian', 'orca'),
    #     help = ('file type for output geometries')
    # )

    args = p.parse_args()

    return args

###############################################################################
################################ MAIN FUNCTION ################################
###############################################################################

def main():

    args = parse_args()

    super_G = compose_bullvalene_supergraph_from_smiles(
        sub_smiles = args.sub_smiles
    )

    # func_group_smile = args.func_group_smile
    # print(f'>> functional group SMILE: {func_group_smile}')
    # func_group = Chem.MolFromSmiles(func_group_smile)

    # print(f'>> {args.n_func_groups} functional groups')

    # confcodes = bullviso.confcodes.gen_confcodes(args.n_func_groups)
    # print(f'>> {len(confcodes)} unique structural isomer(s)\n')

    # print('>> constructing structural isomer(s)...')
    # for confcode in tqdm.tqdm(confcodes, ncols = 80):       
    #     func_bullvalene = bullviso.geoms.functionalise(
    #         bullvalene,
    #         func_group,
    #         confcode,
    #         func_group_attach_idx = args.func_group_attach_idx
    #     )
    #     func_bullvalene = bullviso.geoms.generate_confs(
    #         func_bullvalene,
    #         forcefield = args.forcefield,
    #         prune_rms_thresh = args.prune_rms_thresh,
    #         num_threads = args.num_threads
    #     )
    #     m_confs = min(
    #         func_bullvalene.GetNumConformers(), args.m_confs
    #     )
    #     if m_confs > 0:
    #         confcode_str = bullviso.utils.tuple_to_str(confcode)
    #         d = Path(f'./{confcode_str}')
    #         if not d.is_dir():
    #             d.mkdir()
    #         for conf_idx in range(m_confs):
    #             out_d = d / f'./{confcode_str}_{conf_idx+1:03d}'
    #             if not out_d.is_dir():
    #                 out_d.mkdir()
    #             out_f = out_d / f'./{confcode_str}_{conf_idx+1:03d}'
    #             bullviso.io.mol_to_out_f(
    #                 out_f,
    #                 args.out_f_type,
    #                 func_bullvalene,
    #                 conf_idx = conf_idx
    #             )

################################################################################
############################## PROGRAM STARTS HERE #############################
################################################################################

if __name__ == '__main__':
    main()

################################################################################
############################### PROGRAM ENDS HERE ##############################
################################################################################
