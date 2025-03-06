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

import bullviso as bv
import tqdm
import importlib
import datetime
import ast
from . import utils
from argparse import ArgumentParser, ArgumentTypeError, Namespace
from pathlib import Path
from typing import Union

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
        '--n_subs', '-n', type = _int_or_list_of_ints, default = [1],
        help = ('number of each unique substituent to add')
    )
    p.add_argument(
        '--sub_attach_idx', '-a', type = _int_or_list_of_ints, default = [1],
        help =('atomic index of the substituent-bullvalene attachment point '
            'for each unique substituent')
    )
    p.add_argument('--m_confs', '-m', type = int, default = 1,
        help = ('maximum number of conformational isomers to generate')
    )
    p.add_argument('--forcefield', '-ff', type = str, default = 'uff',
        choices = ('uff', 'mmff'),
        help = ('forcefield for optimising conformational isomers')
    )
    p.add_argument('--prune_rms_thresh', '-rmsd', type = float, default = 0.5,
        help = ('RMSD threshold for pruning conformational isomers')
    )
    p.add_argument('--num_threads', '-nt', type = int, default = 1,
        help = ('number of threads for optimising conformational isomers '
            'in multithreaded/parallel processes')
    )
    p.add_argument('--out_f_type', '-o', type = str, default = 'xyz',
        choices = ('xyz', 'sdf', 'gaussian', 'orca'),
        help = ('file type for outputting geometries')
    )

    args = p.parse_args()

    if len(args.sub_smiles) > 1:
        if len(args.n_subs) == 1:
            args.n_subs = args.n_subs * len(args.sub_smiles)
        if len(args.sub_attach_idx) == 1:
            args.sub_attach_idx = args.sub_attach_idx * len(args.sub_smiles)

    _validate_args(args)

    return args

def _int_or_list_of_ints(
    input: Union[int, list[Union[int, list[int]]]]
) -> list[Union[int, list[int]]]:
    """
    Custom argument type for a command line argument to allow the user to
    pass an integer, list of integers, or list of integers with nested
    (sub)lists of integers as an input value.

    Args:
        input (Union[int, list[Union[int, list[int]]]]): Input value.

    Raises:
        ArgumentTypeError: If `input` is not an integer, list of integers, or
            list of integers with nested (sub)lists of integers.

    Returns:
        list[Union[int, list[int]]: List of integers with or without nested
            (sub)lists of integers; a length-1 list if `input` is an integer,
            else a length-n list (where n = len(`input`)).
    """
    
    try:
        input = ast.literal_eval(input)
        if isinstance(input, int):
            return [input]
        elif isinstance(input, list):
            return _validate_list(input)
        else:
            raise ArgumentTypeError(
                'argument should be an integer, list of integers, or list of '
                'integers with nested (sub)lists of integers'
            )
    except (ValueError, SyntaxError):
        raise ArgumentTypeError(
            f'invalid argument: {input}'
        )
    
def _validate_list(
    input: list
) -> list[Union[int, list[int]]]:
    """
    Validates an input list recursively to ensure that it contains only integer
    elements or nested (sub)lists of integer elements.

    Args:
        input (list): Input list.

    Raises:
        ArgumentTypeError: If `input` contains elements that are not integers
            or nested (sub)lists, or if any nested sublist contains
            non-integer elements.

    Returns:
        list[Union[int, list[int]]]: List of integers and/or nested (sub)lists
            of integers; maintains the initial structure of `input`.
    """
    
    validated_list = []

    for i in input:
        if isinstance(i, int):
            validated_list.append(i)
        elif isinstance(i, list):
            validated_list.append(_validate_list(i))
        else:
            raise ArgumentTypeError(
                'list should contain only integer elements or lists of '
                'integer elements'
            ) 
    
    return validated_list
    
def _validate_args(
    args: Namespace
):
    """
    Validates command line arguments for `bullviso:cli.py`.

    Args:
        args (Namespace): Parsed command line arguments as an
        argparse.Namespace object that holds the arguments as attributes.

    Raises:
        ValueError: If `args.sub_smiles`, `args.n_subs`, and
            `args.sub_attach_idx` are not of equal length; if the sum of the
            integer elements in `args.n_subs` is greater than 9; if there are
            greater than 9 elements in `args.sub_attach_idx` (counting
            elements in the outer list and any nested (sub)lists); or if
            either `args.n_subs` or `args.sub_attach_idx` contained
            null/zero-valued elements.
    """

    if not (
        len(args.sub_smiles) == len(args.n_subs) == len(args.sub_attach_idx)
    ):
        raise ValueError(
            '`args.sub_smiles`, `args.n_subs`, and `args.sub_attach_idx` '
            'should have the same length'
        )
    if sum(args.n_subs) > 9:
        raise ValueError(
            '`args.n_subs` is too large; support is only available for up '
            'to (and including) 9 substituents'
        )
    if utils.count_list_elements(
        utils.repeat_list_elements(args.sub_attach_idx, args.n_subs)
    ) > 9:
        raise ValueError(
            'too many attachment indices in `args.sub_attach_idx`; support is '
            'only available for up to (and including) 9 substituents'
        )
    if 0 in args.n_subs:
        raise ValueError(
            '`args.n_subs` cannot contain null/zero-valued elements'
        )
    if 0 in args.sub_attach_idx:
        raise ValueError(
            '`args.sub_attach_idx` cannot contain null/zero-valued elements'
        )

###############################################################################
################################ MAIN FUNCTION ################################
###############################################################################

def main():

    datetime_ = datetime.datetime.now()
    print(f'launched @ {datetime_.strftime("%H:%M:%S (%Y-%m-%d)")}\n')

    args = parse_args()

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

    super_G = bv.graphs.compose_bullvalene_supergraph_from_smiles(
        sub_smiles = sub_smiles
    )

    connectivity_map = {
        i: f'sub{sub_n+1}_{sub_attach_idx_}'
            for i, (sub_n, sub_attach_idx_) in enumerate(
                utils.iterate_and_index(sub_attach_idx), start = 1
            )
    }

    canonical_barcode = bv.barcodes.create_barcode(
        sub_smiles,
        sub_attach_idx
    )

    print(f'canonical barcode: {canonical_barcode}\n')

    print('identifying inequivalent permutations...')
    barcodes = set(
        barcode for barcode in tqdm.tqdm(
            canonical_barcode.permutations(), ncols = 60, total = 3628800
        )
    )
    print('...done!\n')

    model_bullvalene = bv.io.sdf_to_mol(
        importlib.resources.files('bullviso.structures').joinpath('bv.sdf')
    )

    coord_map = bv.rdkit.get_coord_map(
        model_bullvalene
    )

    print('generating geometries and writing output...')
    for barcode in tqdm.tqdm(barcodes, ncols = 60):
        super_G_ = super_G.copy()
        for i, bit in enumerate(barcode.barcode, start = 1):
            if bit != 0:
                super_G_.add_edge(
                    f'bullvalene_{i}',
                    connectivity_map[bit]
                )
        mol = bv.graphs.graph_to_mol(
            super_G_
        )
        mol = bv.rdkit.generate_confs(
            mol,
            forcefield = args.forcefield,
            prune_rms_thresh = args.prune_rms_thresh,
            coord_map = coord_map,
            num_threads = args.num_threads
        )
        m_confs = min(
            mol.GetNumConformers(), args.m_confs
        )
        if m_confs > 0:
            d = Path(f'./{barcode}')
            if not d.is_dir():
                d.mkdir()
            for conf_idx in range(m_confs):
                out_d = d / f'./{barcode}_{conf_idx+1:03d}'
                if not out_d.is_dir():
                    out_d.mkdir()
                out_f = out_d / f'./{barcode}_{conf_idx+1:03d}'
                bv.io.mol_to_file(
                    out_f,
                    args.out_f_type,
                    mol,
                    conf_idx = conf_idx
                )
    print('...done!\n')

    datetime_ = datetime.datetime.now()
    print(f'finished @ {datetime_.strftime("%H:%M:%S (%Y-%m-%d)")}')

################################################################################
############################## PROGRAM STARTS HERE #############################
################################################################################

if __name__ == '__main__':
    main()

################################################################################
############################### PROGRAM ENDS HERE ##############################
################################################################################
