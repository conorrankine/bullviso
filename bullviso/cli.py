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

import ast
from . import core
from . import utils
from argparse import (
    ArgumentParser,
    ArgumentDefaultsHelpFormatter,
    ArgumentTypeError,
    Namespace
)
from pathlib import Path

###############################################################################
################################## CONSTANTS ##################################
###############################################################################

SUPPORTED_CALCULATORS = ('mmff', 'uff', 'xtb')
SUPPORTED_OUTPUT_FILETYPES = ('xyz', 'sdf', 'gaussian', 'orca')

###############################################################################
############################## ARGUMENT PARSING ###############################
###############################################################################

def parse_args(
    argv: list[str] | None = None
) -> Namespace:
    """
    Parses command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments as an
        argparse.Namespace object that stores the arguments as attributes.
    """

    p = ArgumentParser(
        formatter_class = ArgumentDefaultsHelpFormatter
    )

    p.add_argument(
        'sub_smiles', type = str, nargs = '+',
        help = 'SMILES string representation for each unique substituent'
    )
    p.add_argument(
        '-n', '--n_subs', type = _int_or_list_of_ints, default = [1],
        help = 'number of each unique substituent to add'
    )
    p.add_argument(
        '-a', '--sub_attach_idx', type = _int_or_list_of_ints, default = [1],
        help =('atomic index of the substituent-bullvalene attachment point '
            'for each unique substituent')
    )
    p.add_argument(
        '-m', '--m_confs', type = int, default = 1,
        help = 'maximum number of conformational isomers to generate'
    )

    embedding_group = p.add_argument_group('conformer embedding options')
    embedding_group.add_argument(
        '-en', '--embed_n_confs', type = int, default = None,
        help = 'maximum number of conformational isomers to embed'
    )
    embedding_group.add_argument(
        '-er', '--embed_rmsd_threshold', type = float, default = 0.1,
        help = 'RMSD threshold (Angstroem) for deduplicating embeddings'
    )
    embedding_group.add_argument(
        '-et', '--embed_timeout', type = int, default = None,
        help = 'timout (seconds) for conformer embedding'
    )
    embedding_group.add_argument(
        '-es', '--embed_seed', type = int, default = None,
        help = 'random seed for conformer embedding'
    )

    optimisation_group = p.add_argument_group('conformer optimisation options')
    optimisation_group.add_argument(
        '-c', '--calculator_type', type = str, default = 'mmff',
        choices = SUPPORTED_CALCULATORS,
        help = 'calculator type for conformer optimisation'
    )
    optimisation_group.add_argument(
        '-it', '--max_iter', type = int, default = 600,
        help = 'maximum number of iterations for conformer optimisation'
    )

    cleanup_group = p.add_argument_group('conformer cleanup options')
    cleanup_group.add_argument(
        '-e', '--energy_threshold', type = float, default = 10.0,
        help = 'energy threshold (kcal/mol); '
    )
    cleanup_group.add_argument(
        '-r', '--rmsd_threshold', type = float, default = 0.5,
        help = 'RMSD threshold (Angstroem) for deduplicating optimised conformers'
    )

    system_group = p.add_argument_group('system options')
    system_group.add_argument(
        '-np', '--n_proc', type = int, default = 1,
        help = ('number of parallel processes for conformer embedding and '
            'optimisation')
    )

    output_group = p.add_argument_group('output options')
    output_group.add_argument(
        '-o', '--output_dir', type = Path, default = Path('.'),
        help = 'destination directory for outputting geometries'
    )
    output_group.add_argument(
        '-f', '--output_filetype', type = str, default = 'xyz',
        choices = SUPPORTED_OUTPUT_FILETYPES,
        help = 'filetype (e.g., .xyz, .sdf, etc.) for outputting geometries'
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
    input: int | list[int | list[int]]
) -> list[int | list[int]]:
    """
    Custom argument type for a command line argument to allow the user to
    pass an integer, list of integers, or list of integers with nested
    (sub)lists of integers as an input value.

    Args:
        input (int | list[int | list[int]]): Input value.

    Raises:
        ArgumentTypeError: If `input` is not an integer, list of integers, or
            list of integers with nested (sub)lists of integers.

    Returns:
        list[int | list[int]]: List of integers with or without nested
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
) -> list[int | list[int]]:
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
        list[int | list[int]]: List of integers and/or nested (sub)lists
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
        ValueError: If any of the following conditions are encountered:

            1. `args.sub_smiles`, `args.n_subs`, and `args.sub_attach_idx`
                are not of equal length;

            2. `args.n_subs` has a maximum depth greater than 1;

            3. `args.sub_attach_idx` has a maximum depth greater than 2;

            4. the sum of the integer elements in `args.n_subs` is greater
                than 9;

            5. there are greater than 9 elements in `args.sub_attach_idx`
                (counting elements in the outer list and nested (sub)lists,
                and accounting via multiplication for the number of each
                unique substituent given in `args.n_subs`);

            6. `args.n_subs` and/or `args.sub_attach_idx` contain null/zero-
                valued elements.
    """

    if not (
        len(args.sub_smiles) == len(args.n_subs) == len(args.sub_attach_idx)
    ):
        raise ValueError(
            '`args.sub_smiles`, `args.n_subs`, and `args.sub_attach_idx` '
            'should have the same length'
        )
    if utils.maxdepth(args.n_subs) > 1:
        raise ValueError(
            '`args.n_subs` should be a list of integers with a maximum depth '
            'no greater than 1'
        )
    if utils.maxdepth(args.sub_attach_idx) > 2:
        raise ValueError(
            '`args.sub_attach_idx` should be a list of integers or nested '
            '(sub)lists of integers with a maximum depth no greater than 2'
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

def main(
    argv: list[str] | None = None
) -> None:

    args = parse_args(argv)

    core.run_bullviso(args)
