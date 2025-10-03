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

import click
import typer
import ast
from . import core
from . import utils
from pathlib import Path

# =============================================================================
#                                     APP
# =============================================================================

app = typer.Typer()

# =============================================================================
#                                   CLASSES
# =============================================================================

class NestedIntListParamType(click.ParamType):
    """
    Custom `click.ParamType` that parses a string literal into a (nested) list
    of integers.

    Examples:
        "1" -> [1]
        "[1,2]" -> [1,2]
        "[1,2,3]" -> [1,2,3]
        "[1,[2,3]] -> [1,[2,3]]
    """

    name = '(NESTED) INTEGER LIST'

    def __init__(
        self,
        max_depth: int = 1,
        require_positive: bool = False
    ) -> None:

        self.max_depth = max_depth
        self.require_positive = require_positive

    def convert(
        self,
        input_str: str,
        param: click.Parameter | None,
        ctx: click.Context | None
    ) -> list[int | list]:
        """
        Converts the input string literal into a (nested) list of integers. 

        Args:
            input_str (str): Input string literal to convert.
            param (click.Parameter, Optional): Click parameter (unused;
                required by the `click` API).
            ctx (click.Context, Optional): Click context (unused;
                required by the `click` API).

        Raises:
            typer.BadParameter: If the input string is invalid.

        Returns:
            list[int | list]: (Nested) list of integers.
        """

        try:
            result = ast.literal_eval(input_str.strip())
        except Exception as e:
            raise click.BadParameter(
                f'invalid Python literal syntax: {e}'
            )
        
        if isinstance(result, int):
            result = list([result])
               
        if not isinstance(result, list):
            raise click.BadParameter(
                f'input must be an integer or a (nested) list: got {result} '
                f'(type: {type(result).__name__})'
            )

        if not self._validate_input(
            result, require_positive = self.require_positive
        ):
            if not self.require_positive:
                raise click.BadParameter(
                    f'all elements of the input list and any (nested) '
                    f'sublists must be integers: got {result}'
                )
            else:
                raise click.BadParameter(
                    f'all elements of the input list and any (nested) '
                    f'sublists must be positive integers (> 0): got {result}'
                )
        
        max_depth = utils.maxdepth(result)
        if max_depth > self.max_depth:
            raise click.BadParameter(
                f'input list exceeds max allowed depth ({self.max_depth}): '
                f'got an input list with depth {max_depth}'
            )
        
        return result
        
    def _validate_input(
        self,
        input_: int | list,
        require_positive: bool = False
    ) -> bool:
        """
        Validates that the input is i) an integer, or ii) a list comprising
        only integers and (nested) lists of integers; optionally validates
        that the integers are positive (i.e., > 0).

        Args:
            input_ (int | list): Input to validate.
            require_positive (bool, optional): If True, all integers are
                required to be positive (i.e., > 0). Defaults to False. 

        Returns:
            bool: True if the input passes validation, else False.
        """

        if isinstance(input_, int):
            return input_ > 0 if require_positive else True
        if isinstance(input_, list):
            return all(
                self._validate_input(
                    item, require_positive = require_positive
                )
                for item in input_
            )
        return False

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

def _normalise_substituent_config_vars(
    sub_smiles: list[str],
    n_subs: list[int],
    sub_attach_idx: list[int | list[int]]
) -> tuple[list[int], list[int | list[int]]]:
    """
    Normalises substituent config variables (`n_subs` and `sub_attach_idx`).

    Args:
        sub_smiles (list[str]): SMILES strings for each unique substituent.
        n_subs (list[int]): Count for each unique substituent.
        sub_attach_idx (list[int  |  list[int]]): Attachment indices for each
            unique substituent.

    Raises:
        ValueError: If the lengths of the substituent config variables are
            incompatible, i.e., if `len(var) != 1` when `len(sub_smiles) > 1`,
            and `len(var) != len(sub_smiles)` in this case, for `var` in
            {`n_subs`, `sub_attach_idx`}.

    Returns:
        tuple[list[int], list[int | list[int]]]: Normalised substituent config
            variables (`n_subs` and `sub_attach_idx`).
    """

    n_smiles = len(sub_smiles)

    if n_smiles == 0:
        raise ValueError(
            '`sub_smiles` cannot be empty'
        )

    if len(n_subs) == 1 and n_smiles > 1:
        n_subs = n_subs * n_smiles
    elif len(n_subs) != n_smiles:
        raise ValueError(
            f'length of `n_subs` should be i) 1, or ii) equal to the length '
            f'of `sub_smiles` ({n_smiles}): got {len(n_subs)}'
        )

    if len(sub_attach_idx) == 1 and n_smiles > 1:
        sub_attach_idx = sub_attach_idx * n_smiles
    elif len(sub_attach_idx) != n_smiles:
        raise ValueError(
            f'length of `sub_attach_idx` should be i) 1, or ii) equal to the '
            f'length of `sub_smiles` ({n_smiles}): got {len(sub_attach_idx)}'
        )
    
    sub_attach_idx = [
        [attach_idx] if not isinstance(attach_idx, list) else attach_idx
        for attach_idx in sub_attach_idx
    ]

    sub_attach_idx = _to_zero_based_idxs(sub_attach_idx)

    return n_subs, sub_attach_idx

def _to_zero_based_idxs(
    idxs: list[list[int]]
) -> list[list[int]]:
    """
    Converts a 2D list of 1-based indices into a 2D list of 0-based indices.

    Args:
        idxs (list[list[int]]): (2D) Nested list of 1-based indices.

    Returns:
        list[list[int]]: (2D) Nested list of 0-based indices.
    """
    
    return [[idx - 1 for idx in sublist] for sublist in idxs]

@app.command()
def bullviso_cli(
    sub_smiles: list[str] = typer.Argument(
        ...,
        help = 'SMILES strings for each unique substituent type'
    ),
    n_subs: str = typer.Option(
        "1", '--n_subs', '-n',
        click_type = NestedIntListParamType(
            max_depth = 1, require_positive = True
        ),
        help = 'count for each unique substituent type'
    ),
    sub_attach_idx: str = typer.Option(
        "1", '--sub_attach_idx', '-a',
        click_type = NestedIntListParamType(
            max_depth = 2, require_positive = True
        ),
        help = 'attachment index for each unique substituent type'
    ),
    m_confs: int = typer.Option(
        1, '--m_confs', '-m',
        min = 1, max = None,
        help = 'maximum number of conformers to output'
    ),
    transition_state: bool = typer.Option(
        True,
        '--transition_state/--no_transition_state', '-ts/-no-ts',
        help = 'generate transition state geometries'
    ),
    embed_n_confs: int = typer.Option(
        None, '--embed_n_confs', '-en',
        min = 1, max = None,
        help = 'maximum number of conformers to embed'
    ),
    embed_rmsd_threshold: float = typer.Option(
        0.0, '--embed_rmsd_threshold', '-er',
        min = 0.0, max = None,
        help = 'RMSD threshold (Å) for deduplicating conformer embeddings'
    ),
    embed_timeout: int = typer.Option(
        None, '--embed_timeout', '-et',
        min = 1, max = None,
        help = 'timeout (seconds) for conformer embedding'
    ),
    embed_seed: int = typer.Option(
        None, '--embed_seed', '-es',
        help = 'random seed for conformer embedding'
    ),
    calculator_type: str = typer.Option(
        'mmff', '--calculator_type', '-c',
        click_type = click.Choice(
            ['mmff', 'uff', 'xtb'], case_sensitive = False
        ),
        help = 'calculator type for conformer optimisation'
    ),
    max_iter: int = typer.Option(
        600, '--max_iter', '-it',
        min = 1, max = None,
        help = 'maximum number of iterations for conformer optimisation'
    ),
    energy_threshold: float = typer.Option(
        10.0, '--energy_threshold', '-e',
        min = 0.0, max = None,
        help = 'energy threshold (kcal/mol) for discarding optimised conformers'
    ),
    rmsd_threshold: float = typer.Option(
        0.5, '--rmsd_threshold', '-r',
        min = 0.0, max = None,
        help = 'RMSD threshold (Å) for deduplicating optimised conformers'
    ),
    n_proc: int = typer.Option(
        1, '--n_proc', '-np',
        min = 1, max = None,
        help = 'number of parallel processes to launch'
    ),
    output_dir: Path = typer.Option(
        Path('.'), '--output_dir', '-o',
        exists = True, writable = True, dir_okay = True, file_okay = False,
        help = 'output directory path'
    ),
    output_filetype: str = typer.Option(
        'xyz', '--output_filetype', '-f',
        click_type = click.Choice(
            ['xyz', 'sdf', 'gaussian', 'orca'], case_sensitive = False
        ),
        help = 'output filetype'
    )
) -> None:
    
    try:

        n_subs, sub_attach_idx = _normalise_substituent_config_vars(
            sub_smiles, n_subs, sub_attach_idx
        )
    
        params = core.BullvisoParams(
            sub_smiles = sub_smiles,
            n_subs = n_subs,
            sub_attach_idx = sub_attach_idx,
            transition_state = transition_state,
            m_confs = m_confs,
            embed_n_confs = embed_n_confs,
            embed_rmsd_threshold = embed_rmsd_threshold,
            embed_timeout = embed_timeout,
            embed_seed = embed_seed,
            calculator_type = calculator_type,
            max_iter = max_iter,
            energy_threshold = energy_threshold,
            rmsd_threshold = rmsd_threshold,
            n_proc = n_proc,
            output_dir = output_dir,
            output_filetype = output_filetype,
        )
    
        core._bullviso(params)

    except Exception as e:
        err_prefix = typer.style('Error:', fg = typer.colors.RED)
        typer.echo(f'{err_prefix} {e}', err = True)
        raise typer.Exit(1)

# =============================================================================
#                                     EOF
# =============================================================================
