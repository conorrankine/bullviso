"""
PYBULL
Copyright (C) 2022  Conor D. Rankine

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

from rdkit import Chem
from rdkit.Chem import AllChem

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def functionalise(
    mol: Chem.Mol,
    func_group: Chem.Mol,
    confcode: tuple
) -> Chem.Mol:
    """
    Functionalises a Chem.Mol molecule `mol` with the Chem.Mol functional group
    `func_group` according to a configuration code tuple `confcode`; N 
    (0 <= N <= `mol.GetNumAtoms()`) instances of `func_group` are attached to
    `mol` via single bonds at the positions indicated in `confcode`.

    Args:
        mol (Chem.Mol): A molecule to functionalise.
        func_group (Chem.Mol): A functional group.
        confcode (tuple): A configuration code encoding the functionalisation
            of `mol` as a tuple of binary (eg. 1/0) elements with length
            `mol.GetNumAtoms()`; each binary element indicates whether
            or not that (functionalisable) site in `mol` is functionalised.

    Returns:
        Chem.Mol: A functionalised molecule.
    """

    func_sites = tuple(i for i in range(len(confcode)) if confcode[i])

    func_mol = Chem.Mol(mol)
    
    for i, func_site in enumerate(func_sites):
        func_mol = Chem.EditableMol(
            Chem.CombineMols(func_mol, func_group)
        )
        func_mol.AddBond(
            func_site, 10 + (i * func_group.GetNumAtoms()),
            order = Chem.rdchem.BondType.SINGLE
        )
        func_mol = func_mol.GetMol()

    Chem.SanitizeMol(func_mol)

    func_mol = Chem.AddHs(func_mol)

    return func_mol
