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

import networkx as nx
from rdkit import Chem

###############################################################################
################################## FUNCTIONS ##################################
###############################################################################

def mol_to_graph(mol: Chem.Mol) -> nx.Graph:
    """
    Converts a molecule (RDKit Mol object) into a molecular graph (Network-X
    Graph object).

    Args:
        mol (Chem.Mol): molecule (RDKit Mol object).

    Returns:
        nx.Graph: molecular graph (Network-X Graph object).
    """

    G = nx.Graph()
    
    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            element = atom.GetSymbol()
        )

    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type = bond.GetBondType()
        )

    return G
