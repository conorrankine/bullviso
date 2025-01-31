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

def graph_to_mol(G: nx.Graph) -> Chem.Mol:
    """
    Converts a molecular graph (Network-X Graph object) into a molecule (RDKit
    Mol object).

    Args:
        G (nx.Graph): molecular graph (Network-X Graph object).

    Returns:
        Chem.Mol: molecule (RDKit Mol object).
    """

    mol = Chem.RWMol()

    for node, data in G.nodes(data = True):
        atom = Chem.Atom(data.get('element'))
        mol.AddAtom(atom)

    for node_i, node_j, data in G.edges(data = True):
        bond_type = data.get('bond_type')
        mol.AddBond(node_i, node_j, bond_type)

    mol = mol.GetMol()

    return mol

def smiles_to_graph(smiles: str, sanitize = True) -> nx.Graph:
    """
    Converts the SMILES string for a molecule into a molecular graph (Network-X
    Graph object).

    Args:
        smiles (str): SMILES string.
        sanitize (bool, optional): toggles sanitization of the intermediate
            molecule generated from the SMILES string. Defaults to True.

    Returns:
        nx.Graph: molecular graph (Network-X object).
    """

    mol = Chem.MolFromSmiles(smiles, sanitize = sanitize)
    return mol_to_graph(mol)

def graph_to_smiles(G: nx.Graph, sanitize = True) -> str:
    """
    Converts a molecular graph (Network-X Graph object) into a SMILES string
    for the molecule.

    Args:
        G (nx.Graph): molecular graph (Network-X object).
        sanitize (bool, optional): toggles sanitization of the intermediate
            molecule generated from the molecular graph. Defaults to True.

    Returns:
        str: SMILES string.
    """

    mol = graph_to_mol(G)
    if sanitize:
        Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)
