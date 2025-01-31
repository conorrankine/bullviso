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

def mol_to_graph(mol: Chem.Mol, node_label_prefix: str = None) -> nx.Graph:
    """
    Converts a molecule (RDKit Mol object) into a molecular graph (Network-X
    Graph object).

    Args:
        mol (Chem.Mol): molecule (RDKit Mol object).
        node_label_prefix (str, optional): prefix to use for labelling nodes in
            the molecular graph.

    Returns:
        nx.Graph: molecular graph (Network-X Graph object).
    """

    G = nx.Graph()

    if node_label_prefix is None:
        node_label_prefix = ''
    
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        G.add_node(
            node_label_prefix + f'{i}',
            element = atom.GetSymbol()
        )

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        G.add_edge(
            node_label_prefix + f'{i}',
            node_label_prefix + f'{j}',
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

    node_to_atom_mapping = {}

    for node_i, data in G.nodes(data = True):
        element = data.get('element')
        node_to_atom_mapping[node_i] = mol.AddAtom(
            Chem.Atom(element)
        )

    for node_i, node_j, data in G.edges(data = True):
        bond_type = data.get('bond_type')
        mol.AddBond(
            node_to_atom_mapping[node_i],
            node_to_atom_mapping[node_j],
            bond_type
        )

    mol = mol.GetMol()

    return mol

def smiles_to_graph(
    smiles: str,
    sanitize: bool = True,
    node_label_prefix: str = None
) -> nx.Graph:
    """
    Converts the SMILES string for a molecule into a molecular graph (Network-X
    Graph object).

    Args:
        smiles (str): SMILES string.
        sanitize (bool, optional): toggles sanitization of the intermediate
            molecule generated from the SMILES string. Defaults to True.
        node_label_prefix (str, optional): prefix to use for labelling nodes in
            the molecular graph.

    Returns:
        nx.Graph: molecular graph (Network-X object).
    """

    mol = Chem.MolFromSmiles(smiles, sanitize = sanitize)
    G = mol_to_graph(mol, node_label_prefix = node_label_prefix)
    
    return G

def graph_to_smiles(G: nx.Graph, sanitize: bool = True) -> str:
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
    smiles = Chem.MolToSmiles(mol)

    return smiles
