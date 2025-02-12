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

def mol_to_graph(
    mol: Chem.Mol,
    node_label_prefix: str = None
) -> nx.Graph:
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
            node_label_prefix + f'{i + 1}',
            element = atom.GetSymbol()
        )

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        G.add_edge(
            node_label_prefix + f'{i + 1}',
            node_label_prefix + f'{j + 1}',
            bond_type = bond.GetBondType()
        )

    return G

def graph_to_mol(
    G: nx.Graph,
    sanitize: bool = True
) -> Chem.Mol:
    """
    Converts a molecular graph (Network-X Graph object) into a molecule (RDKit
    Mol object).

    Args:
        G (nx.Graph): molecular graph (Network-X Graph object).
        sanitize (bool, optional): If `True`, the molecular structure is
            sanitized by RDKit using `Chem.SanitizeMol()` before the molecule is
            returned. Defaults to `True`.

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
            bond_type if bond_type is not None
                else Chem.rdchem.BondType.SINGLE
        )

    mol = mol.GetMol()

    if sanitize:
        Chem.SanitizeMol(mol)

    return mol

def smiles_to_graph(
    smiles: str,
    node_label_prefix: str = None
) -> nx.Graph:
    """
    Converts the SMILES string for a molecule into a molecular graph (Network-X
    Graph object).

    Args:
        smiles (str): SMILES string.
        node_label_prefix (str, optional): prefix to use for labelling nodes in
            the molecular graph.

    Returns:
        nx.Graph: molecular graph (Network-X object).
    """

    mol = Chem.MolFromSmiles(smiles)
    G = mol_to_graph(mol, node_label_prefix = node_label_prefix)
    
    return G

def graph_to_smiles(
    G: nx.Graph,
    sanitize: bool = True
) -> str:
    """
    Converts a molecular graph (Network-X Graph object) into a SMILES string
    for the molecule.

    Args:
        G (nx.Graph): molecular graph (Network-X object).
        sanitize (bool, optional): If `True`, the molecular structure is
            sanitized by RDKit using `Chem.SanitizeMol()` before the SMILES
            string is returned. Defaults to `True`.

    Returns:
        str: SMILES string.
    """

    mol = graph_to_mol(G, sanitize = sanitize)
    smiles = Chem.MolToSmiles(mol)

    return smiles

def compose_bullvalene_supergraph_from_smiles(
    bullvalene_smile: str = 'C12C=C4.C13C=C5.C23C=CC45',
    sub_smiles: list[str] = None
) -> nx.Graph:
    """
    Composes a disconnected molecular supergraph from SMILES strings that
    comprises the molecular graph of bullvalene and the molecular graphs of
    the supplied substituents; the nodes of the bullvalene molecular graph
    are prefixed with 'bullvalene_' and the nodes of the substituent molecular
    graphs are prefixed with 'sub[i]_' (where [i] is the number of the
    substituent, starting from 1) in the disconnected molecular supergraph.

    Args:
        bullvalene_smile (str, optional): SMILES string for bullvalene.
            Defaults to 'C12C=C4.C13C=C5.C23C=CC45'.
        sub_smiles (list[str], optional): SMILES strings for substituents.
            Defaults to None.

    Returns:
        nx.Graph: Disconnected molecular supergraph.
    """
    
    bullvalene_G = smiles_to_graph(
        bullvalene_smile, node_label_prefix = 'bullvalene_'
    )

    if sub_smiles is not None:
        sub_G = [
            smiles_to_graph(
                sub_smile, node_label_prefix = f'sub{i}_'
            ) for i, sub_smile in enumerate(sub_smiles)
        ]
    else:
        sub_G = []
        
    return nx.compose_all([bullvalene_G, *sub_G])
