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

import networkx as nx
from rdkit import Chem

# =============================================================================
#                                  FUNCTIONS
# =============================================================================

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
            element = atom.GetSymbol(),
            chiral_tag = atom.GetChiralTag(),
            charge = atom.GetFormalCharge(),
            unpaired_electrons = atom.GetNumRadicalElectrons()
        )

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        G.add_edge(
            node_label_prefix + f'{i + 1}',
            node_label_prefix + f'{j + 1}',
            bond_type = bond.GetBondType(),
            bond_dir = bond.GetBondDir(),
            stereo_tag = bond.GetStereo()
        )

    return G

def graph_to_mol(
    G: nx.Graph,
    sanitize: bool = True,
    assign_stereochemistry: bool = True
) -> Chem.Mol:
    """
    Converts a molecular graph (Network-X Graph object) into a molecule (RDKit
    Mol object).

    Args:
        G (nx.Graph): molecular graph (Network-X Graph object).
        sanitize (bool, optional): If `True`, the molecule is sanitized by
            RDKit using `Chem.SanitizeMol()` before return. Defaults to `True`.
        assign_stereochemistry (bool, optional): If `True`, the stereochemistry
            of the molecule is assigned by RDKit using
            `Chem.AssignStereochemistry()` before return. Defaults to `True`.

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
        atom = mol.GetAtomWithIdx(
            node_to_atom_mapping[node_i]
        )
        atom.SetChiralTag(
            data.get('chiral_tag')
        )
        atom.SetFormalCharge(
            data.get('charge')
        )
        atom.SetNumRadicalElectrons(
            data.get('unpaired_electrons')
        )

    for node_i, node_j, data in G.edges(data = True):
        bond_type = data.get('bond_type')
        mol.AddBond(
            node_to_atom_mapping[node_i],
            node_to_atom_mapping[node_j],
            bond_type if bond_type is not None
                else Chem.rdchem.BondType.SINGLE
        )
        bond = mol.GetBondBetweenAtoms(
            node_to_atom_mapping[node_i],
            node_to_atom_mapping[node_j]
        )
        bond_dir = data.get('bond_dir')
        bond.SetBondDir(
            bond_dir if bond_dir is not None
                else Chem.rdchem.BondDir.NONE
        )
        stereo_tag = data.get('stereo_tag')
        bond.SetStereo(
            stereo_tag if stereo_tag is not None
                else Chem.rdchem.BondStereo.STEREONONE
        )

    mol = mol.GetMol()

    if sanitize:
        Chem.SanitizeMol(mol)

    if assign_stereochemistry:
        Chem.AssignStereochemistry(mol, cleanIt = True)

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

def set_atom_stereochemistry(
    G: nx.Graph,
    atom_stereo_map: dict
) -> nx.Graph:
    """
    Sets the R/S stereochemistry of the specified atom(s) in a molecular graph.

    Args:
        G (nx.Graph): Molecular graph.
        atom_stereo_map (dict): Dictionary mapping node labels (str) to chiral
            tags (str); valid chiral tags are {'ccw', 'cw', 'none'}.

    Raises:
        ValueError: If a node label is not a node in the molecular graph.
        ValueError: If a chiral tag is not one of {'ccw', 'cw', 'none'}.

    Returns:
        nx.Graph: Molecular graph with updated R/S (atom) stereochemistry.
    """
    
    stereo_map = {
        'cw': Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
        'ccw': Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
        'none': Chem.rdchem.ChiralType.CHI_UNSPECIFIED
    }

    for node_label, chiral_tag in atom_stereo_map.items():
        node_i = str(node_label)
        if node_i not in G.nodes:
            raise ValueError(
                f'node `{node_i}` is not a node in the graph'
            )
        if chiral_tag.lower() not in stereo_map:
            raise ValueError(
                f'`{chiral_tag.lower()}` is not a valid chiral tag: valid '
                f'chiral tags are {{{", ".join(stereo_map)}}}'
            )
        G.nodes[node_i]['chiral_tag'] = stereo_map[chiral_tag.lower()]

    return G

def set_bond_stereochemistry(
    G: nx.Graph,
    bond_stereo_map: dict
) -> nx.Graph:
    """
    Sets the E/Z stereochemistry of the specified bond(s) in a molecular graph.

    Args:
        G (nx.Graph): Molecular graph.
        bond_stereo_map (dict): Dictionary mapping pairs of node labels
            (tuple[str|int, str|int]) to stereo tags (str); valid stereo tags
            are {'cis', 'trans', 'none'}.

    Raises:
        ValueError: If a key in `bond_stereo_map` is not a tuple (or list) with
            two elements.
        ValueError: If the edge defined by a pair of node labels is not an edge
            in the molecular graph.
        ValueError: If a stereo tag is not one of {'cis', 'trans', 'none'}.

    Returns:
        nx.Graph: Molecular graph with updated E/Z (bond) stereochemistry.
    """
    
    stereo_map = {
        'cis': Chem.rdchem.BondStereo.STEREOCIS,
        'trans': Chem.rdchem.BondStereo.STEREOTRANS,
        'none': Chem.rdchem.BondStereo.STEREONONE        
    }

    for node_labels, stereo_tag in bond_stereo_map.items():
        if not isinstance(node_labels, (tuple, list)) or len(node_labels) != 2:
            raise ValueError(
                f'keys in `bond_stereo_map` should be tuples or lists of node '
                f'labels: `{node_labels}` is not a valid key'
            )
        node_i, node_j = map(str, node_labels)
        if not G.has_edge(node_i, node_j):
            raise ValueError(
                f'edge `{node_i}<->{node_j}` is not an edge in the graph'
            )
        if stereo_tag.lower() not in stereo_map:
            raise ValueError(
                f'`{stereo_tag.lower()}` is not a valid stereo tag: valid '
                f'stereo tags are {{{", ".join(stereo_map)}}}'
            )
        G[node_i][node_j]['stereo_tag'] = stereo_map[stereo_tag.lower()]

    return G

# =============================================================================
#                                     EOF
# =============================================================================