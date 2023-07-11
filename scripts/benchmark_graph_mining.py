#!/usr/bin/env python3
import argparse 
import typing as ty
from collections import defaultdict
from itertools import chain, combinations

from tqdm import tqdm 

import numpy as np

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem

from comet import bit_enrichment, mol_to_barcode, get_bit_smarts


Graph = ty.Dict[int, ty.List[int]]


# argparse with two subparsers that determine mode:
def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments.

    Arguments
    ---------
    None

    Returns
    -------
    argparse.Namespace -- parsed command line arguments.
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode")
    mode1_parser = subparsers.add_parser("naive_mode", help="naive mode")
    mode2_parser = subparsers.add_parser("fingerprint_mode", help="fingerprint mode")
    return parser.parse_args()


def smi_to_mol(smi: str) -> Chem.Mol:
    """
    Convert a SMILES string to a molecule.

    Args:
        smi (str): The SMILES string to convert.
    
    Returns:
        Chem.Mol: The molecule.
    """
    mol = Chem.MolFromSmiles(smi, sanitize=False)

    # Sanitize the molecule, but don't set aromaticity.
    rdmolops.SanitizeMol(mol, rdmolops.SanitizeFlags.SANITIZE_ALL ^ rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY)

    return mol


def mol_to_graph(mol: Chem.Mol) -> Graph:
    """
    Convert a molecule to a graph representation.

    Args:
        mol (Chem.Mol): The molecule to convert.
    
    Returns:
        Graph: The graph representation of the molecule.
    """
    # Initialize the graph as a dictionary of lists.
    graph = defaultdict(list)

    # Add all atoms to the graph.
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()
        graph[atom1].append(atom2)
        graph[atom2].append(atom1)

    return graph


def smi_to_graph(smi: str) -> Graph:
    """
    Convert a SMILES string to a graph representation.

    Args:
        smi (str): The SMILES string to convert.
    
    Returns:
        Graph: The graph representation of the molecule.
    """
    mol = smi_to_mol(smi)
    graph = mol_to_graph(mol)
    return graph


def dfs(graph: Graph, start: int, visited: ty.Optional[ty.Set[int]] = None) -> ty.Set[int]:
    """
    Perform a depth-first search (DFS) starting from a given node.

    Args:
        graph (Graph): The graph to perform DFS on.
        start (int): The node to start DFS from.
        visited (ty.Optional[ty.Set[int]], optional): The set of visited nodes. Defaults to None.

    Returns:
        ty.Set[int]: The set of visited nodes.
    """
    # If this is the first call, initialize visited to an empty set.
    if visited is None:
        visited = set()
    
    # Add the starting node to the visited set.
    visited.add(start)
    
    # For every node that is connected to the starting node and not yet visited, perform DFS.
    for next in set(graph[start]) - visited:
        dfs(graph, next, visited)
    
    # Return the set of visited nodes.
    return visited


def connected_components(graph: Graph) -> ty.Generator[ty.Set[int], None, None]:
    """
    Find all connected components in a graph.
    
    Args:
        graph (Graph): The graph to find connected components in.
        
    Returns:
        ty.Generator[ty.Set[int], None, None]: A generator of connected components.
    """
    # Initialize visited to an empty set.
    visited = set()
    
    # For every node in the graph.
    for node in graph:
        # If the node hasn't been visited yet, perform DFS from this node and yield the result.
        if node not in visited:
            yield dfs(graph, node, visited)


def powerset(iterable: ty.Iterable) -> ty.Iterable[ty.Iterable]:
    """
    Generate all possible subsets of an iterable (including the empty set).

    Args:
        iterable (ty.Iterable): The iterable to generate subsets from.

    Returns:
        ty.Iterable: The powerset of the iterable.
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


def is_connected(graph: Graph) -> bool:
    """
    Check if a graph is connected.  
    
    Args:
        graph (Graph): The graph to check.
    
    Returns:
        bool: True if the graph is connected, False otherwise.
    """
    # If the graph has more than one connected component, it is not connected.
    return len(list(connected_components(graph))) == 1


def all_subgraphs(graph: Graph) -> ty.Generator[Graph, None, None]:
    """
    Generate all subgraphs of a graph that are connected.
    
    Args:
        graph (Graph): The graph to generate subgraphs from.
    
    Returns:
        ty.Generator[Graph, None, None]: A generator of all subgraphs of the graph that are connected.
    """
    # Find all connected components in the graph.
    components = list(connected_components(graph))
    
    # For each connected component.
    for component in components:
        # Generate all possible subsets of the component (i.e., all subgraphs).
        subgraphs = list(powerset(component))
        
        # For each subgraph.
        for subgraph in subgraphs:
            # If the subgraph isn't empty, yield it.
            if subgraph:
                # Generate the subgraph as a dictionary.
                subgraph = {node: [n for n in graph[node] if n in subgraph] for node in subgraph}

                # If the subgraph is connected, yield it (final check to prevent disconnected subgraphs).
                if is_connected(subgraph):
                    yield subgraph


def get_submol(mol: Chem.Mol, atoms: ty.List[int]) -> Chem.Mol:
    """
    Get a submolecule from a molecule. Replace atoms that are not in the submolecule but are still connected to it with dummy atoms.

    Args:
        mol (Chem.Mol): The molecule to get the submolecule from.
        atoms (ty.List[int]): The atoms to include in the submolecule.

    Returns:
        Chem.Mol: The submolecule.
    """
    # Create a submolecule from the molecule.
    submol = Chem.RWMol(mol)

    # List to keep track of the dummy atoms.
    dummy_atoms = []
    
    # Iterate over all atoms in the molecule. Start batch edit.
    submol.BeginBatchEdit()
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        
        # If the atom is not part of the given atom list.
        if idx not in atoms:
            # If the atom is connected to an atom that is in the atom list, mark it as a dummy atom.
            if any(neighbor.GetIdx() in atoms for neighbor in atom.GetNeighbors()):
                dummy_atoms.append(idx)
            # Else, remove the atom from the submolecule.
            else:
                submol.RemoveAtom(idx)
                
    # Replace atoms that are not in the atom list but are connected to the submolecule with dummy atoms.
    for idx in dummy_atoms:
        atom = submol.GetAtomWithIdx(idx)
        atom.SetAtomicNum(0)  # 0 is the atomic number for a dummy atom.

    # Remove bonds between wildcards.
    for bond in submol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 0 and bond.GetEndAtom().GetAtomicNum() == 0:
            submol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # End batch edit.
    submol.CommitBatchEdit()
    
    # Return submol.
    return submol.GetMol()


def mol_to_fp(mol: Chem.Mol, radius: int, nbits: int) -> np.ndarray:
    """
    Convert a molecule to a Morgan fingerprint.

    Args:
        mol (Chem.Mol): Molecule.
        radius (int): Morgan fingerprint radius.
        nbits (int): Number of bits in the fingerprint.

    Returns:
        np.ndarray: Morgan fingerprint.
    """
    arr = np.zeros((0,), dtype=int)
    morgan_bits = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    DataStructs.ConvertToNumpyArray(morgan_bits, arr)
    return arr


def mol_to_encoding(mol: Chem.Mol, radius: int, nbits: int) -> int:
    """
    Convert a molecule to an encoding.  

    Args:
        mol (Chem.Mol): Molecule.
        radius (int): Morgan fingerprint radius.
        nbits (int): Number of bits in the fingerprint.

    Returns:
        int: Encoding.
    """
    fp = mol_to_fp(mol, radius, nbits)
    return hash(fp.data.tobytes())


def encode_smi(smi: str, radius: int = 2, nbits: int = 2048) -> int:
    """
    Encode a SMILES string.

    Args:
        smi (str): SMILES string.
        radius (int, optional): Morgan fingerprint radius. Defaults to 2.
        nbits (int, optional): Number of bits in the fingerprint. Defaults to 2048.

    Returns:
        int: Encoding.
    """
    mol = smi_to_mol(smi)
    return mol_to_encoding(mol, radius, nbits)


def main() -> None:
    RDLogger.DisableLog("rdApp.*")

    args = parse_args()
    smi = r"CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"
    # smi = r"CCC"

    submol_encs, submols = set(), set()
    success, failed = 0, 0

    if args.mode == "naive_mode":
        mol = smi_to_mol(smi)
        graph = mol_to_graph(mol)

        for subgraph in tqdm(all_subgraphs(graph), leave=True):
            try:
                atom_inds = list(subgraph.keys())
                submol = get_submol(mol, atom_inds)
                submol_smi = Chem.MolToSmiles(submol)
                enc = encode_smi(submol_smi)
                if not enc in submol_encs:
                    submols.add(submol_smi)
                    submol_encs.add(enc)
                success += 1
            except:
                failed += 1

    elif args.mode == "fingerprint_mode":
        mol = Chem.MolFromSmiles(smi)
        bits, bit_info = mol_to_barcode(mol, num_bits=2048, radius=2)
     
        for _, smarts in tqdm(get_bit_smarts(mol, bit_info, bits), leave=True):
            try:
                enc = encode_smi(smarts)
                if not enc in submol_encs:
                    submols.add(smarts)
                    submol_encs.add(enc)
                success += 1
            except:
                try:
                    # Replace aromatic bonds with wildcard bonds when mol cannot
                    # be created due to failing kekulization:
                    smarts = smarts.replace(":", "~")
                    enc = encode_smi(smarts)
                    if not enc in submol_encs:
                        submols.add(smarts)
                        submol_encs.add(enc)
                    success += 1
                except:
                    failed += 1
    
    else:
        raise ValueError(f"Unknown mode: {args.mode}")
    
    print(f"Success: {success}/{success + failed} ({success / (success + failed):.2%})")
    print(f"Unique submols: {len(submols)}")

    exit(0)


if __name__ == "__main__":
    main()
