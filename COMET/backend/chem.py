import typing as ty 
import numpy as np 
from rdkit import RDLogger, Chem, DataStructs
from rdkit.Chem import AllChem 


def smiles_to_mol(smiles: str) -> Chem.Mol:
    """
    Construct molecular graph of SMILES string.

    Arguments
    ---------
    smiles (str) -- SMILES string representation of molecular graph.

    Returns
    -------
    Chem.Mol -- molecular graph.
    """
    return Chem.MolFromSmiles(smiles)


def mol_to_barcode(
    mol: Chem.Mol, 
    num_bits: int, 
    radius: int
) -> ty.Tuple[np.array, ty.Dict[int, ty.Tuple[int, int]]]:
    """
    Construct Morgan fingerprint barcode for molecular graph.

    Arguments
    ---------
    mol (Chem.Mol) -- molecular graph.
    num_bits (int) -- number of bits in Morgan fingerprint barcode to construct.
    radius (int) -- radius to use to construct Morgan fingerprint barcode.

    Returns
    -------
    np.array -- Morgan fingerprint barcode for molecular graph.
    dict -- bit setting info as { int: (atom_idx (int), radius (int)), ... }.
    """
    barcode = np.zeros((0,), dtype=int)
    bit_info = {}
    bit_vector = AllChem.GetMorganFignerprintAsBitVect(mol, radius, num_bits, bitInfo=bit_info)
    DataStructs.ConvertToNumpyArray(bit_vector, barcode)
    return barcode, bit_info


def mol_list_to_barcode(
    mols: ty.List[Chem.Mol],
    num_bits: int,
    radius: int
) -> np.array:
    """
    Construct barcode for a list of molecular graphs.
    
    Arguments
    ---------
    mols (Chem.Mol list) -- list of molecular graphs.
    num_bits (int) -- number of bits per Morgan fingerprint barcode to construct.
    radius (int) -- radius to use to construct per Morgan fingerprint barcode.

    Returns
    -------
    np.array -- barcode for list of molecular graphs.
    """
    return np.sum([mol_to_barcode(mol, num_bits, radius) for mol in mols], axis=0)