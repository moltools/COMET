import typing as ty 
from copy import deepcopy

import numpy as np 
from rdkit import RDLogger, Chem, DataStructs
from rdkit.Chem import AllChem 


def mute_rdkit_warnings() -> None:
    """
    Mute RDKit warnings.

    Arguments
    ---------
    None

    Returns
    -------
    None
    """
    RDLogger.DisableLog("rdApp.*")


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
    mol = Chem.MolFromSmiles(smiles, )
    Chem.RemoveStereochemistry(mol)
    return mol


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
    bit_vector = AllChem.GetMorganFingerprintAsBitVect(mol, radius, num_bits, bitInfo=bit_info)
    DataStructs.ConvertToNumpyArray(bit_vector, barcode)
    return barcode, bit_info


def mol_list_to_barcode(
    mols: ty.List[Chem.Mol],
    num_bits: int,
    radius: int
) -> ty.Tuple[np.array, ty.List[ty.Dict[int, ty.Tuple[int, int]]]]:
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
    dict list -- bit setting info as { int: (atom_idx (int), radius (int)), ... } per mol.
    """ 
    barcodes, bit_infos = [], []

    for mol in mols:
        barcode, bit_info = mol_to_barcode(mol, num_bits, radius)
        barcodes.append(barcode)
        bit_infos.append(bit_info)

    return np.sum(barcodes, axis=0), bit_infos


def get_bit_smarts(
    mol: Chem.Mol,
    bit_info: ty.Dict[int, ty.Tuple[int, int]],
    bits: np.array
) -> ty.Any:    
    """
    Extract corresponding SMARTS from mol when given a bit index for barcode.

    Arguments
    ---------
    mol (Chem.Mol) -- molecular graph.
    bit_info (dict) -- bit setting info for mol as { int: (atom_idx (int), radius (int)), ... }.
    bits (int list) -- bits to retrieve SMARTS for, if present in bit_info. 

    Returns
    -------
    (bit index (int), SMARTS (str)) tuple -- resulting SMARTS for every bit.
    """
    for bit in np.nonzero(bits)[0]:
        if bit_setters := bit_info.get(bit, None):
            for atom_idx, radius in bit_setters:
                temp_mol = deepcopy(mol)
                env = Chem.FindAtomEnvironmentOfRadiusN(temp_mol, radius, atom_idx)
                atoms_to_use = set((atom_idx,))
                for bond in env:
                    bond_to_use = temp_mol.GetBondWithIdx(bond)
                    atoms_to_use.add(bond_to_use.GetBeginAtomIdx())
                    atoms_to_use.add(bond_to_use.GetEndAtomIdx())

                enlarged_env = set()
                for atom_idx in atoms_to_use:
                    atom = temp_mol.GetAtomWithIdx(atom_idx)    
                    for bond in atom.GetBonds():
                        bond_idx = bond.GetIdx()
                        if bond_idx not in env:
                            enlarged_env.add(bond_idx)

                enlarged_env = list(enlarged_env)
                enlarged_env += env 

                neighboring_atoms = []
                for atom_idx in atoms_to_use:
                    neighbors = temp_mol.GetAtomWithIdx(atom_idx).GetNeighbors()
                    for neighbor in neighbors:
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor_idx not in atoms_to_use:
                            neighboring_atoms.append(neighbor_idx)

                # Replace atom numbers with zero (i.e., highlight atoms as 
                # `any atom` a.k.a wildcard).
                for atom_idx in neighboring_atoms:
                    temp_mol.GetAtomWithIdx(atom_idx).SetAtomicNum(0)

                submol = Chem.PathToSubmol(temp_mol, enlarged_env)
                smarts = Chem.MolToSmarts(submol)
                smarts = smarts.replace("[#0]", "*")
                smarts = smarts.replace("[#0H]", "*")

                yield bit, smarts
