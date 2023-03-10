import typing as ty 

import numpy as np 

from .chem import mol_list_to_barcode, get_bit_smarts
from .stats import bit_enrichment 


def run_comet(
    data: ty.Dict[str, ty.Dict[str, ty.List]],
    num_bits: int, 
    radius: int, 
    alpha: float, 
    mtc: str,
    **kwargs
) -> ty.Generator[ty.Tuple[str, str, int, float, str], None, None]:
    """
    Run COMET.

    Arguments
    ---------
    data (dict of dicts) -- as { group_id (str): { mols: Chem.Mol list, mol_ids: str list } } 
        with mol_ids as optional.
    num_bits -- number of bits in barcodes.
    radius -- radius to use to construct barcodes.
    alpha -- significance threshold.
    mtc -- multiple testing correction type.

    Returns
    -------
    generates tuple of (group_id (str), mol_id (str) or None, bit (int), p-value (float), SMARTS (str))
    """
    # Create group barcodes and save bit setting info per mol.
    group_labels = list(data.keys())
    group_sizes = [len(data[label]["mols"]) for label in group_labels]
    group_barcodes, group_bit_infos = zip(*[
        mol_list_to_barcode(data[label]["mols"], num_bits, radius) 
        for label in group_labels])
    group_barcodes = np.array(group_barcodes)

    # Get significant bits.
    adj_pvals, significant_bits = bit_enrichment(group_barcodes, group_sizes, alpha, mtc)

    # Extract significant SMARTS and write them out to out file.
    for label_idx, label in enumerate(group_labels):
        significant_bits_for_group = significant_bits[label_idx, :]
        
        for mol_idx, mol in enumerate(data[label]["mols"]):
            bit_info = group_bit_infos[label_idx][mol_idx]
            found_for_datum = set()

            if "mol_ids" in data[label].keys():
                data_label = data[label]["mol_ids"][mol_idx]
            else:
                data_label = None

            for bit, smarts in get_bit_smarts(mol, bit_info, significant_bits_for_group):
                adj_pval = adj_pvals[label_idx, bit]
                found = frozenset((bit, smarts))

                if found not in found_for_datum:
                    found_for_datum.add(found)
                    yield label, data_label, bit, adj_pval, smarts
