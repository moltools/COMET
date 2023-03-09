#!/usr/bin/env python3
import os
from collections import defaultdict

import numpy as np

from .cli import cli, parse_input_file
from .backend import mute_rdkit_warnings, smiles_to_mol, mol_list_to_barcode, bit_enrichment, get_bit_smarts


def main() -> None:
    """
    Driver code for COMET.
    """
    mute_rdkit_warnings()
    args = cli()

    # Parse data from input file.
    _, file_extension = os.path.splitext(args.input)
    sep = { ".csv": ",", ".tsv": "\t" }[file_extension]
    data, data_labels = defaultdict(list), defaultdict(list)
    for smiles_id, group_id, smiles in parse_input_file(args.input, args.header, sep):
        data[group_id].append(smiles_to_mol(smiles))
        data_labels[group_id].append(smiles_id)

    # Create group barcodes and save bit setting info per mol.
    group_labels = list(data.keys())
    group_sizes = [len(data[label]) for label in group_labels]
    group_barcodes, group_bit_infos = zip(*[
        mol_list_to_barcode(data[label], args.num_bits, args.radius) 
        for label in group_labels])
    group_barcodes = np.array(group_barcodes)

    # Get significant bits.
    adj_pvals, significant_bits = bit_enrichment(group_barcodes, group_sizes, args.alpha, args.mtc)

    # Extract significant SMARTS and write them out to out file.
    with open(args.output, "w") as out_fo:
        out_fo.write("group_id\titem_id\tbit_index\tpval\tsmarts\n")

        for label_idx, label in enumerate(group_labels):
            significant_bits_for_group = significant_bits[label_idx, :]
            
            for mol_idx, mol in enumerate(data[label]):
                bit_info = group_bit_infos[label_idx][mol_idx]
                data_label = data_labels[label][mol_idx]
                found_for_datum = set()

                for bit, smarts in get_bit_smarts(mol, bit_info, significant_bits_for_group):
                    adj_pval = adj_pvals[label_idx, bit]
                    found = frozenset((bit, smarts))

                    if found not in found_for_datum:
                        found_for_datum.add(found)
                        result_str = f"{label}\t{data_label}\t{bit}\t{round(adj_pval, 4)}\t{smarts}\n"
                        out_fo.write(result_str)

    exit("Done.")


if __name__ == "__main__":
    main()
