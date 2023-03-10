#!/usr/bin/env python3
import os
from collections import defaultdict

from tqdm import tqdm

from .cli import cli, parse_input_file
from .comet import mute_rdkit_warnings, run_comet, smiles_to_mol


def main() -> None:
    """
    Driver code for COMET.
    """
    mute_rdkit_warnings()
    args = cli()

    # Parse data from input file.
    _, file_extension = os.path.splitext(args.input)
    sep = { ".csv": ",", ".tsv": "\t" }[file_extension]
    data = defaultdict(lambda: defaultdict(list))
    for smiles_id, group_id, smiles in tqdm(parse_input_file(args.input, args.header, sep)):
        data[group_id]["mols"].append(smiles_to_mol(smiles))
        data[group_id]["mol_ids"].append(smiles_id)

    # Extract significant SMARTS and write them out to out file.
    with open(args.output, "w") as out_fo:
        out_fo.write("group_id\titem_id\tbit_index\tpval\tsmarts\n")

        for label, data_label, bit, adj_pval, smarts in tqdm(run_comet(data, **vars(args))):
            result_str = f"{label}\t{data_label}\t{bit}\t{round(adj_pval, 4)}\t{smarts}\n"
            out_fo.write(result_str)

    exit("done")


if __name__ == "__main__":
    main()
