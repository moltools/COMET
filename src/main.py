#!/usr/bin/env python3
import argparse 
import os
import typing as ty 
from collections import defaultdict

from tqdm import tqdm

from comet import mute_rdkit_warnings, run_comet, smiles_to_mol


def cli() -> argparse.Namespace:
    """
    Command line interface for COMET.

    Arguments
    ---------
    None

    Returns
    -------
    argparse.Namespace -- parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        prog="comet",
        description="COMET: COmparative Molecular Enrichment Tool")

    # Required arguments.
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Path to input file with lines as `smiles_id\\tgroup_id\\tsmiles\\n`.")
    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="Path to output tsv file.")
    
    # Optional arguments.
    parser.add_argument(
        "-a", "--alpha", type=float, required=False, default=0.05,
        help="Significance threshold.")
    parser.add_argument(
        "-n", "--num_bits", type=int, required=False, default=2048,
        help="Number of bits in molecular barcode.")
    parser.add_argument(
        "-r", "--radius", type=int, required=False, default=3,
        help="Node radius to use for constructing molecular barcode.")
    parser.add_argument(
        "-c", "--mtc", type=str, required=False, choices=["bonferroni"], default=None,
        help="Multiple testing correction to use.")
    
    # Flags.
    parser.add_argument(
        "--header", action="store_true",
        help="Flag to indicate that input file contains header.")
    
    return parser.parse_args()


def parse_input_file(
    file_path: str, 
    header: bool, 
    sep: str
) -> ty.Generator[ty.Tuple[str, str, str], None, None]:
    """
    Parse input file.

    Arguments
    ---------
    file_path (str) -- path to input file.
    header (bool) -- ignores first line if True, otherwise parses first line.
    sep (str) -- seperator to split line in input file on.

    Returns
    -------
    (smiles_id (str), group_id( str), smiles (str) list -- yields parsed lines.
    """
    with open(file_path, "r") as fo:
        if header:
            fo.readline()

        for line in fo:
            smiles_id, group_id, smiles = line.strip().split(sep)
            yield smiles_id, group_id, smiles


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
