#!/usr/bin/env python3
import argparse
import typing as ty
from dataclasses import dataclass

from rdkit import Chem


@dataclass
class Record:
    smiles_id: str
    group_id: str 
    smiles: str

    def as_csv(self) -> str:
        return f"{self.smiles_id},{self.group_id},{self.smiles}\n"


def check_smiles(smiles: str) -> bool:
    if not len(smiles):
        return False 
    else:
        try:
            if Chem.MolFromSmiles(smiles) is not None:
                return True
            else:
                return False
        except:
            return False


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Input data file from Eros")
    parser.add_argument("-o", type=str, required=True, help="Output file (csv)")
    return parser.parse_args()


def parse_input_file_Eros(input_file: str, header: bool = True) -> ty.Generator[Record, None, None]:
    fo = open(input_file, "r")
    header = fo.readline().strip().split("\t")
    for line in fo:
        line = line.strip().split("\t")
        named_items = dict(zip(header, line))
        smiles_id = f"{named_items['fragment_id']}_{named_items['fragment_index']}"
        group_id = f"{named_items['index_symbol']}_{named_items['methylated']}" 
        smiles = named_items["fragment_smiles"]
        if check_smiles(smiles):
            yield Record(smiles_id, group_id, smiles)   
    fo.close()


def main() -> None:
    args = cli()
    out_fo = open(args.o, "w")
    out_fo.write("smiles_id,group_id,smiles\n")
    total = 0
    for record in parse_input_file_Eros(args.i):
        out_fo.write(record.as_csv())
        total += 1  
    out_fo.close()
    exit(f"Parsed {total} records from {args.i}")


if __name__ == "__main__":
    main()
