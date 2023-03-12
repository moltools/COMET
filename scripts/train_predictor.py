#!/usr/bin/env python3
from comet import smiles_to_mol


def main() -> None:
    print(smiles_to_mol("CCC"))
    exit(0)


if __name__ == "__main__":
    main()
