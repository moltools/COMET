#!/usr/bin/env python3
import argparse 

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def cli() -> argparse.Namespace:
    """
    Command line interface for script.

    Arguments
    ---------
    None

    Returns
    -------
    argparse.Namespace -- command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="SMARTS string to draw.")
    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="Path to output svg file.")
    return parser.parse_args()


def main() -> None:
    """
    Driver code for script.
    """
    args = cli()

    try:    
        mol = Chem.MolFromSmarts(args.input)
        Chem.SanitizeMol(mol)
        drawing = rdMolDraw2D.MolDraw2DSVG(*(300, 300))
        options = drawing.drawOptions()
        options.useBWAtomPalette()
        drawing.DrawMolecule(mol)
        drawing.FinishDrawing()
        svg_str = drawing.GetDrawingText().replace("svg:", "")
        
        with open(args.output, "w") as out_fo:
            out_fo.write(svg_str)
    
    except:
        exit(f"Could not draw SMARTS for `{args.input}`.")

    exit(f"Success. Written out svg to `{args.output}`.")


if __name__ == "__main__":
    main()
