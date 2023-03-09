import argparse 


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
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Path to input tsv file with lines as `id\tgroup_id\tsmiles\n`")
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
    return parser.parse_args()
