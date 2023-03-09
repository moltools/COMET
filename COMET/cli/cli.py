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
    parser.add_argument("-i", "--input", type=str, required=True)
    return parser.parse_args()
