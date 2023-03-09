import typing as ty


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
