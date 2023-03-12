import typing as ty 
from enum import Enum, auto

import numpy as np
from scipy.stats import hypergeom


class MultipleTestingCorrection(Enum):
    Bonferroni = auto()

    @staticmethod
    def has_item(mtc: str) -> ty.List[str]:
        """
        Check if MultipleTestingCorrection item exists.

        Arguments
        ---------
        mtc (str) -- name of multiple testing correction.

        Returns
        -------
        bool -- True if name exists, otherwise False.
        """
        return mtc in [x.name.lower() for x in MultipleTestingCorrection]


def bit_enrichment(
    group_barcodes: np.array, 
    group_sizes: ty.List[int],
    alpha: float,
    mtc: str
) -> np.array:
    """
    Perform bit enrichment on a set of group barcodes.

    Arguments
    ---------
    group_barcodes (np.array of shape (N, M)) -- set of group barcodes.
    group_sizes (int list of length N) -- group sizes for all barcodes.
    alpha (float) -- significance threshold.
    mtc (str) -- multiple testing correction type.

    Returns
    -------
    (adjusted p-values (np.array), significant bits (np.array)) tuple
    """
    assert (len(group_barcodes.shape) == 2), (
        f"Group barcodes array should be of shape (N, M): found `{group_barcodes.shape}`")
    assert (group_barcodes.shape[0] ==  len(group_sizes)), (
        f"Number of group sizes should correspond to number of group barcodes: `{group_barcodes.shape[0]}` != `{len(group_sizes)}`")
    assert (mtc is None or MultipleTestingCorrection.has_item(mtc)), (
        f"Got unknown multiple testing correction type: `{mtc}`")

    # Count total number of bit occurences per bit.
    ns = {i: sum(group_barcodes[:, i]) for i in range(group_barcodes.shape[1])}

    # Calculate population size.
    M = sum(group_sizes)

    # Create storage array for storing all calculated p-values in.
    pvals = np.zeros(group_barcodes.shape)

    for barcode_idx in range(group_barcodes.shape[0]):
        for bit_idx in range(group_barcodes.shape[1]):
            # Retrieve number of successes in group/sample.
            k = group_barcodes[barcode_idx, bit_idx]

            # Retrieve number of successes in population.
            n = ns[bit_idx]

            # Retrieve size of group/sample.
            N = group_sizes[barcode_idx]
            
            # Calculate and store pval.
            pvals[barcode_idx, bit_idx] = hypergeom.sf(k - 1, M, n, N) 

    # Create storage array for storing all adjusted p-values in.
    adj_pvals = np.zeros(pvals.shape)

    # Create storage array for storing all significance decisions in.
    significant_bits = np.zeros(pvals.shape)

    for barcode_idx in range(pvals.shape[0]):
        for bit_idx in range(pvals.shape[1]):
            pval = pvals[barcode_idx, bit_idx]
            
            # Apply multiple testing correction.
            if mtc == "bonferroni":
                adj_pval = min(pval * pvals.size, 1.0)
            else:
                adj_pval = pval
                
            # Decide on significance.
            significant_bits[barcode_idx, bit_idx] = int(adj_pval <= alpha) 

            # Adjust p-value.
            adj_pvals[barcode_idx, bit_idx] = adj_pval

    return adj_pvals, significant_bits 
