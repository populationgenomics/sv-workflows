#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""

This script computes gene-level p-values using ACAT.

"""
import numpy as np
from scipy.stats import cauchy

import click

def CCT(pvals, weights=None):
    """
    Code adapted from the STAR package https://github.com/xihaoli/STAAR/blob/dc4f7e509f4fa2fb8594de48662bbd06a163108c/R/CCT.R wtih a modifitcaiton: when indiviudal p-value = 1, use minimum p-value
    #' An analytical p-value combination method using the Cauchy distribution
    #'
    #' The \code{CCT} function takes in a numeric vector of p-values, a numeric
    #' vector of non-negative weights, and return the aggregated p-value using Cauchy method.
    #' @param pvals a numeric vector of p-values, where each of the element is
    #' between 0 to 1, to be combined.
    #' @param weights a numeric vector of non-negative weights. If \code{NULL}, the
    #' equal weights are assumed.
    #' @return the aggregated p-value combining p-values from the vector \code{pvals}.
    #' @examples pvalues <- c(2e-02,4e-04,0.2,0.1,0.8)
    #' @examples CCT(pvals=pvalues)
    #' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
    #' with analytic p-value calculation under arbitrary dependency structures.
    #' \emph{Journal of the American Statistical Association 115}(529), 393-402.
    #' (\href{https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485}{pub})
    #' @export
    R code is implemented in python
    """

    #check if there is NA
    if np.isnan(pvals).sum() >0:
        raise ValueError("Cannot have NAs in the p-values!")

    # check if all p-values are between 0 and 1
    if ((pvals < 0).sum() + (pvals > 1).sum()) > 0:
        raise ValueError("All p-values must be between 0 and 1!")

    # check if there are p-values that are either exactly 0 or 1.
    is_zero = (pvals == 0).sum() >= 1
    is_one = (pvals == 1).sum() >= 1
    if is_zero:
        return 0
    if is_one:
        print('There are p-values exactly equal to 1!')
        return min(1, min(pvals) * len(pvals))

    # check the validity of weights (default: equal weights) and standardize them.
    if weights is None:
        weights = np.repeat(1/len(pvals), len(pvals))
    elif len(weights) != len(pvals):
        raise ValueError("The length of weights should be the same as that of the p-values!")
    elif (weights < 0).sum() > 0:
        raise ValueError("All the weights must be positive!")
    else:
        weights = weights / np.sum(weights)

    # check if there are very small non-zero p-values
    is_small = pvals < 1e-16
    if is_small.sum() == 0:
        cct_stat = np.sum(weights * np.tan((0.5 - pvals) * np.pi))
    else:
        cct_stat = np.sum((weights[is_small] / pvals[is_small]) / np.pi)
        cct_stat += np.sum(weights[~is_small] * np.tan((0.5 - pvals[~is_small]) * np.pi))

    # check if the test statistic is very large.
    if cct_stat > 1e+15:
        pval = (1 / cct_stat) / np.pi
    else:
        pval = 1 - cauchy.cdf(cct_stat)

    return pval

@click.option(
    '--input-dir',
    help='GCS path to the raw results of associaTR',
    type=str,
)
@click.option(
    '--cell-types',
    help='Name of the cell type, comma separated if multiple',
)
@click.option(
    '--chromosomes',
    help='Chromosome number eg 1, comma separated if multiple',
)
@click.command()
def main(input_dir, cell_types, chromosomes):
    for cell_type in cell_types.split(','):
        for chromosome in chromosomes.split(','):
            # read the raw results
            raw_results = pd.read_csv(f"{input_dir}/{cell_type}/chr{chromosome}.assoc", sep='\t')
            # compute gene-level p-values
            gene_pvals = raw_results.groupby('gene')['pval'].apply(CCT)
            # write the gene-level p-values to a file
            gene_pvals.to_csv(f"{input_dir}/{cell_type}/chr{chromosome}_gene_pvals.tsv", sep='\t', header=False)

