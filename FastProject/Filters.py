# -*- coding: utf-8 -*-
"""Functions that are used to select genes

Functions here use a variety of criteria to reduce the number
of genes to a more manageable size - ideally extracting the
genes that are more biologically informative.

"""
from __future__ import absolute_import, print_function, division;

import numpy as np
from . import Utils
from .Global import FP_Output;


def apply_filters(data, threshold, nofilter, lean):
    """Applies filters to the ExpressionData object
    May remove rows.
    Populates the data.filters field

    Parameters
    ----------
    data : DataTypes.ExpressionData
        Data object to be filtered
    threshold : int
        Minimum number of samples gene must be detected in to pass
    nofilter : boolean
        if true, Only filter rows that have all identical values
    lean : boolean
        if true, skip extra filter methods (HDT, Fano)


    Returns
    -------
    data_out : DataTypes.ExpressionData
        Filtered data object

    """

    if(nofilter):
        filter_dict = {};
        data = filter_genes_novar(data);

        filter_dict.update({'No_Filter': set(data.row_labels)});
        data.filters = filter_dict;

        return data;

    else:
        filter_dict = {};
        data = filter_genes_threshold(data, threshold);

        filter_dict.update({
            'Threshold': set(data.row_labels),
        });

        if(not lean):
            for name, method in _filter_methods.items():
                FP_Output("Applying filter method:", name);

                mask = method(data);

                if(np.array(mask).sum() > 10):  # Only add these filters if they have enough genes
                    filter_dict.update({
                        name: set([data.row_labels[i] for i, x in enumerate(mask) if x]),
                    });

        data.filters = filter_dict;

        return data;


# This filter is run first, does not need to be registered at the bottom
def filter_genes_threshold(data, threshold):
    """Filters out genes that are at least active in <threshold>
    number of samples.
    """

    keep_indices = np.where((data > 0).sum(axis=1) > threshold)[0];

    return data.subset_genes(keep_indices);


# This filter is run when the --nofilter option is selected
# Needed for stability during other methods.
def filter_genes_novar(data):
    """
    Filters out genes with 0 variance across samples
    :param data: ExpressionData from DataTypes
    :return: Subset of data 0-variance rows removed
    """

    row_var = np.var(data, axis=1);
    keep_indices = np.where(row_var != 0)[0];

    return data.subset_genes(keep_indices);


# --------------------------------------------------------------------------- #
#                                                                             #
#                         Define Filter Methods Here                          #
#                                                                             #
# --------------------------------------------------------------------------- #

# Built-in Methods


def filter_genes_fano(data, num_mad=2):
    """
    Uses fano filtering on the genes.  Splits into quantiles first.
    Only retains genes that have fano factor <num_mad> median absolute
        deviations in each quantile.

    :param data: numpy.ndarray (Num_Genes x Num_Samples)
    :param num_mad: float
    :return: numpy.ndarray, bool, (Num_Genes)
        True for rows that pass the filter.  Otherwise, False;
    """
    mu = np.mean(data, axis=1);
    sigma = np.std(data, axis=1);

    # slice up mu and sigma into bins, by mu
    aa = np.argsort(mu);
    mu_sort = mu[aa];
    sigma_sort = sigma[aa];

    N_QUANTS = 30;
    m = mu_sort.size // N_QUANTS;

    gene_passes = np.zeros(data.shape[0]) == 1;

    for i in range(N_QUANTS):
        if(i == N_QUANTS - 1):
            rr = np.arange(i * m, mu_sort.size)
        else:
            rr = np.arange(i * m, (i + 1) * m);

        mu_quant = mu_sort[rr];
        mu_quant[mu_quant == 0] = 1;  # so we don't divide by zero later
        sigma_quant = sigma_sort[rr];
        fano_quant = sigma_quant**2 / mu_quant;
        mad_quant = np.median(np.abs(fano_quant - np.median(fano_quant)));
        gene_passes_quant = fano_quant > np.median(fano_quant) + num_mad * mad_quant;
        gene_passes_quant_i = np.nonzero(gene_passes_quant)[0];
        gene_passes_i = gene_passes_quant_i + i * m;
        gene_passes[gene_passes_i] = True;

    # gene_passes are in sorted mu order, need to revert back to original order
    original_ii = np.argsort(aa);
    gene_passes = gene_passes[original_ii];

    return gene_passes;


def filter_genes_hdt(data, p_val=0.05):
    """Filters out genes that pass the Hartigans Dip Test for bimodality
    with at least p < p_val"""
    # perform Hartigans dip test on the rest of the rows

    (dips, ps, xlows, xups) = Utils.HDT_Sig_batch(data, 1000);

    return ps <= p_val;


# Add additional filtering methods here

# Register filter methods
_filter_methods = dict();
_filter_methods["HDT"] = filter_genes_hdt;
_filter_methods["Fano"] = filter_genes_fano;
