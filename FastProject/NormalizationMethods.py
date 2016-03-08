# -*- coding: utf-8 -*-
"""Methods to normalize data

Different options to normalize the data before the signature
values are calculated.

Determined using the command-line argument: --sig_norm_method

"""
from __future__ import absolute_import, print_function, division;
import numpy as np;
from scipy.stats import rankdata;

# Normalization methods


def no_normalization(data):
    """
    Does nothing. Returns original data
    """
    return data;


def col_normalization(data):
    """
    Perform z-normalization on all columns
    """
    data_mu = data.mean(axis=0, keepdims=True);
    data_sigma = data.std(axis=0, keepdims=True);
    data_sigma[data_sigma == 0] = 1.0;
    return (data - data_mu) / data_sigma;


def row_normalization(data):
    """
    Perform z-normalization on all rows
    """
    data_mu = data.mean(axis=1, keepdims=True);
    data_sigma = data.std(axis=1, keepdims=True);
    data_sigma[data_sigma == 0] = 1.0;
    return (data - data_mu) / data_sigma;


def row_and_col_normalization(data):
    """
    Normalize rows, then columns
    """
    data = row_normalization(data);
    data = col_normalization(data);
    return data;


def col_rank_normalization(data):
    """
    Create new version of data that has ranks (column-wise) instead of values
    """
    rdata = np.zeros(data.shape);
    for i in range(data.shape[1]):
        rdata[:, i] = rankdata(data[:, i], method="min");

    return rdata;
