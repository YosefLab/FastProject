"""
Functions to facilitate interactive use of the library
"""
import pandas as pd
from . import FileIO
from . import Transforms
from .DataTypes import ExpressionData
from .Projections import perform_weighted_PCA
from .Projections import permutation_wPCA as _permutation_wPCA


__all__ = ["estimate_weights", "weighted_PCA", "permutation_wPCA"]


def estimate_weights(data, housekeeping_genes=None):
    """
    Data is a pandas dataframe

    Returns: DataFrame of equivalent size with weights
    """

    expressionMatrix = ExpressionData(
        data.values,
        list(data.index),
        list(data.columns)
    )

    if housekeeping_genes is None:
        housekeeping_genes = FileIO.load_housekeeping_genes()

    (fit_func, params) = Transforms.create_false_neg_map(
        expressionMatrix, housekeeping_genes)

    weights = Transforms.compute_weights(fit_func, params, expressionMatrix)

    weights = pd.DataFrame(weights, index=expressionMatrix.row_labels,
                           columns=expressionMatrix.col_labels)

    return weights


def weighted_PCA(data, weights, max_components=200):
    """
    Wrapper around Projections.perform_weighted_PCA

    Paramters
    ---------

    data: pandas dataframe, genes x samples

    weights: pandas dataframe, genes x samples

    max_components: Number of components to compute
        - if max_components is smaller than either dimension,
          the results are truncated

    Returns
    -------

    pc_data : (Num_Samples x Num_Components) pandas dataframe
        Data transformed using PCA
        Same sample labels as input

    e_val:
        eigenvalues

    e_vec:
        eigenvectors
    """

    if (data.shape[0] != weights.shape[0] or
            data.shape[1] != weights.shape[1]):
        raise ValueError("Input data and weights must be same dimension")

    pc_data, e_val, e_vec = perform_weighted_PCA(data.values,
                                     weights.values, max_components)

    pc_data = pd.DataFrame(pc_data.T, index=data.columns)
    e_val = pd.Series(e_val, index=pc_data.columns)
    e_vec = pd.DataFrame(e_vec, index=data.index, columns=pc_data.columns)

    return pc_data, e_val, e_vec


def permutation_wPCA(data, weights, max_components=50,
                     p_threshold=0.05, verbose=False):
    """
    Computes weighted PCA on data.  Returns only significant components.

    After performing PCA on the data matrix, this method then uses
    a permutation procedure based on Buja A and Eyuboglu N (1992)
    to asses components for significance.

    :param data: pandas.Dataframe of shape Num_Features x Num_Samples
    :param weights: pandas.DataFrame of shape Num_Features x Num_Samples
    :param max_components: int, Max components to calculate
    :param p_threshold: float, P-value to cutoff components at
    :param verbose: bool

    Returns
    -------

    pc_data : (Num_Samples x Num_Components) pandas dataframe
        Data transformed using PCA
        Same sample labels as input

    e_val:
        eigenvalues

    e_vec:
        eigenvectors

    """

    pc_data, e_val, e_vec = _permutation_wPCA(data.values, weights.values, max_components,
                                p_threshold, verbose)

    pc_data = pd.DataFrame(pc_data.T, index=data.columns)
    e_val = pd.Series(e_val, index=pc_data.columns)
    e_vec = pd.DataFrame(e_vec, index=data.index, columns=pc_data.columns)

    return pc_data, e_val, e_vec
