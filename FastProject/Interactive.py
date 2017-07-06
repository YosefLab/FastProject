"""
Functions to facilitate interactive use of the library
"""
import pandas as pd
from . import FileIO
from . import Transforms
from .DataTypes import ExpressionData
from .Projections import perform_weighted_PCA
from .Projections import permutation_wPCA as _permutation_wPCA
from .Signatures import (generate_random_sigs,
                         calculate_sig_scores as _calculate_sig_scores,
                         sigs_vs_projections as _sigs_vs_projections)

from .SigScoreMethods import SignatureScores


__all__ = ["estimate_weights", "weighted_PCA",
           "permutation_wPCA", "generate_random_sigs",
           "calculate_sig_scores"]


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
                     p_threshold=0.05, verbose=False, debug=False):
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
                                p_threshold, verbose, debug=debug)

    pc_data = pd.DataFrame(pc_data.T, index=data.columns)
    e_val = pd.Series(e_val, index=pc_data.columns)
    e_vec = pd.DataFrame(e_vec, index=data.index, columns=pc_data.columns)

    return pc_data, e_val, e_vec


def calculate_sig_scores(data, weights, signatures, method='weighted_avg',
                         zero_locations=None, min_signature_genes=1e99):
    """
    Interactive version of Signatures.calculate_sig_scores that
    uses DataFrames instead

    Parameters
    ----------
    data: pandas.DataFrame
        Data matrix to use to calculate signature scores
    weights: pandas.DataFrame
        Weight matrix to use for weighted sig-score evaluation
    signatures: list of Signature
        Signatures of which scores are computed
    method: str
        Which scoring method to use
    zero_locations: boolean numpy.ndarray
        Same size as data
        True where 'data' was originally zero
        Used for normalization methods that transform zeros to some other value
    min_signature_genes:  numeric
        Signatures that match less that `min_signature_genes` are discarded

    Returns
    -------
    pandas.DataFrame
        Size is NUM_SAMPLES x NUM_SIGNATURES
    """

    expressionMatrix = ExpressionData(
        data.values,
        list(data.index),
        list(data.columns)
    )

    expressionMatrix.weights = weights.values

    result_dict = _calculate_sig_scores(
        expressionMatrix, signatures,
        method=method,
        zero_locations=zero_locations,
        min_signature_genes=min_signature_genes)

    # Transform result_dict to dataframe

    all_series = [pd.Series(x.scores, index=x.sample_labels, name=name)
                  for name, x in result_dict.items()]

    all_scores = pd.concat(all_series, axis=1)

    num_genes = []
    names = []
    for name, sigscores in result_dict.items():
        names.append(name)
        num_genes.append(sigscores.numGenes)

    num_genes = pd.Series(num_genes, index=names)

    return all_scores, num_genes


def sigs_vs_projection(projection, sig_scores, num_genes,
                       random_sig_scores, random_num_genes,
                       numerical_sig_scores=None, factor_sig_scores=None,
                       NEIGHBORHOOD_SIZE=0.33):
    """
    Interactive wrapper for Signatures.sigs_vs_projections

    Uses dataframes for inputs/outputs

    projection: pandas.DataFrame
        Projection coordinates per sample.  Shape NUM_SAMPLES x 2

    sig_scores: pandas.DataFrame
        Signatures scores per sample.  Shape NUM_SAMPLES x NUM_SIGNATURES

    num_genes: pandas.Series
        Number of genes that matched for each signature.
        Size NUM_SIGNATURES

    random_sig_scores: pandas.DataFrame
        Background Signatures scores per sample.
        Shape NUM_SAMPLES x NUM_RANDOM_SIGNATURES

    random_num_genes: pandas.Series
        Number of genes that matched for each random signature.
        Size NUM_RANDOM_SIGNATURES

    numerical_sig_scores: pandas.DataFrame
        Precomputed numerical scores per sample.
        These use a permutation test instead of background signatures
        Shape NUM_SAMPLES x NUM_NUMERICAL_SIGNATURES

    factor_sig_scores: pandas.DataFrame
        Precomputed factor levels per sample.
        These use a permutation test instead of background signatures
        Shape NUM_SAMPLES x NUM_FACTOR_SIGNATURES

    NEIGHBORHOOD_SIZE: float
        Size in the projection to use for the kernel

    Returns
    -------
    consistences: pandas.DataFrame
        Consistency scores for each sig/projection pair
        Shape TOTAL_NUM_SIGNATURES x 1

    pvals: pandas.DataFrame
        P-vales for each sig/projection pair
        Shape TOTAL_NUM_SIGNATURES x 1

    """

    # Ensure all dataframes are aligned and same size
    assert len(projection.index & sig_scores.index) == len(projection.index)
    sig_scores = sig_scores.loc[projection.index]

    assert (len(projection.index & random_sig_scores.index) ==
            len(projection.index))
    random_sig_scores = random_sig_scores.loc[projection.index]

    if numerical_sig_scores is not None:
        assert (len(projection.index & numerical_sig_scores.index) ==
                len(projection.index))
        numerical_sig_scores = numerical_sig_scores.loc[projection.index]

    if factor_sig_scores is not None:
        assert (len(projection.index & factor_sig_scores.index) ==
                len(projection.index))
        factor_sig_scores = factor_sig_scores.loc[projection.index]

    # Convert values
    projections = {'proj_name': projection.values}
    sig_scores_dict = {}

    for name, column in sig_scores.iteritems():
        x = SignatureScores(column.values, name, column.index.tolist(),
                            isFactor=False, isPrecomputed=False,
                            numGenes=num_genes[name])
        sig_scores_dict[name] = x

    if numerical_sig_scores is not None:
        for name, column in numerical_sig_scores.iteritems():
            x = SignatureScores(column.values, name, column.index.tolist(),
                                isFactor=False, isPrecomputed=True,
                                numGenes=-1)
            sig_scores_dict[name] = x

    if factor_sig_scores is not None:
        for name, column in factor_sig_scores.iteritems():
            x = SignatureScores(column.values.tolist(), name,
                                column.index.tolist(),
                                isFactor=True, isPrecomputed=True,
                                numGenes=-1)
            sig_scores_dict[name] = x

    random_sig_scores_dict = {}
    for name, column in random_sig_scores.iteritems():
        x = SignatureScores(column.values, name, column.index.tolist(),
                            isFactor=False, isPrecomputed=False,
                            numGenes=random_num_genes[name])
        random_sig_scores_dict[name] = x

    result = _sigs_vs_projections(projections, sig_scores_dict,
                                  random_sig_scores_dict,
                                  NEIGHBORHOOD_SIZE=NEIGHBORHOOD_SIZE)

    row_labels, col_labels, sig_proj_matrix, sig_proj_matrix_p = result

    consistencies = pd.DataFrame(sig_proj_matrix, index=row_labels,
                                 columns=col_labels)

    pvals = pd.DataFrame(sig_proj_matrix_p, index=row_labels,
                         columns=col_labels)

    return consistencies, pvals
