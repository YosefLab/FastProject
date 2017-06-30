# -*- coding: utf-8 -*-
"""Functions for generating projections

This module handles the generation of lower-dimensional
projections from the higher-dimensional data objects

"""
from __future__ import absolute_import, print_function, division;
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA
from sklearn.decomposition import KernelPCA
from sklearn.manifold import TSNE;
from sklearn.manifold import Isomap
from sklearn.manifold import MDS
from sklearn.manifold import SpectralEmbedding
from sklearn.cluster import MiniBatchKMeans
from scipy.stats import norm;
import scipy.stats;

from .Utils import ProgressBar;
from .DataTypes import PCData;
from .Global import FP_Output, RANDOM_SEED;
from . import _tsne_fix

import numpy as np;


def generate_projections(data, filter_name=None, input_projections=None, lean=False):
    """
    Projects data into 2 dimensions using a variety of linear and non-linear methods

    Parameters
    ----------
    data : (Num_Features x Num_Samples) numpy.ndarray
        Matrix containing data to project into 2 dimensions
    filter_name : String
        Filter to apply to the signatures.  Should match a filter in data.filters
        If not specified, no filtering is applied.
    input_projections : dict
        Keys are of type str, representing projection names
        Values are of type 2xN pandas.DataFrame with column names matching
        sample names in `expressionMatrix`

    Returns
    -------
    projections : dict(string, (2 x Num_Samples) numpy.ndarray)
        dictionary mapping the projection type (e.g. "tSNE") to an array containing
        the two-dimensional coordinates for each sample in the projection.

    PC_Data : (Num_Components x Num_Samples) numpy.ndarray
        Weighted PCA of the original data object

    """
    if(input_projections is None):
        input_projections = {}

    IS_PCDATA = type(data) is PCData;

    method_list = get_projection_methods(lean, IS_PCDATA)

    pbar = ProgressBar(1 + len(method_list));

    projections = dict();

    proj_data = data.projection_data(filter_name);
    proj_weights = data.projection_weights(filter_name);

    if(not IS_PCDATA):
        wpca_data, e_val, e_vec = permutation_wPCA(proj_data, proj_weights, components=50, p_threshold=0.05, verbose=True);
        pcdata = PCData(wpca_data, e_val, e_vec, data);
    else:
        pcdata = data;

    result = pcdata.T; #Now rows are samples, columns are components

    result12 = result[:,[0,1]]
    projections['PCA: 1,2'] = result12;

    result23 = result[:,[1,2]]
    projections['PCA: 2,3'] = result23;

    result13 = result[:,[0,2]]
    projections['PCA: 1,3'] = result13;
    pbar.update();

    for method in method_list:
        result = method_list[method](proj_data, proj_weights);
        projections[method] = result;
        pbar.update();

    # Input projections
    # Load any projections that were supplied as an input file
    projections.update(apply_input_projections(input_projections, data.col_labels));
    

    #Normalize projections
    #Mean-center X
    #Mean-center Y
    #Scale so that the 90th percentile radius = 1
    for p in projections:
        coordinates = projections[p];
        coordinates[:,0] = coordinates[:,0] - np.mean(coordinates[:,0]);
        coordinates[:,1] = coordinates[:,1] - np.mean(coordinates[:,1]);

        r = np.sum(coordinates**2, axis=1)**(0.5);
        r90 = np.percentile(r, 90);

        if(r90 > 0):
            coordinates /= r90;
        
        coordinates = coordinates.T;        
        
        projections[p] = coordinates;
    
    pbar.complete();

    return projections, pcdata;


def perform_PCA(data, N=0, variance_proportion=1.0):
    """
    Performs PCA on the data
    
    Parameters
    ----------
    data : (Num_Features x Num_Samples) numpy.ndarray 
        Matrix containing data to project into 2 dimensions
    N : int
        Number of Principle Components to retain
    variance_proportion: float (0.0 - 1.0)
        Retain top X principal components such that a total of <variance_proportion>
        of the variance is retained.

    Returns
    -------
    pca_data : (Num_Components x Num_Samples) numpy.ndarray
        Data transformed using PCA.  Num_Components = Num_Samples

    """
    # PCA
    # Note, by default, PCA will subtract out the mean, but will not divide by
    # the variance.  Use the parameter 'whiten = True' for this.
    pca = PCA();
    pca_data = pca.fit_transform(data.T);
    
    if(N > pca_data.shape[1]): N = pca_data.shape[1];
    
    if(N != 0): #If N is specified, then return top N PC's
        pca_data = pca_data[:, 0:N];
    else: #Otherwise, if variance_proportion is specified, then return top PCs until variance proportion is reached
        total_variance = np.cumsum(pca.explained_variance_ratio_);
        last_i = np.nonzero(total_variance <= variance_proportion)[0][-1];
        pca_data = pca_data[:, 0:last_i + 1];

    pca_data = pca_data.T

    output = PCData(pca_data, pca.explained_variance_, data);
    return output;

def perform_weighted_PCA(data, weights, max_components=200):
    """
    Performs Weighted PCA on the data

    Parameters
    ----------
    data : (Num_Features x Num_Samples) numpy.ndarray (or subclass)
        Matrix containing data to project into 2 dimensions

    weights : (Num_Features x Num_Samples) numpy.ndarray (or subclass)
        Matrix containing weights to use for each coordinate in data

    max_components: int
        Maximum number of components to calculate

    Returns
    -------
    pca_data : (Num_Components x Num_Samples) numpy.ndarray
        Data transformed using PCA.  Num_Components = Num_Samples

    """
    proj_data = data;

    #Weighted means
    wmean = np.sum(proj_data * weights, axis=1) / np.sum(weights, axis=1);
    wmean = wmean.reshape((wmean.size, 1));

    data_centered = proj_data - wmean;
    weighted_data_centered = data_centered * weights;

    denom = np.dot(weights,weights.T)
    denom[denom == 0] = 1.0 # Avoid 0/0 : just make it 0/1

    wcov = np.dot(weighted_data_centered, weighted_data_centered.T) / denom

    model = PCA(n_components=min(proj_data.shape[0], proj_data.shape[1], max_components), svd_solver='randomized');
    model.fit(wcov);
    e_vec = model.components_;

    wpca_data = np.dot(e_vec, data_centered);

    # Try a different determination of e_val
    #e_val = np.var(wpca_data, axis=1);
    #total_var = np.sum(np.var(proj_data, axis=1));
    #e_val /= total_var;

    e_val = model.explained_variance_ratio_

    return wpca_data, e_val, e_vec.T;


def filter_PCA(data, scores=None, N=0, variance_proportion=1.0, min_components = 0):
    """
    Removes PC's that correlate with scores across samples

    :param data: Data that has been PCA transformed.  ndarray (Num Components x Num Samples)
    :param scores: Values (1 per sample)
    :param N: Int.  Number of PCs to retain (performed first)
    :param variance_proportion: Float.  Proportion of variance to retain (performed last)
    :param min_components: Int.  Minimum # of components to retain if using variance_proportion
        Default is no minimum.
    :return: The data with some PC's (rows) removed
    """

    #Spearman correlation for each Component with each Sample

    ##Fast version
    # data_indices = np.argsort(data).argsort(); #Two argsorts give rank
    # score_indices = np.argsort(scores).argsort();
    #
    # n = data_indices.shape[1];
    #
    # rho = 1 - 6*np.sum((data_indices - score_indices)**2, axis=1)/(n*(n**2-1));
    #
    # #Use t-like statistic for significance
    # t = rho * np.sqrt((n-2) / (1-rho**2));
    # p_plus = scipy.stats.t.cdf(t*-1, n-2)*2;
    # p_minus = scipy.stats.t.cdf(t, n-2)*2;
    # p = p_plus;
    # p[t < 0] = p_minus[t < 0];
    # print("Done");

    #Not as fast, but fast enough for a few thousand genes
    if(N > 0):
        if(N > data.shape[0]):
            N = data.shape[0];
        data = data.subset_components(np.arange(N));

    if(scores is not None):
        rho = np.zeros(data.shape[0]);
        p = np.zeros(data.shape[0]);
        for i in range(data.shape[0]):
           rho[i], p[i] = scipy.stats.spearmanr(data[i,:], scores)

        good_pcs = np.nonzero(p > 1e-5)[0];

        data = data.subset_components(good_pcs);

    if(variance_proportion < 1.0):
        explained_variance_ratio = data.variance / np.sum(data.variance);
        total_variance = np.cumsum(explained_variance_ratio);
        good_pcs = np.nonzero(total_variance <= variance_proportion)[0];
        if(good_pcs.size == 0):  #This means the first PC is already more than total_variance
            last_i = 0;
        else:
            last_i = good_pcs[-1];
        if(last_i+1 < min_components): last_i = min_components-1;

        data = data.subset_components(np.arange(last_i+1));

    return data

def define_clusters(projections):
    """
    Creates several different clusterings of the data in projections.

    :param projections: dict(string, (2 x Num_Samples) numpy.ndarray)
        dictionary mapping the projection type (e.g. "tSNE") to an array containing
        the two-dimensional coordinates for each sample in the projection.

    :return: dict of string (projection name) =>
        (dict of string (cluster technique) => np.ndarray of size N_Samples (cluster assignments))
    """

    pbar = ProgressBar(4 * len(projections));

    out_clusters = dict();

    for key in projections:

        proj_data = projections[key];
        proj_clusters = dict();

        # K-means for k = 2-5
        for k in range(2, 6):
            clust_name = "K-Means, k=" + str(k);
            kmeans = MiniBatchKMeans(n_clusters=k);
            clust_assignments = kmeans.fit_predict(proj_data.T);
            proj_clusters.update({clust_name: clust_assignments});
            pbar.update();

        out_clusters.update({key: proj_clusters});

    pbar.complete();

    return out_clusters;


def apply_input_projections(input_projections, sample_names):
    """
    Returns projection coordinates supplied as file input when running FastProject

    Parameters
    ----------
    input_projections : dict
        Keys are of type str, representing projection names
        Values are of type 2xN pandas.DataFrame with column names matching
        sample names in `expressionMatrix`
    sample_names : list of str
        Labels for samples used to subset/order the coordinates

    Returns
    -------
    dict
        proj_name(string) -> coordinates(numpy.ndarray size n_samples x 2)
    """

    output = dict()
    for proj_name in input_projections:
        coordinates = input_projections[proj_name];
        coordinates = coordinates.loc[sample_names, :];
        output.update({proj_name:
                       coordinates.values});

    return output;


def permutation_wPCA(data, weights, components=50, p_threshold=0.05, verbose=False, debug=False):
    """
    Computes weighted PCA on data.  Returns only significant components.

    After performing PCA on the data matrix, this method then uses a permutation
    procedure based on Buja A and Eyuboglu N (1992) to asses components for significance.


    :param data: numpy.ndarray of shape Num_Features x Num_Samples
    :param weights: numpy.ndarray of shape Num_Features x Num_Samples
    :param components: int, Max components to calculate
    :param p_threshold: float, P-value to cutoff components at
    :param verbose: bool,
    :return: reduced data, numpy.ndarray of shape Num_Components x Num_Samples
    """

    np.random.seed(RANDOM_SEED);

    # Make sure components isn't too big for matrix dimension
    components = min(components, data.shape[0], data.shape[1]);


    NUM_REPEATS = 20;

    wpca_data, e_val, e_vec = perform_weighted_PCA(data, weights, components);

    bg_vals = np.zeros((NUM_REPEATS, components));
    bg_data = np.zeros(data.shape);
    bg_weights = np.zeros(data.shape);

    for i in range(NUM_REPEATS):
        for j in range(data.shape[0]):
            random_i = np.random.permutation(data.shape[1]);
            bg_data[j,:] = data[j,random_i];
            bg_weights[j,:] = weights[j,random_i];

        wpca_bg, bg_val, bg_vec = perform_weighted_PCA(bg_data, bg_weights, components);

        bg_vals[i,:] = bg_val;

    # Compute distributions (normal) for each column in bg_vals
    mu = bg_vals.mean(axis=0);
    sigma = bg_vals.std(axis=0);

    sigma[sigma == 0] = 1.0e-19

    p_vals = norm.sf((e_val - mu) / sigma);
    threshold_component_i = np.nonzero(p_vals > p_threshold)[0]

    if(threshold_component_i.size == 0):  # True if ALL PCs deemed significant
        threshold_component = wpca_data.shape[0];
    else:
        threshold_component = threshold_component_i[0];

    if(verbose):
        FP_Output("Permutation test on wPCA: ", str(threshold_component), " components retained.");

    if(threshold_component < 5):
        threshold_component = 5;

        if(verbose):
            FP_Output("Less than 5 components identified as significant.  Preserving top 5.");

    if(debug):
        import seaborn as sns;
        import pandas as pd;
        import matplotlib.pyplot as plt;
        fig = plt.figure()
        sns.boxplot(pd.DataFrame(bg_vals[:,0:20]));
        plt.plot(e_val);

    wpca_data = wpca_data[0:threshold_component, :];
    e_val = e_val[0:threshold_component];
    e_vec = e_vec[:, 0:threshold_component];

    return wpca_data, e_val, e_vec;


# --------------------------------------------------------------------------- #
#                                                                             #
#                    Define Projection Methods Here                           #
#                                                                             #
# --------------------------------------------------------------------------- #

# Built-in Methods


# ICA
def apply_ICA(proj_data, proj_weights=None):
    ica = FastICA(n_components=2, random_state=RANDOM_SEED);
    result = ica.fit_transform(proj_data.copy().T);  # Copy needed because ICA whitens the input matrix
    return result;


# tSNE
def apply_tSNE10(proj_data, proj_weights=None):
    model = TSNE(n_components=2, perplexity=10.0, metric="euclidean",
                 learning_rate=200, early_exaggeration=4.0,
                 random_state=RANDOM_SEED);
    result = model.fit_transform(proj_data.T);
    return result;


# tSNE
def apply_tSNE30(proj_data, proj_weights=None):
    model = TSNE(n_components=2, perplexity=30.0, metric="euclidean",
                 learning_rate=200, early_exaggeration=4.0,
                 random_state=RANDOM_SEED);
    result = model.fit_transform(proj_data.T);
    return result;


# ISOMap
def apply_ISOMap(proj_data, proj_weights=None):
    model = Isomap(n_neighbors=4, n_components=2);
    result = model.fit_transform(proj_data.T);
    return result;


# PCA with RBF Kernel
def apply_rbf_PCA(proj_data, proj_weights=None):
    model = KernelPCA(n_components=2, kernel='rbf');
    result = model.fit_transform(proj_data.T);
    return result;


# MDS
def apply_MDS(proj_data, proj_weights=None):
    model = MDS(n_components=2, dissimilarity="euclidean",
            random_state=RANDOM_SEED)
    result = model.fit_transform(proj_data.T);
    return result;


# Spectral Embedding
def apply_spectral_embedding(proj_data, proj_weights=None):
    model = SpectralEmbedding(n_components=2, random_state=RANDOM_SEED)
    result = model.fit_transform(proj_data.T);
    return result;

# Add New Methods Here


# Gets the current projection methods to use
def get_projection_methods(lean, pca):
    """
    lean: bool
        Set by CLI argument.  Used to remove methods that don't scale well
        to many samples
    pca: bool
        Only return a subset of methods when PCA has been called first as
        it doesn't make sense to call PCA before some methods
    """

    proj_methods = {}

    if not pca:
        proj_methods['ICA'] = apply_ICA

        if not lean:
            proj_methods['Spectral Embedding'] = apply_spectral_embedding
            proj_methods['MDS'] = apply_MDS

        proj_methods['RBF Kernel PCA'] = apply_rbf_PCA
        proj_methods['ISOMap'] = apply_ISOMap
        proj_methods['tSNE30'] = apply_tSNE30
        proj_methods['tSNE10'] = apply_tSNE10

    if pca:

        proj_methods['ISOMap'] = apply_ISOMap
        proj_methods['tSNE30'] = apply_tSNE30
        proj_methods['tSNE10'] = apply_tSNE10

        if not lean:
            proj_methods['Spectral Embedding'] = apply_spectral_embedding
            proj_methods['MDS'] = apply_MDS

    return proj_methods
