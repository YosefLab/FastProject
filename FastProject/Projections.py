# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 13:21:47 2015

@author: David
"""
from sklearn.decomposition import PCA
from sklearn.decomposition import RandomizedPCA
from sklearn.decomposition import FastICA
from sklearn.decomposition import KernelPCA
from sklearn.manifold import TSNE;
from sklearn.manifold import Isomap
from sklearn.manifold import MDS
from sklearn.manifold import SpectralEmbedding
from sklearn.cluster import MiniBatchKMeans
import scipy.stats;
import pandas as pd;
from pandas.parser import CParserError;
import os;

from .Utils import ProgressBar;
from .DataTypes import PCData;
from .Global import args;

import numpy as np;

# Input projection coordinates are read once (at start) and cached here
# Variable stores a list of pandas.DataFrame objects
_input_projections = dict();
def load_input_projections(sample_names):
    """
    Loads input projection coordinates into module variable _input_projections

    :param sample_names: List of sample names.  Used to validate input
    :return: None
    """
    global _input_projections;

    sample_names_norm = [x.upper() for x in sample_names];

    if(args.projections):
        for projection_file in args.projections:
            if(not os.path.isfile(projection_file)):
                raise ValueError("Option Error: projection file " + projection_file + " not found.");

            try:
                loaded_coords = pd.read_table(projection_file, header=None);
            except CParserError:
                # if Header row has 2 entries and rest has 3, get this error
                # Try to load data anyway and rework it back into the correct initial format
                loaded_coords = pd.read_table(projection_file);
                loaded_coords = pd.DataFrame({0: loaded_coords.index,
                                              1: loaded_coords.iloc[:,0].values,
                                              2: loaded_coords.iloc[:,1].values},
                                             index=np.arange(loaded_coords.shape[0]));


            # Verify the file
            # Check for 3 columns
            if(loaded_coords.shape[1] != 3):
                raise ValueError("Format Error: projection file " + projection_file +
                                 " should have 3 columns (sample name, x coordinate, y coordinate).");

            # Fix in case header was included with first cell empty
            if(pd.isnull(loaded_coords.iloc[0,0])):
                loaded_coords.iloc[0,0] = "_____"; #Dummy value, won't match gene

            # Normalize names to upper-case temporarily
            loaded_coords[0] = [x.upper() for x in loaded_coords[0]];

            # Did they accidentally include a heading?
            if(loaded_coords.iloc[0,0] not in sample_names_norm):
                loaded_coords = loaded_coords.drop(loaded_coords.index[0], axis='index'); # drop first row

            # Check for duplicate names
            if(loaded_coords[0].duplicated().any()):
                duplicated_names = loaded_coords.loc[loaded_coords[0].duplicated(),0].values;
                raise ValueError("Format Error: projection file " + projection_file +
                                 " contains duplicate sample names:\n" + "\n".join(duplicated_names));

            # Set as the index
            loaded_coords = loaded_coords.set_index(0);

            # Check that all data matrix sample names are in the projection file
            not_matched = np.array([x not in loaded_coords.index for x in sample_names_norm]);
            if(not_matched.any()):
                missing_samples = [sample_names[i] for i in np.nonzero(not_matched)[0]];
                raise ValueError("Format Error: projection file " + projection_file +
                                 " missing entries for:\n" + "\n".join(missing_samples));

            # If we got here, it's safe to re-order the file
            loaded_coords = loaded_coords.loc[sample_names_norm,:];
            loaded_coords.index = sample_names;  # Restore original case
            loaded_coords["X"] = loaded_coords[1];
            loaded_coords["Y"] = loaded_coords[2];
            loaded_coords = loaded_coords.drop(1, axis="columns").drop(2, axis="columns");

            # Check for null values and throw error
            if(loaded_coords.isnull().any().any()):
                bad_rows = loaded_coords.index[loaded_coords.isnull().any(axis=1)];
                raise ValueError("Format Error: projection file " + projection_file +
                                 " contains NaNs for samples:\n" + "\n".join(bad_rows));

            proj_name = os.path.splitext(projection_file)[0];
            _input_projections.update({proj_name: loaded_coords});


def generate_projections(data, filter_name = None):
    """
    Projects data into 2 dimensions using a variety of linear and non-linear methods
    
    Parameters
    ----------
    data : (Num_Features x Num_Samples) numpy.ndarray 
        Matrix containing data to project into 2 dimensions
    filter_name : String
        Filter to apply to the signatures.  Should match a filter in data.filters
        If not specified, no filtering is applied.

    Returns
    -------
    projections : dict(string, (2 x Num_Samples) numpy.ndarray)
        dictionary mapping the projection type (e.g. "tSNE") to an array containing
        the two-dimensional coordinates for each sample in the projection.

    PC_Data : (Num_Components x Num_Samples) numpy.ndarray
        Weighted PCA of the original data object
          
    """

    pbar = ProgressBar(7);
    
    projections = dict();
    


    proj_data = data.projection_data(filter_name);
    proj_weights = data.projection_weights(filter_name);

    if(type(data) is not PCData):
        wpca_data, e_val, e_vec = perform_weighted_PCA(proj_data, proj_weights, max_components = 50);
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

    # ICA
    
    ica = FastICA(n_components = 2);
    result = ica.fit_transform(proj_data.T.copy()); #Copy needed because ICA whitens the input matrix
    
    projections['ICA'] = result;
    pbar.update();
    
    # tSNE
    model = TSNE(n_components=2, perplexity=10.0, metric="euclidean", learning_rate = 100, early_exaggeration=4.0);
    result = model.fit_transform(proj_data.T);
    
    projections['tSNE10'] = result;
    pbar.update();

    # tSNE
    model = TSNE(n_components=2, perplexity=30.0, metric="euclidean", learning_rate = 100, early_exaggeration=4.0);
    result = model.fit_transform(proj_data.T);
    
    projections['tSNE30'] = result;
    pbar.update();

    # ISOMap
    
    model = Isomap(n_neighbors = 4, n_components = 2);
    result = model.fit_transform(proj_data.T);
    
    projections['ISOMap'] = result;
    pbar.update();
    
    # PCA with RBF Kernel

    model =  KernelPCA(n_components=2, kernel='rbf');
    result = model.fit_transform(proj_data.T);
    
    projections['RBF Kernel PCA'] = result;
    pbar.update();
    
    # MDS

    model = MDS(n_components=2, dissimilarity="euclidean")
    result = model.fit_transform(proj_data.T);

    projections['MDS'] = result;
    pbar.update();

    # Spectral Embedding
    # Issues using precomputed affinity matrix.  Need to understand how to construct it better
    model = SpectralEmbedding(n_components=2)
    result = model.fit_transform(proj_data.T);
    
    projections['Spectral Embedding'] = result;
    pbar.update();

    # Input projections
    # Load any projections that were supplied as an input file
    projections.update(apply_input_projections(data.col_labels));
    
    # Add new projections here!



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
        pca_data = pca_data[:,range(N)];
    else: #Otherwise, if variance_proportion is specified, then return top PCs until variance proportion is reached
        total_variance = np.cumsum(pca.explained_variance_ratio_);
        last_i = np.nonzero(total_variance <= variance_proportion)[0][-1];
        pca_data = pca_data[:,range(last_i+1)];

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

    wcov = np.dot(weighted_data_centered, weighted_data_centered.T) / np.dot(weights,weights.T);
    model = RandomizedPCA(n_components=min(proj_data.shape[0], proj_data.shape[1], max_components));
    model.fit(wcov);
    e_vec = model.components_;

    wpca_data = np.dot(e_vec, data_centered);
    e_val = np.var(wpca_data, axis=1);
    total_var = np.sum(np.var(proj_data, axis=1));
    e_val /= total_var;

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
        for i in xrange(data.shape[0]):
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

    for key in projections.keys():

        proj_data = projections[key];
        proj_clusters = dict();

        # K-means for k = 2-5
        for k in xrange(2, 6):
            clust_name = "K-Means, k=" + str(k);
            kmeans = MiniBatchKMeans(n_clusters=k);
            clust_assignments = kmeans.fit_predict(proj_data.T);
            proj_clusters.update({clust_name: clust_assignments});
            pbar.update();

        out_clusters.update({key: proj_clusters});

    pbar.complete();

    return out_clusters;

def apply_input_projections(sample_names):
    """
    Returns projection coordinates supplied as file input when running FastProject

    :param sample_names: Labels for samples used to subset/order the coordinates
    :return: dict or proj_name(string) -> coordinates(numpy.ndarray size n_samples x 2)
    """

    output = dict()
    for proj_name in _input_projections:
        coordinates = _input_projections[proj_name];
        coordinates = coordinates.loc[sample_names,:];
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

    NUM_REPEATS = 20;

    wpca_data, e_val, e_vec = perform_weighted_PCA(data, weights, components);

    bg_vals = np.zeros((NUM_REPEATS, components));
    bg_data = np.zeros(data.shape);
    bg_weights = np.zeros(data.shape);

    for i in xrange(NUM_REPEATS):
        for j in xrange(data.shape[0]):
            random_i = np.random.permutation(data.shape[1]);
            bg_data[j,:] = data[j,random_i];
            bg_weights[j,:] = weights[j,random_i];

        print("starting wpca");
        wpca_bg, bg_val, bg_vec = perform_weighted_PCA(bg_data, bg_weights, components);
        print("finished wpca!");

        bg_vals[i,:] = bg_val;

    # Compute distributions (normal) for each column in bg_vals
    mu = bg_vals.mean(axis=0);
    sigma = bg_vals.std(axis=0);

    p_vals = norm.sf((e_val - mu) / sigma);
    threshold_component = np.nonzero(p_vals > p_threshold)[0][0];

    wpca_data = wpca_data[0:threshold_component,:];

    if(verbose):
        print("Permutation test on wPCA: ", str(wpca_data.shape[0]), " components retained.");

    if(debug):
        import seaborn as sns;
        import pandas as pd;
        import matplotlib.pyplot as plt;
        sns.violinplot(pd.DataFrame(bg_vals[:,0:20]));
        plt.plot(e_val);

    return wpca_data;


