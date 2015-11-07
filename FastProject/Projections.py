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

from .Utils import ProgressBar;
from .DataTypes import PCData;

import numpy as np;


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

    pbar = ProgressBar(8);
    
    projections = dict();
    


    proj_data = data.projection_data(filter_name);
    proj_weights = data.projection_weights(filter_name);
    dist_matrix = data.distance_matrix(filter_name);

    # PCA - Uses subsampling inside perform_weighted_pca
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
    
    # tSNE with higher perplexity
    model = TSNE(n_components=2, perplexity=30.0, metric="precomputed", learning_rate = 100, early_exaggeration=4.0);
    result = model.fit_transform(dist_matrix);
    
    projections['tSNE30'] = result;
    pbar.update();

    # tSNE with lower perplexity
    model = TSNE(n_components=2, perplexity=1.0, metric="precomputed", learning_rate = 100, early_exaggeration=4.0);
    result = model.fit_transform(dist_matrix);

    projections['tSNE1'] = result;
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
    
    model = MDS(n_components=2, dissimilarity="precomputed")
    result = model.fit_transform(dist_matrix);
    
    projections['MDS'] = result;
    pbar.update();
        
    # Spectral Embedding
    # Issues using precomputed affinity matrix.  Need to understand how to construct it better
    model = SpectralEmbedding(n_components=2)
    result = model.fit_transform(proj_data.T);
    
    projections['Spectral Embedding'] = result;
    pbar.update();
    
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





