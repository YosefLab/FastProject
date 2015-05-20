# -*- coding: utf-8 -*-
"""
Created on Wed Feb 04 13:21:47 2015

@author: David
"""
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA
from sklearn.manifold import TSNE;
from sklearn.manifold import Isomap
from sklearn.manifold import LocallyLinearEmbedding
from sklearn.manifold import MDS
from sklearn.manifold import SpectralEmbedding
import scipy.stats;

from .Utils import ProgressBar;
from . import DataTypes

import numpy as np;

def write_projection_file(filename, sample_labels, projections):
    """
    Outputs the coordinates for each projection to a file.
    
    Parameters
    ----------
    filename : String
        File to write
    sample_labels : list(String) len=Num_Samples
        Labels for each sample
    projections : dict(string, (2 x Num_Samples) numpy.ndarray)
        dictionary mapping the projection type (e.g. "tSNE") to an array containing
        the two-dimensional coordinates for each sample in the projection.
          
    """
    
    ff = open(filename,'w');
    
    for proj in projections.keys():
        coordinates = projections[proj];
        for i in range(coordinates.shape[1]):
            ff.write(proj + '\t');
            ff.write(sample_labels[i] + '\t');
            ff.write('{:.5f}'.format(coordinates[0,i]) + '\t');
            ff.write('{:.5f}'.format(coordinates[1,i]) + '\t');
            ff.write('\n');
    
    ff.close();
    

def generate_projections(data):
    """
    Projects data into 2 dimensions using a variety of linear and non-linear methods
    
    Parameters
    ----------
    data : (Num_Features x Num_Samples) numpy.ndarray 
        Matrix containing data to project into 2 dimensions

    Returns
    -------
    projections : dict(string, (2 x Num_Samples) numpy.ndarray)
        dictionary mapping the projection type (e.g. "tSNE") to an array containing
        the two-dimensional coordinates for each sample in the projection.
          
    """
    
    pbar = ProgressBar(7);
    
    projections = dict();
    
    dist_matrix = data.distance_matrix();    
    
    # PCA
    
    pca = PCA();
    result = pca.fit_transform(data.T);
    result = result[:,[0,1]]
  
    projections['PCA'] = result;
    pbar.update();
      
    # ICA
    
    ica = FastICA(n_components = 2);
    result = ica.fit_transform(data.T);
    
    projections['ICA'] = result;
    pbar.update();
    
    # tSNE, but with built-in from sklearn
    model = TSNE(n_components=2, perplexity=1, metric="precomputed", learning_rate = 100, early_exaggeration=1);
    result = model.fit_transform(dist_matrix);
    
    projections['tSNE'] = result;
    pbar.update();
    
    # ISOMap
    
    model = Isomap(n_neighbors = 4, n_components = 2);
    result = model.fit_transform(data.T);
    
    projections['ISOMap'] = result;
    pbar.update();
    
    # LLE
    
    model = LocallyLinearEmbedding(n_neighbors = 2, n_components=2)
    result = model.fit_transform(data.T);
    
    projections['LLE'] = result;
    pbar.update();
    
    # MDS
    
    model = MDS(n_components=2, dissimilarity="precomputed")
    result = model.fit_transform(dist_matrix);
    
    projections['MDS'] = result;
    pbar.update();
        
    # Spectral Embedding
    if(type(data) is DataTypes.ProbabilityData):
        model = SpectralEmbedding(n_components=2, affinity="precomputed")
        result = model.fit_transform(1-dist_matrix);
    else:
        model = SpectralEmbedding(n_components=2)
        result = model.fit_transform(data.T);
    
    projections['Spectral Embedding'] = result;
    pbar.update();
    
    # Add new projections here!


    #Normalize projections 
    #Mean-center X
    #Mean-center Y
    #Scale so that E(R^2) = 1
    for p in projections:
        coordinates = projections[p];
        coordinates[:,0] = coordinates[:,0] - np.mean(coordinates[:,0]);
        coordinates[:,1] = coordinates[:,1] - np.mean(coordinates[:,1]);
        
        ave_r2 = (coordinates*coordinates).mean(axis=0).sum();
        
        coordinates = coordinates / np.sqrt(ave_r2);
        
        coordinates = coordinates.T;        
        
        projections[p] = coordinates;
    
    pbar.complete();

    return projections; 

def perform_PCA(data, N=0, variance_proportion=1.0):
    """
    Performs PCA on the data
    
    Parameters
    ----------
    data : (Num_Features x Num_Samples) numpy.ndarray 
        Matrix containing data to project into 2 dimensions
    row_labels : list(String) len=Num_Features
        Former row labels
    N : int
        Number of Principle Components to retain

    Returns
    -------
    pca_data : (Num_Components x Num_Samples) numpy.ndarray
        Data transformed using PCA.  Num_Components = Num_Samples
    row_labels : list(String) len=N
        New row labels since former labels no longer make sense
          
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
    
    return pca_data.T;

def filter_PCA(data, scores):
    """
    Removes PC's that correlate with scores across samples

    :param data: Data that has been PCA transformed.  ndarray (Num Components x Num Samples)
    :param scores: Values (1 per sample)
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
    rho = np.zeros(data.shape[0]);
    p = np.zeros(data.shape[0]);
    for i in xrange(data.shape[0]):
       rho[i], p[i] = scipy.stats.spearmanr(data[i,:], scores)

    good_pcs = np.nonzero(p > 1e-5)[0];

    data = data.subset_components(good_pcs);
    return data
