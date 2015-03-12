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

from .Utils import ProgressBar;

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
    
    model = TSNE(n_components=2, perplexity=1, metric="euclidean", learning_rate = 100, early_exaggeration=1);
    result = model.fit_transform(data.T);
    
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
    
    model = MDS(n_components=2)
    result = model.fit_transform(data.T);
    
    projections['MDS'] = result;
    pbar.update();
        
    # Spectral Embedding
    
    model = SpectralEmbedding(n_components=2)
    result = model.fit_transform(data.T);
    
    projections['Spectral Embedding'] = result;
    pbar.update();
    
    # Add new projections here!


    #Normalize projection 
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

def perform_PCA(data, N=0):
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
    
    if(N != 0):
        pca_data = pca_data[:,range(N)];
    
    row_labels = ["PC"+str(i+1) for i in range(N)];    
    
    return pca_data.T, row_labels;
    