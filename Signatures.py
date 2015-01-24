# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 16:34:49 2015

@author: David
"""

import numpy as np;
from sklearn.neighbors import NearestNeighbors

def read_signature(filename=''):
    """Reads in a signature file.  Returns an instance of the Signature class"""
    
    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();
    
    ff = open(filename, 'r');
    sig_dict = dict();    
    
    try:    
        
        firstline = ff.readline();
        
        num_cols = len(firstline.split('\t'));
        
        if(num_cols == 2):
            signed = False;
        elif(num_cols == 3):
            signed = True;
        else:
            raise ValueError("Signature file should contain 2 (unsigned) or 3 (signed) columns")
        
        ff.seek(0);
        
        
        for i, line in enumerate(ff):
            row_data = line.strip().split('\t');
            if(len(row_data) != 2 and len(row_data) != 3):
                print i
            if signed:
                sig_sign = row_data[1].lower();
                if(sig_sign == 'plus'):
                    sig_val = 1;
                elif(sig_sign == 'minus'):
                    sig_val = -1;
                elif(sig_sign == 'both'):
                    sig_val = 0;
                elif(sig_sign == 'mius'):  #LOL, close enough
                    sig_val = -1;
                ## Other spelling variants go here
                else:
                    print "Error on line ", str(i), " Couldn't read signature value."
                    print "   :", line;
                    continue;
                    
                sig_dict[row_data[2]] = sig_val;
                
            else:  #unsigned case
                sig_val = 0;
                sig_dict[row_data[1]] = sig_val;
     
    except:
        raise;
    finally:
        ff.close();
        
    
    sig = Signature(sig_dict, signed, filename);
    
    return sig;

def conformity(data_loc, sig, n_neighbors):
    """
    Score each sample based on how similar its signature score is to its neighborhood
    
    Parameters
    ---------
    data_loc : array-like, shape (Num_Dimensions, Num_Samples)
        locations of points used to calculate neareset neighbors
    sig : array-like, 1D, shape (Num_Samples)
        Signature value for each sample.  Get using Signature.eval_data
    n_neighbors : int
        Number of neighbors to use in defining the neighborhood
    
    Returns
    -------
    dissimilarity : array-like, 1D, shape (Num_Samples)
        Score for each sample as to how dissimilar it is from neighboring points
        
    """
    
    #Incremenet n_neighbors since it counts a point as it's own neighbor
    nbrs = NearestNeighbors(n_neighbors=n_neighbors+1, algorithm='ball_tree');
    nbrs.fit(data_loc.T);
    
    distances, indices = nbrs.kneighbors(data_loc.T);
    
    neighborhood = sig[indices];


    ##Weighted mean of neighborhood point signatures defines a prediction
    ##Weights are 1/distance
    weights = distances ** -1;

    neighborhood_prediction = np.sum(neighborhood[:,1:] * weights[:,1:]) \
                / np.sum(weights[:,1:]);
    
    
    ##Neighborhood dissimilarity score = |actual - predicted|
    dissimilarity = np.abs(sig - neighborhood_prediction);
    
    return dissimilarity;

    
class Signature:
    
    def __init__(self, sig_dict, signed, filename):
        self.sig_dict = sig_dict;
        self.signed = signed;
        self.source = filename;
    
    def eval_data(self, data, genes):
        """For each sample, calculate a score against the signature
        
        Parameters
        ----------
        data :   array-like, shape (Num_Genes, Num Samples)
            expression data to evaluate the signature against
        genes :  list of strings
            list containing gene names (symbol) for each row in data
        
        Returns
        -------
        out : 1D ndarray, length = Num_Samples 
            result of evaluating this signature on each sample in data
        
        """
        
        sig_vector = self._sig_indices(genes);
        
        if(data.ndim == 1):
            data.shape = (data.shape[0], 1);
            
        pdata = data * sig_vector;
        
        return pdata.sum(axis = 0);
        
    def _sig_indices(self, genes):
        """Helper method
        
        Returns an array with length = len(genes)
        Entries in the array are 0 if the gene is not in the signature
        Otherwise, value is determined by the signature type"""
        
        out = np.zeros((len(genes), 1));
        
        for i, gene in enumerate(genes):
            if(self.sig_dict.has_key(gene)):
                val = self.sig_dict[gene];
                if(val == 1 or val == 0):
                    out[i] = 1;
                if(val == -1):
                    out[i] = -1;
            
        return out;
    