# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 16:34:49 2015

@author: David
"""

import numpy as np;
from sklearn.neighbors import NearestNeighbors

def read_signature(filename='', match_terms=[]):
    """Reads in a signature file.  Returns a list of Signature objects
    
    Parameters
    ----------
    filename : string
        Name (and path if not in working directory) of signature file to read
        If no file is entered, opens a file dialog
    match_terms : list(String)
        List of terms to be matched against signature names.
        Signature is retained if any of the terms match.
    
    Returns
    -------
    signatures : list(DimReduce.Signature)
        The filtered signature list
    """
    
    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();
    
    ff = open(filename, 'r');
    found_signatures = dict();    ## dictionary of signatures
    
    if(type(match_terms) is str):
        match_terms = [match_terms];
    
    match_terms = map(lambda term: term.lower(), match_terms);    
    
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
                
            name = row_data[0];
            if(not found_signatures.has_key(name)):  ## Add the signature if we haven't seen it yet
                
                matched = False if len(match_terms) > 0 else True                ## Only add signatures if they match one of the supplied terms
                lname = name.lower();                
                for term in match_terms:
                    if(lname.find(term) >= 0):
                        matched = True;
                        break;
                
                if(matched):
                    sig = Signature(dict(), signed, filename, name);
                    found_signatures.update({name: sig});
                else:
                    continue;
            
            sig_dict = found_signatures[name].sig_dict;
            
            
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
    
    return [found_signatures[key] for key in found_signatures.keys()]  #dict to list

def filter_sig_list(signatures, match_terms):
    """
    Filters the list of signature objects to retain only those that contain
    one of the specified terms
    
    Parameters
    ----------
    signatures : list(Signature)
        List of signatures to be filtered
    match_terms : list(String)
        List of terms to be matched against signature names.
        Signature is retained if any of the terms match.
    
    Returns
    -------
    filtered_signatures : list(DimReduce.Signature)
        The filtered signature list
        
    
    """
    
    if(type(match_terms) is str):       #So you dont need to wrap a single string in []
        match_terms = [match_terms];    
    
    match_terms = map(lambda term: term.lower(), match_terms);    
    
    filtered_signatures = list();
    for sig in signatures:
        name = sig.name.lower();
        
        for term in match_terms:            
            if(name.find(term) >= 0):
                filtered_signatures.append(sig);
                break;
    
    
    return filtered_signatures;


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
    
    def __init__(self, sig_dict, signed, filename, name):
        self.sig_dict = sig_dict;
        self.signed = signed;
        self.source = filename;
        self.name = name;
    
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
    