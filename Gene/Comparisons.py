# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15

@author: David DeTomaso
"""


#Compute the distance matrix between samples
#Input:  Data = matrix of NumGenes x NumSamples
#Output:  Similarity matrixs of size NumSamples x NumSamples

import numpy as np

def euclidean_distance(data):
    return np.cov(data.T);
    
    
#Entries of data represent the probability of some occurrence
#Diff(i,j) = pi*(1-pj) + (1-pj)*pi
def probability_difference(data):
    p = data;
    pi = 1-data;
    
    return np.dot(p.T, pi) + np.dot(pi.T,p);
    

#Same as computing:
#Diff(i,j) = pi*pj + (1-pi)*(1-pj)
def probability_similarity(data):
    return 1 - probability_difference(data);
