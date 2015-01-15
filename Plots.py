# -*- coding: utf-8 -*-
"""
Created on Sat Jan 03 21:29:39 2015

@author: daved_000
"""
import numpy as np;
import re;

def gates_from_sample_labels(sample_names, pattern):
    
    #Iterate through cells, building a map of unique matches
    matches = dict();
    num = 1;
    
    gates = np.zeros(len(sample_names));
    
    for i,name in enumerate(sample_names):
        out = re.search(pattern,name);
        if out is not None:
            if out.group() not in matches.keys():
                matches[out.group()] = num;
                num = num + 1;
            gates[i] = matches[out.group()];
        else:
            print "Error: no match in element ", i, " , ", name;
        
    return matches, gates
    
    
    #Sample way to find a pattern like xxx_xxx
    #out = re.search('^[^_]+_[^_]+',test)
    #if out is not None
    #match = out.group()    
    
    #Print unique matches (for QC purposes)
    
    #Assign each cell a match
    
    
    #return a vector of numbers that can be used for shading
    

    
    

    