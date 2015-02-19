# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15

@author: David DeTomaso
"""
from __future__ import print_function;

#Compute the distance matrix between samples
#Input:  Data = matrix of NumGenes x NumSamples
#Output:  Similarity matrixs of size NumSamples x NumSamples

import numpy as np
import os;
from Utils import HDT_Sig
from Utils import ProgressBar

this_directory = os.path.dirname(os.path.abspath(__file__));


def filter_genes_threshold(data, genes, threshold):
    """Filters out genes that are at least active in <threshold> 
    of samples.  Threshold ranges from 0 to 1"""
    
    row_length = data.shape[1];
    cutoff = row_length * threshold;
    
    keep_indices = np.where((data > 0).sum(axis=1) > cutoff)[0];
    
    return filter_genes_indices(data,genes,keep_indices);
    

def filter_genes_hdt(data, genes, p_val):
    """Filters out genes that pass the Hartigans Dip Test for bimodality
    with at least p < p_val"""
    #perform Hartigans dip test on the rest of the rows
    
    if(p_val > .5):
        print("Error, p_val must be less than 0.5")
        return;
    
    print();
    print('Filtering using HDT with p<' +str(p_val));        
    
    #first with cutoff p=.5 and 1 iteration
    print();
    print('First Pass');
    p_cut = 0.5;
    hdt_p = np.zeros(data.shape[0]);
    pp = ProgressBar(data.shape[0]);
    for i in np.arange(data.shape[0]):
      (dip, p, xlow, xup) = HDT_Sig(data[i,:],1);
      hdt_p[i] = p;
      pp.update();
    
    pp.complete();

    
    keep_indices = np.nonzero(hdt_p <= p_cut)[0];
    (data, genes) = filter_genes_indices(data,genes,keep_indices);
      
    
    
    #second with cutoff p=p_val and 10 iterations
    print();
    print('Second Pass');
    pp = ProgressBar(data.shape[0]);
    p_cut = p_val;
    hdt_p = np.zeros(data.shape[0]);
    for i in np.arange(data.shape[0]):
      (dip, p, xlow, xup) = HDT_Sig(data[i,:],10);
      hdt_p[i] = p;
      pp.update();

    pp.complete();    
    
    keep_indices = np.nonzero(hdt_p <= p_cut)[0];
    (data, genes) = filter_genes_indices(data,genes,keep_indices);
    
    #third with cutoff p=p_val and 40 iterations
#    if(p_val < .2):
#        p_cut = p_val;
#        print "Third pass";
#        hdt_p = np.zeros(data.shape[0]);
#        for i in np.arange(data.shape[0]):
#          (dip, p, xlow, xup) = HDT_Sig(data[i,:],40);
#          hdt_p[i] = p;
#          if(np.mod(i,data.shape[0]/20)==0):
#              print round(float(i)/data.shape[0]*100.0), '%'
#        
#        keep_indices = np.nonzero(hdt_p <= p_cut)[0];
#        (data, genes) = filter_genes_indices(data,genes,keep_indices);
    
    return (data, genes);
    
def remove_from_file(data,genes,filename):
    """trim out rows that match genes found in file <filename>
    First row is ignored
    if multiple entries per line, seperated by comma, take the first"""

    #load gene names from file
    ff = open(filename);
    xx = ff.readline();  #discard first line
    
    hk_genes = list();
    for line in ff:
      entries = line.split(',');
      hk_genes.append(entries[0].strip());
    
    ff.close();
    
    #match hk genes to gene indices
    missing = 0;
    hk_indices = list();
    for hk_gene in hk_genes:
      try:
        ii = genes.index(hk_gene);
        hk_indices.append(ii);
      except ValueError:
        missing+=1;
    
    #remove rows that match
    all_indices = np.arange(data.shape[0]);
    if(len(hk_indices) != 0):    
        keep_indices = np.setxor1d(all_indices,hk_indices); 
    else:
        keep_indices = all_indices;
    
    return filter_genes_indices(data,genes,keep_indices);

def filter_housekeeping(data, genes, housekeeping_file=""):
    """
    Filters out Housekeeping genes.  Uses specified file if provided.
    Otherwise, uses every file in the Housekeeping folder with the module.
    """
    housekeeping_dir = os.path.join(this_directory,'Housekeeping Genes');
    
    housekeeping_files = list();
    
    if(housekeeping_file != ""): #If file specified, use that file
        housekeeping_files.append(housekeeping_file);
    else:  #Otherwise, use all the files!
        files = os.listdir(housekeeping_dir);
        for ff in files:
            housekeeping_files.append(os.path.join(housekeeping_dir, ff));
    
    for hkf in housekeeping_files:
        (data, genes) = remove_from_file(data, genes, hkf);
    
    
    return (data, genes);
        


def filter_genes_indices(data, genes, indices):
    """Utility method.  Given a data matrix (genes x samples) and a list of
    gene names, filter the two only retaining entries at locations <indices>"""
    
    data = data[indices,:];
    genes = [genes[i] for i in indices];
    
    return (data, genes);
    
def remove_genes(data, genes, gene_indices):
    """Utility method.  Given a data matrix (genes x samples) and a list of
    gene names, remove entries at locations <indices>"""
    all_indices = np.arange(data.shape[0]);    
    keep_indices = np.setxor1d(all_indices,gene_indices); 
    return filter_genes_indices(data,genes,keep_indices);
    
    
def save_filter(genes, filename):
    #save names of genes to a file
    ff = open(filename,'w');
    for gene in genes:
      ff.write(gene + '\n');
    
    ff.close();
    
def load_from_file(data, genes, filename):
    loaded_genes = list();
    
    ff = open(filename,'r')
    for line in ff.readlines():
      loaded_genes.append(line.replace('\n',''));
    
    ff.close();
    
    keep_indices = list();
    for lgene in loaded_genes:
      for i, gene in enumerate(genes):
        if(gene == lgene):
          keep_indices.append(i);
          continue;
    
    (data, genes) = filter_genes_indices(data,genes,keep_indices);
    
    return (data, genes)