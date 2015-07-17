# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15

@author: David DeTomaso
"""
from __future__ import print_function, division;

#Compute the distance matrix between samples
#Input:  Data = matrix of NumGenes x NumSamples
#Output:  Similarity matrices of size NumSamples x NumSamples

import numpy as np
import os;
from . import Utils

this_directory = os.path.dirname(os.path.abspath(__file__));

def filter_genes_fano(data, num_mad):
    """
    Uses fano filtering on the genes.  Splits into quantiles first.
    Only retains genes that have fano factor <num_mad> median absolute
        deviations in each quantile.

    :param data: numpy.ndarray (Num_Genes x Num_Samples)
    :param num_mad: float
    :return: numpy.ndarray, bool, (Num_Genes)
        True for rows that pass the filter.  Otherwise, False;
    """
    mu = np.mean(data, axis=1);
    sigma = np.std(data, axis=1);

    #slice up mu and sigma into bins, by mu
    aa = np.argsort(mu);
    mu_sort = mu[aa];
    sigma_sort = sigma[aa];

    N_QUANTS = 30;
    m = mu_sort.size//N_QUANTS;

    gene_passes = np.zeros(data.shape[0]) == 1;

    for i in xrange(N_QUANTS):
        if(i == N_QUANTS-1):
            rr = np.arange(i*m, mu_sort.size)
        else:
            rr = np.arange(i*m, (i+1)*m);

        mu_quant = mu_sort[rr];
        mu_quant[mu_quant == 0] = 1; #so we don't divide by zero later
        sigma_quant = sigma_sort[rr];
        fano_quant = sigma_quant**2 / mu_quant;
        mad_quant = np.median(np.abs(fano_quant - np.median(fano_quant)));
        gene_passes_quant = fano_quant > np.median(fano_quant) + num_mad*mad_quant;
        gene_passes_quant_i = np.nonzero(gene_passes_quant)[0];
        gene_passes_i = gene_passes_quant_i + i*m;
        gene_passes[gene_passes_i] = True;

    #gene_passes are in sorted mu order, need to revert back to original order
    original_ii = np.argsort(aa);
    gene_passes = gene_passes[original_ii];

    return gene_passes;



def filter_genes_threshold(data, threshold):
    """Filters out genes that are at least active in <threshold> 
    of samples.  Threshold ranges from 0 to 1"""
    
    row_length = data.shape[1];
    cutoff = row_length * threshold;
    
    keep_indices = np.where((data > 0).sum(axis=1) > cutoff)[0];
    
    return data.subset_genes(keep_indices);
    

def filter_genes_hdt(data, p_val):
    """Filters out genes that pass the Hartigans Dip Test for bimodality
    with at least p < p_val"""
    #perform Hartigans dip test on the rest of the rows

    (dips, ps, xlows, xups) = Utils.HDT_Sig_batch(data, 1000);

    return ps <= p_val;

def remove_from_file(data, filename):
    """trim out rows that match genes found in file <filename>
    First row is ignored
    if multiple entries per line, seperated by comma, take the first"""

    #load gene names from file
    ff = open(filename);
    xx = ff.readline();  #discard first line
    
    hk_genes = list();
    for line in ff:
      entries = line.split(',');
      hk_genes.append(entries[0].strip().lower());
    
    ff.close();
    
    #match hk genes to gene indices
    missing = 0;
    hk_indices = list();
    lower_row_labels = [gene.lower() for gene in data.row_labels];
    for hk_gene in hk_genes:
      try:
        ii = lower_row_labels.index(hk_gene);
        hk_indices.append(ii);
      except ValueError:
        missing+=1;
    
    #remove rows that match
    all_indices = np.arange(data.shape[0]);
    if(len(hk_indices) != 0):    
        keep_indices = np.setxor1d(all_indices,hk_indices); 
    else:
        keep_indices = all_indices;
    
    return data.subset_genes(keep_indices);

def filter_housekeeping(data, housekeeping_file=""):
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
        data = remove_from_file(data, hkf);
    
    
    return data;
        
        
def save_filter(data, filename):
    #save names of genes to a file
    ff = open(filename,'w');
    for gene in data.row_labels:
      ff.write(gene + '\n');
    
    ff.close();
    
def load_from_file(data, filename):
    loaded_genes = list();
    
    ff = open(filename,'r')
    for line in ff.readlines():
      loaded_genes.append(line.strip().lower());
    
    ff.close();
    
    keep_indices = list();
    for lgene in loaded_genes:
      for i, gene in enumerate(data.row_labels):
        if(gene.lower() == lgene):
          keep_indices.append(i);
          continue;
    
    data = data.subset_genes(keep_indices);
    
    return data;