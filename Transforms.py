# -*- coding: utf-8 -*-
"""
Created on Wed Jan 07 15:12:02 2015

@author: David
"""



from DimReduce.Stats import em;
from DimReduce import Filters;
import numpy as np;
import os;

this_directory = os.path.dirname(os.path.abspath(__file__));

#from numba import jit;



def probability_of_expression(data, nozero=True):
    cutoffs = np.ones(data.shape[0]);
    if(nozero):
        (gamma, mu_l, mu_h, st_l, st_h, Pi, L) = em.em_exp_norm_mixture2(data,cutoffs);
    else:
        (gamma, mu_l, mu_h, st_l, st_h, Pi, L) = em.em_exp_norm_mixture(data,cutoffs);
        
    return (gamma, mu_h);

def make_monotonic(gamma, data):
    """For each row in gamma, finds the first corresponding value in data
    in which gamma hits it's maximum.  Then, for all values of data > than that
    value, set gamma = gamma_max
    
    Dimensions of Gamma should be equal to Dimensions of Data"""
    
    for i in np.arange(gamma.shape[0]):
        max_g = np.max(gamma[i,:]);
        locs = np.flatnonzero(gamma[i,:] == max_g);
        max_d = data[i,locs[0]];
        
        locs_to_change = data[i,:] > max_d;
        gamma[i,locs_to_change] = max_g;

    return gamma;        
    
def create_false_neg_map(data, genes, housekeeping_file=""):
    """Uses gene names in <filename> to create a mapping of false negatives. 
    This assumes all genes in <filename> are actually active, despite measured
    expression level.  Should use housekeeping genes.  If filename is blank,
    use all stored housekeeping gene names
    
    Returns out_func: a function that maps mu_h -> prob_false_negative    
    """
    
    housekeeping_dir = os.path.join(this_directory,'Housekeeping Genes');
    
    housekeeping_files = list();
    
    if(housekeeping_file != ""): #If file specified, use that file
        housekeeping_files.append(os.path.join(housekeeping_dir, housekeeping_file));
    else:  #Otherwise, use all the files!
        files = os.listdir(housekeeping_dir);
        for ff in files:
            housekeeping_files.append(os.path.join(housekeeping_dir, ff));
    
    data_hk = np.zeros((0, data.shape[1]));
    genes_hk = list();
    for hkf in housekeeping_files:
        (data_t, genes_t) = Filters.load_from_file(data, genes, hkf);        
        data_hk = np.vstack((data_hk, data_t));
        genes_hk.extend(genes_t);


    (data_hk, genes_hk) = Filters.filter_genes_threshold(data_hk, genes_hk, 0.2);
        
    #calculate distributions for hk gene
    cutoffs = np.ones(data_hk.shape[0]);
    (gamma, mu_l, mu_h, st_l, st_h, Pi, L) = em.em_exp_norm_mixture2(data_hk,cutoffs);
        
    
    #neg_vals = np.sum(data_hk==0, axis=1) / float(data_hk.shape[1]);    
    neg_vals = np.mean(1-gamma, axis=1); 

    #Fit a function mapping mu to neg_vals

    x = mu_h.flatten();
    y = neg_vals;
    
    def func(xvals, x0, a):
        return 1/(1 + np.exp((xvals-x0)*a));
    
    from scipy.optimize import curve_fit;
    param, cov = curve_fit(func, x, y);
    
    #Uncomment to plot result
#    from pylab import figure, scatter, plot, ion;
#    ion();
#    figure();
#    domain = np.linspace(0,10,1000);
#    scatter(x,y);
#    plot(domain, func(domain, param[0], param[1]));
    
    def out_func(xvals):
        return func(xvals, param[0], param[1]);
        
    return out_func;
    
    
def set_range(gamma, low, high):
    """Changes the range of gamma to stretch between low and high
    Used after finding the false-negative rate to raise the floor on the
    prediction of whether a gene is active

    gamma:  probability of expression in Genes x Samples
    low:    lower value of adjusted gamma, Genes x 1 or 1x1
    high:   upper value of adjusted gamma, Genes x 1 or 1x1
    
    """
    
    if(gamma.shape[0] > 1):
        if(type(low) is not np.ndarray):
            low = np.ones(gamma.shape[0])*low;
        if(type(high) is not np.ndarray):
            high = np.ones(gamma.shape[0])*high;
    
    high = high.flatten();
    low = low.flatten();    
    
    #transpose array so we can use matrix / vector operations
    gamma_out = gamma.copy().T;    
    
    old_min = np.min(gamma_out, axis=0);
    old_max = np.max(gamma_out, axis=0);
    
    gamma_out = gamma_out - old_min;
    gamma_out = gamma_out/old_max * (high-low);
    gamma_out = gamma_out + low;
    
    return gamma_out.T;
    

def plot_em_norm_distribution(gamma, mu_l, mu_h, st_l, st_h, data, i):

    mu_lx = mu_l[i];
    mu_hx = mu_h[i];
    st_hx = st_h[i];
    
    domain = np.linspace(0,10,10000);
    p_low = np.exp(-1*domain/mu_lx)/mu_lx;
    p_low[isnan(p_low)] = 0;    
    
    p_high = np.exp(-1 * (domain - mu_hx)**2 / (2*st_hx**2)) / st_hx / np.sqrt(2*np.pi);
    
    
    hold(False);
    (n, bins, patches) = hist(data[i,:], range=(0,10),bins=100, normed=True);
    hold(True);    
    plot(domain, p_low, color='green');
    plot(domain, p_high, color='green');
    scatter(data[i,:], gamma[i,:], color='red');
    ylim(0, 1.1);
    


#Could be used in make_monotonic.  However, I'm avoiding reliance on numba
#@jit   
#def find_max(gamma):
#    max_vals = np.zeros(gamma.shape[0]);
#    max_j = np.zeros(gamma.shape[0]);    
#    nR = gamma.shape[0];
#    nC = gamma.shape[1];    
#    
#    for i in range(nR):
#        for j in range(nC):
#            if(j == 0):
#                max_vals[i] = gamma[i,j];
#                max_j[i] = j;
#            else:
#                if(gamma[i,j] > max_vals[i]):
#                    max_vals[i] = gamma[i,j];
#                    max_j[i] = j;
#                    
#    return max_j;

def do_the_thing(i, cutoff):
    #cutoff = 5;
    #i = 1;
    vals = data[i,:];
    vals = vals[vals != 0];
    
    (gamma, mu_l, mu_h, st_l, st_h, Pi, L) = em.em_exp_norm_mixture(vals,cutoff);
    mu_l.shape = 1;
    mu_h.shape = 1;
    st_h.shape = 1;
    
    domain = np.linspace(0,10,10000);
    p_low = exp(-1*domain/mu_l)/mu_l;
    
    p_high = exp(-1 * (domain - mu_h)**2 / (2*st_h**2)) / st_h / np.sqrt(2*np.pi);
    
    
    f = figure();
    (n, bins, patches) = hist(vals, range=(0,10),bins=100, normed=True);
    plot(domain, p_low, color='green');
    plot(domain, p_high, color='green');
    scatter(vals, gamma, color='red');
    ylim(0, 1.1);
    
    display(f)