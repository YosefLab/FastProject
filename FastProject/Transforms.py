# -*- coding: utf-8 -*-
"""
Created on Wed Jan 07 15:12:02 2015

@author: David
"""
from __future__ import print_function;


from .Utils import em_exp_norm_mixture;
from . import Filters;
import numpy as np;
import os;

this_directory = os.path.dirname(os.path.abspath(__file__));

#from numba import jit;



def probability_of_expression(data, nozero=True):
    cutoffs = np.mean(data,axis=1)/4;  #Empirically found to be good mosy of the time
    
    (gamma, mu_l, mu_h, st_l, st_h, Pi, L) = em_exp_norm_mixture(data,cutoffs);
        
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
    
    Creates a functional fit for each sample based on that samples HK genes

    Returns
    ----------
    fit_func : function 
        Used to fit expression values to FN rate
    params : (Num_Params x Num_Samples) numpy.ndarray
        Sample-specific parameters to use with fit_func

    """
    
    housekeeping_dir = os.path.join(this_directory,'Housekeeping Genes');
    
    housekeeping_files = list();
    
    if(housekeeping_file != ""): #If file specified, use that file
        housekeeping_files.append(housekeeping_file);
    else:  #Otherwise, use all the files in housekeeping directory
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
    (gamma, mu_l, mu_h, st_l, st_h, Pi, L) = em_exp_norm_mixture(data_hk,cutoffs);
    
    #create bins based on average gene expression in population
#    bin_start = 0.0;
#    bin_end = np.ceil(np.max(mu_h));
#    
#    bin_starts = np.arange(bin_start, bin_end);
#    bin_ends = bin_starts + 1;
#    bin_ends[-1] = bin_ends[-1] + .1;
#    
#    binned_gammas = np.zeros((bin_starts.shape[0], gamma.shape[1]));
#    
#    for i in np.arange(binned_gammas.shape[0]):
#        start = bin_starts[i];
#        end = bin_ends[i];
#        
#        indices = np.logical_and(mu_h >= start, mu_h < end).flatten();
#        
#        if(indices.any()):
#            binned_gammas[i,:] = gamma[indices,:].mean(axis=0);
#        
#    
        
    #Fit a function mapping mu to gammas

    from scipy.optimize import curve_fit;
    def func(xvals, x0, a):
        return 1/(1 + np.exp((xvals-x0)*a));

    params = np.zeros((2,gamma.shape[1]));
    x = mu_h.flatten();

    count_fails = 0;
    for i in range(gamma.shape[1]):

        y = 1-gamma[:,i]

        try:
            param, cov = curve_fit(func, x, y);
        except RuntimeError:  #This occurs if the optimizer can't fit the function
            param = np.array([-5,5]);  
            #If this occurs, these parameters will just return zero for every
            #possible mu value.  In a sense, it will disable the false-negative
            #correction for that sample.
            count_fails = count_fails + 1;
            
        params[:,i] = param;
    
    if(count_fails > 0):
        print();
        print("Warning:  Could not fit False-Negative correction function for " + str(count_fails) + " samples");
        print("          False-Negative correction disabled for these samples");
        print();
#        #Uncomment to plot result
#        from pylab import figure, scatter, plot, ion;
#        ion();
#        figure();
#        domain = np.linspace(0,10,1000);
#        scatter(x,y);
#        plot(domain, func(domain, param[0], param[1]));
#        
        
    return func, params;

    
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
    p_low[np.isnan(p_low)] = 0;    
    
    p_high = np.exp(-1 * (domain - mu_hx)**2 / (2*st_hx**2)) / st_hx / np.sqrt(2*np.pi);
    
    from matplotlib.pyplot import hold, hist, plot, scatter, ylim    
    
    hold(False);
    (n, bins, patches) = hist(data[i,:], range=(0,10),bins=100, normed=True);
    hold(True);    
    plot(domain, p_low, color='green');
    plot(domain, p_high, color='green');
    scatter(data[i,:], gamma[i,:], color='red');
    ylim(0, 1.1);
    
def probability_transform(data, original_data, original_genes, housekeeping_filename):
    """ Process data to evaluate probability of expression """
    
    print()
    print('Fitting expression data to exp/norm mixture model');
    (prob, mu_h) = probability_of_expression(data);
    prob = make_monotonic(prob, data);
    
    print();
    print('Correcting for false-negatives using housekeeping gene levels');
    (fit_func, params) = create_false_neg_map(original_data, original_genes, housekeeping_filename);
    
    #fit_func is the fitting function of the form fit_func(mu_h, param[0], param[1], etc)
    #params is a matrix with parameters for each sample.  Size is Num_Params x Num_samples
    
    fn_prob = np.zeros(prob.shape)
    
    for i in range(data.shape[1]):
        fn_prob[:,i] = fit_func(mu_h, params[0,i], params[1,i]).ravel();
    
    prob2 = prob + (1-prob)*fn_prob;
    
    return prob2, fn_prob

#def utility_plotting_routine(i, cutoff):
#    #cutoff = 5;
#    #i = 1;
#    vals = data[i,:];
#    vals = vals[vals != 0];
#    
#    (gamma, mu_l, mu_h, st_l, st_h, Pi, L) = em.em_exp_norm_mixture(vals,cutoff);
#    mu_l.shape = 1;
#    mu_h.shape = 1;
#    st_h.shape = 1;
#    
#    domain = np.linspace(0,10,10000);
#    p_low = exp(-1*domain/mu_l)/mu_l;
#    
#    p_high = exp(-1 * (domain - mu_h)**2 / (2*st_h**2)) / st_h / np.sqrt(2*np.pi);
#    
#    
#    f = figure();
#    (n, bins, patches) = hist(vals, range=(0,10),bins=100, normed=True);
#    plot(domain, p_low, color='green');
#    plot(domain, p_high, color='green');
#    scatter(vals, gamma, color='red');
#    ylim(0, 1.1);
#    
#    display(f)