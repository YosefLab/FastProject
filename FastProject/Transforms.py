# -*- coding: utf-8 -*-
"""
Created on Wed Jan 07 15:12:02 2015

@author: David
"""
from __future__ import print_function;


from .Utils import em_exp_norm_mixture;
from . import Filters;
from .DataTypes import ExpressionData, ProbabilityData, PCData;
import numpy as np;
import os;

this_directory = os.path.dirname(os.path.abspath(__file__));


def probability_of_expression(data, nozero=True):
    cutoffs = np.mean(data,axis=1)/4;  #Empirically found to be good mosy of the time
    
    (gamma, mu_l, mu_h, st_l, st_h, Pi, L) = em_exp_norm_mixture(data,cutoffs);
    
    gamma = make_monotonic(gamma, data);    
    
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
    

def create_false_neg_map(data, housekeeping_file=""):
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
        data_t = Filters.load_from_file(data, hkf);        
        data_hk = np.vstack((data_hk, data_t));
        genes_hk.extend(data_t.row_labels);

    data_hk = ExpressionData(data_hk, genes_hk, data.col_labels);

    data_hk = Filters.filter_genes_threshold(data_hk, 0.2);
        
    #calculate distributions for hk gene
    cutoffs = np.mean(data_hk,axis=1)/4;
    (gamma, mu_l, mu_h, st_l, st_h, Pi, L) = em_exp_norm_mixture(data_hk,cutoffs);
           
    #Fit a function mapping mu to gammas

    def func(xvals, x0, a, L, H):
        return L + (H-L)/(1 + np.exp((xvals-x0)*a));

    def efun(x,y, args):
        out = func(x, args[0], args[1], args[2], args[3]);
        return np.sum((out-y)**2);

    params = np.zeros((4,gamma.shape[1]));
    x = mu_h.flatten();
    
    sort_i = np.argsort(x);
    x_sorted = x[sort_i];
    
    initial_guess = [3.5, 1.26, 0, 1];
    bounds = [(0, np.inf),(0, 2),(0,1), (0,1)];

    if(len(x_sorted) > 30):
        q_indices = len(x_sorted)//31 * np.arange(31);
        q_indices = q_indices[1:];
    else:
        q_indices = np.arange(30);
    q_indices = q_indices.astype(np.int64);
    x_quant = x_sorted[q_indices];
    
    x_unique, q_i = np.unique(x_quant, return_index=True);
    q_indices = q_indices[q_i];

    from scipy.optimize import minimize;

    for i in range(gamma.shape[1]):
        y = 1-gamma[:,i]
        y_sorted = y[sort_i];

        y_unique = y_sorted[q_indices];

        res = minimize(lambda args: efun(x_unique,y_unique,args), initial_guess, bounds=bounds);
        param = res.x;         
        params[:,i] = param;
    
    

#        #Uncomment to plot result - useful for debugging
#        from pylab import figure, scatter, plot, ion;
#        ion();
#        figure();
#        domain = np.linspace(0,10,1000);
#        scatter(x_unique,1-gamma[sort_i,i][q_indices]);
#        plot(domain, func(domain, params[0,i], params[1,i], params[2,i], params[3,i]));
#        ylabel("P(gene not expressed in " + cells[i] + ")");
#        xlabel("Gene average in samples expressing gene")
#        print(params[:,i])
#        i = i+1;

        
    return func, params;

def quality_check(params):
    """Integrates the logistic false-negative curves.  Flags samples whose 
    integral is more than 1.6 MAD lower than the population.  
    
    Parameters
    ----------
    params : (4 x Num_Samples) numpy.ndarray 
        Matrix containing parameters for the false-negative fit function

    Returns
    -------
    sample_passes : (Num_Samples) boolean numpy.ndarray
        Vector containing True for samples that pass this quality check
          
    """
    
    #Logistic parameters
    x0 = params[0,:];
    a = params[1,:];
    L = params[2,:];
    H = params[3,:];
    
    #Bounds of integration
    low = 0;
    high = 9;
    
    a[a == 0] = 1e-6;  #Fix so that integral isn't mis-calculated as inf or nan
    
    #Evaluate integral
    int_low = H*low   - (H-L)/a * np.log(np.exp(a*(low -x0)) + 1)
    int_high = H*high - (H-L)/a * np.log(np.exp(a*(high-x0)) + 1)
    
    int_val = int_high - int_low;
    int_val_mean = np.mean(int_val);    
    
    abs_dev = np.abs(int_val - int_val_mean);
    
    MAD = np.mean(abs_dev);
    
    sample_passes = int_val < (int_val_mean + 1.6*MAD);
    
    return sample_passes;
    

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

    
def correct_for_fn(prob, mu_h, fit_func, params):
    """
    Uses the estimated false_negative, fn(mu) curves to correct probability values.
        
    Parameters
    ----------
    prob : (Num_Genes x Num_Samples) numpy.ndarray
        Matrix containing estimate for probability of expression of each gene in each sample
    mu_h : (Num_Genes x 1) numpy.ndarray
        Average expression value of each gene across all samples in which gene is expressed
    fit_func : function (mu_h, params)
        Function, parameterized by params, that maps each mu_h to a false negative estimate 
    params : (4 x Num_Samples) numpy.ndarray 
        Matrix containing parameters for the false-negative fit function (fit_func)

    Returns
    -------
    out_prob : (Num_Genes x Num_Samples) numpy.ndarray     
        Adjusted probability values
    fn_prob : (Num_Genes x Num_Samples) numpy.ndarray
        Estimated False-negative probability for each gene in each sample
        
    """

    fn_prob = np.zeros(prob.shape)
    
    for i in range(prob.shape[1]):
        fn_prob[:,i] = fit_func(mu_h, *params[:,i]).ravel();
    
    out_prob = prob + (1-prob)*fn_prob;
    
    return out_prob, fn_prob


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

def z_normalize(data):
    """
    Z-normalizes the rows of the matrix in data.
    No return as operation is done in place.
    Mean is subtracted out result is scaled so standard deviation = 1
    :param data: numpy.ndarray, 2 dimensions
    :return: None
    """
    if(data is ProbabilityData or data is PCData):
        raise TypeError("Should not be z-normalizing Probability or PCData, exiting");

    mu = data.mean(axis=1, keepdims=True);
    sigma = data.std(axis=1, keepdims=True);

    #Shouldn't be necessary for expression data that's been filtered, but good practice
    sigma[sigma == 0] = 1;

    #Note, operations are in place, no return needed
    data -= mu;
    data /= sigma;
