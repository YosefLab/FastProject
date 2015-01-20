# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 12:29:19 2014

@author: daved_000
"""

from numpy import *

#zmat is Genes x Samples
#cutoff is Genes x 1 - represents the initial cutoff between normal and exponential distributions

def em_exp_norm_mixture(zmat, cutoff):
    #promote to 2d if single gene given
    if(zmat.ndim == 1):
        zmat = zmat.reshape(1, zmat.shape[0]);
        cutoff = array([cutoff]);
    
    cutoff = cutoff.reshape(cutoff.shape[0],1); 
    cutoffs = tile(cutoff, (1,zmat.shape[1]));
    gamma = zmat > cutoffs;
    Pi = mean(gamma,1).reshape((zmat.shape[0],1));
    
    mu_l=weighted_mean(zmat,1-gamma);
    mu_h=weighted_mean(zmat,gamma);
    st_h=weighted_std(zmat,gamma);
    st_l=weighted_std(zmat,1-gamma);
    st_l[isnan(st_l)]=1e-3;
    
    to_stop = zeros((zmat.shape[0],1)) > 1;
    
    p_low = zeros(zmat.shape);
    p_high = zeros(zmat.shape);    
    
    for niter in arange(0,150)+1:
        gamma_prev = gamma;
        #E
        mu_l2 = tile(mu_l, (1,zmat.shape[1]));
        p_low = exp(-1*zmat/mu_l2)/mu_l2;
        
        mu_h2 = tile(mu_h, (1,zmat.shape[1]));
        st_h2 = tile(st_h, (1,zmat.shape[1]));
        p_high = exp(-1 * (zmat - mu_h2)**2 / (2*st_h2**2)) / st_h2 / sqrt(2*pi);
        
        p_low[isnan(p_low)] = 1e-5;
        p_low[p_low < 1e-5] = 1e-5;
        
        p_high[isnan(p_high)] = 1e-5;
        p_high[p_high<1e-5]   = 1e-5;
        
        Pix = tile(Pi,(1,zmat.shape[1]));
        
        gamma = (Pix * p_high) / ( ((1-Pix)*p_low) + (Pix*p_high));
        #gamma[to_stop] = gamma_prev[to_stop]; #positions marked (done) are frozen

        #M
        Pi = mean(gamma,1).reshape((zmat.shape[0],1));
        mu_l=weighted_mean(zmat,1-gamma);
        mu_h=weighted_mean(zmat,gamma);
        st_h=weighted_std(zmat,gamma);
        st_l=weighted_std(zmat,1-gamma);
        st_l[isnan(st_l)]=1e-3;
        
        L = sum(gamma * log(p_high) + (1-gamma)*log(p_low),axis=1);
        
        #Stop
        to_stop = (gamma_prev > 0.5) == (gamma > 0.5);
        d=sum(to_stop) / float(size(to_stop));
        if niter == 1: print mu_l, mu_h, st_l, st_h    
        print 'Iteration: ', niter, ' L: ', sum(L);
        #if d>0.95:
            #break;

    return (gamma, mu_l, mu_h, st_l, st_h, Pi, L);


#y and gamma same size
#returns average of each row, weighted by gammas
def weighted_mean(y, gamma):
    mu = sum(gamma * y, axis=1) / sum(gamma, axis=1);
    mu = mu.reshape(mu.shape[0], 1);  #make it 2 dimensional    
    return mu

#y and gamma same size
#returns std dev of each row, weighted by gammas
def weighted_std(y, gamma):
    mu = weighted_mean(y,gamma);
    mu_tiled = tile(mu, (1,y.shape[1])); 
    wstd = sqrt(sum((gamma * (y-mu_tiled)**2), axis=1) / sum(gamma, axis=1));
    wstd = wstd.reshape(wstd.shape[0],1);    
    return wstd