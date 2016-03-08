# -*- coding: utf-8 -*-
"""EM Algorithm

Used to fit data to an exponential/normal mixture of
distributions

Derived from Nir Yosef's MATLAB code

"""
from __future__ import absolute_import, print_function, division;
import numpy as np;
from . import ProgressBar

#zmat is Genes x Samples
#cutoff is Genes x 1 - represents the initial cutoff between normal and exponential distributions

def em_exp_norm_mixture(zmat, cutoff, progressbar = True):
    with np.errstate(divide='ignore', invalid='ignore'):  #nans and infs are handled explicitly, don't need warning
        max_iter = 150;

        if(progressbar): pbar = ProgressBar(max_iter);

        #promote to 2d if single gene given
        if(zmat.ndim == 1):
            zmat = zmat.reshape((1, zmat.shape[0]));
            cutoff = np.array([cutoff]);

        #Make sure cutoff is 2d
        cutoff = cutoff.reshape((cutoff.shape[0],1));

        cutoffs = np.tile(cutoff, (1,zmat.shape[1]));
        gamma = zmat > cutoffs;
        Pi = np.mean(gamma,axis=1).reshape((zmat.shape[0],1));
        Pi[Pi == 0] = 1 / zmat.shape[1];

        mu_l=weighted_mean(zmat,1-gamma);
        mu_h=weighted_mean(zmat,gamma);
        st_h=weighted_std(zmat,gamma, mu_h);

        for niter in range(max_iter):
            #E
            prev_gamma = gamma;
            p_low = np.exp(-1*zmat/mu_l)/mu_l;

            p_high = np.exp(-1 * (zmat - mu_h)**2 / (2*st_h**2)) / st_h / np.sqrt(2*np.pi);

            p_low[~np.isfinite(p_low)] = 1e-5;
            p_low[p_low < 1e-5] = 1e-5;

            p_high[~np.isfinite(p_high)] = 1e-5;
            p_high[p_high<1e-5]   = 1e-5;

            gamma = (Pi * p_high) / ( ((1-Pi)*p_low) + (Pi*p_high));

            #M
            Pi = np.mean(gamma,axis=1).reshape((zmat.shape[0],1));
            mu_l=weighted_mean(zmat,1-gamma);
            mu_h=weighted_mean(zmat,gamma);
            st_h=weighted_std(zmat,gamma, mu_h);


            if(progressbar): pbar.update();

            if(niter % 10 == 0):
                biggest_change = np.max(np.abs(gamma - prev_gamma));
                if(biggest_change < 0.01):
                    break;


        #if niter == 1: print mu_l, mu_h, st_l, st_h
        #print 'Iteration: ', niter, ' L: ', sum(L);
        #if d>0.95:
        #break;
        if(progressbar): pbar.complete();

        L = np.sum(gamma * np.log(p_high) + (1-gamma)*np.log(p_low),axis=1);
        st_l=weighted_std(zmat,1-gamma, mu_l);

        #Final E
        prev_gamma = gamma;
        p_low = np.exp(-1*zmat/mu_l)/mu_l;

        p_high = np.exp(-1 * (zmat - mu_h)**2 / (2*st_h**2)) / st_h / np.sqrt(2*np.pi);

        p_low[~np.isfinite(p_low)] = 1e-5;
        p_low[p_low < 1e-5] = 1e-5;

        p_high[~np.isfinite(p_high)] = 1e-5;
        p_high[p_high<1e-5]   = 1e-5;

        gamma = (Pi * p_high) / ( ((1-Pi)*p_low) + (Pi*p_high));

        return (gamma, mu_l, mu_h, st_l, st_h, Pi, L);

#y and gamma same size
#returns average of each row, weighted by gammas
def weighted_mean(y, gamma):
    mu = np.sum(gamma * y, axis=1) / np.sum(gamma, axis=1);
    mu = mu.reshape((mu.shape[0], 1));  #make it 2 dimensional
    return mu

#y and gamma same size
#returns std dev of each row, weighted by gammas
def weighted_std(y, gamma, mu=None):
    if(mu is None):
        mu = weighted_mean(y,gamma);

    wstd = np.sqrt(np.sum((gamma * (y-mu)**2), axis=1) / np.sum(gamma, axis=1));
    wstd = wstd.reshape((wstd.shape[0],1));
    return wstd
