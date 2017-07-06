# -*- coding: utf-8 -*-
""" This module handles the splitting and joining of the data if sub-sampling is enabled

This process is likely out of date - needs attention
before re-enabling this feature.
"""
from __future__ import absolute_import, print_function, division;
from .Utils import ProgressBar;
from . import Signatures;
from . import SigScoreMethods;
from . import Projections;
from . import Transforms;
from .DataTypes import ExpressionData;
from .Global import FP_Output, RANDOM_SEED;
from scipy.spatial.distance import cdist;
import numpy as np;

def split_samples(data, subsample_size):
    np.random.seed(RANDOM_SEED); #Set seed so outputs are repeatable
    sub_ii = np.random.choice(data.shape[1], subsample_size, replace=False);
    holdouts_ii = np.setxor1d(np.arange(data.shape[1]), sub_ii);

    subset = data.subset_samples(sub_ii);
    holdouts = data.subset_samples(holdouts_ii);

    return holdouts, subset;

def merge_samples(all_data, Models, sigs, prob_params, args):
    #Remove genes that aren't in the Expression data set

    edata = Models["Expression"]["Data"];
    keep_genes = set(edata.row_labels);
    keep_genes_i = [i for i,gene in enumerate(all_data.row_labels) if gene in keep_genes];
    edata_all = all_data.subset_genes(keep_genes_i);
    edata_all.filters = edata.filters;

    #First need to recompute prob info
    #TODO:  Ensure order of genes is the same as the ordering in mu, st, etc
    #M Step of EM Algorithm
    if(not args["nomodel"]):
        FP_Output("Extending probability model to holdouts");
        (mu_h, mu_l, st_h, Pi) = prob_params;
        zmat = edata_all.base;
        p_low = np.exp(-1*zmat/mu_l)/mu_l;

        p_high = np.exp(-1 * (zmat - mu_h)**2 / (2*st_h**2)) / st_h / np.sqrt(2*np.pi);

        p_low[~np.isfinite(p_low)] = 1e-5;
        p_low[p_low < 1e-5] = 1e-5;

        p_high[~np.isfinite(p_high)] = 1e-5;
        p_high[p_high<1e-5]   = 1e-5;

        gamma = (Pi * p_high) / ( ((1-Pi)*p_low) + (Pi*p_high));
        pdata_all = Transforms.make_monotonic(gamma, zmat);
        (fit_func, params) = Transforms.create_false_neg_map(edata_all, args["housekeeping"]);
        (pdata_all, fn_prob) = Transforms.correct_for_fn(pdata_all, mu_h, fit_func, params);
        fn_prob[edata_all > 0] = 0;

        edata_all.weights = 1-fn_prob;
        sample_passes_qc, sample_qc_scores = Transforms.quality_check(params);

        #Need to make sure we retain all the original labels
        sample_passes_labels = {label for label,val in zip(edata_all.col_labels, sample_passes_qc) if val};
        sample_passes_labels = sample_passes_labels.union(edata.col_labels);
        sample_passes_qc = np.array([i for i,x in enumerate(edata_all.col_labels) if x in sample_passes_labels]);

        #If specified, remove items that did not pass qc check
        if(args["qc"]):
                edata_all = edata_all.subset_samples(sample_passes_qc);

    Transforms.z_normalize(edata_all);

    model = Models["Expression"];
    sub_data = model["Data"];
    all_data = edata_all;

    sub_data_cols = set(sub_data.col_labels);
    all_cols = set(all_data.col_labels);
    holdout_cols = all_cols.difference(sub_data_cols);
    holdout_indices = [i for i,x in enumerate(all_data.col_labels) if x in holdout_cols];

    holdout_data = all_data.subset_samples(holdout_indices);
    model["Holdout_Data"] = holdout_data;

    #Merge in projections
    #Re calculate clusters
    for name, model in Models.items():
        sub_data = model["Data"];
        holdout_data = model["Holdout_Data"];

        FP_Output("Calculating projection coordinates and cluster assignments for Holdouts -", name);

        for projData in model["projectionData"]:
            #Merge projection back in
            projData["projections"] = merge_projections(projData["projections"], sub_data, holdout_data);

            #Re-cluster
            clusters = Projections.define_clusters(projData["projections"]);
            projData["clusters"] = clusters;


    #Then merge back in data
    FP_Output("Merging holdout data matrices back into sub-sample");
    if(not args["nomodel"]):
        model = Models["Probability"];
        merge_pdata = model["Data"].merge_data(model["Holdout_Data"]);
        merge_edata = merge_pdata.expression_data;
        Models["Expression"]["Data"] = merge_edata;
        Models["Probability"]["Data"] = merge_pdata;
    else:
        model = Models["Expression"];
        merge_edata = model["Data"].merge_data(model["Holdout_Data"]);
        Models["Expression"]["Data"] = merge_edata;


    #Some cleanup and fix the labels
    for name, model in Models.items():
        model.pop("Holdout_Data");
        model["sampleLabels"] = model["Data"].col_labels;


    #Merge sig scores back in by re-calculating them
    for name, model in Models.items():
        data = model["Data"];

        #Evaluate Signatures
        FP_Output("Calculating Signature Scores for Holdouts -", name);

        sig_scores_dict = dict();
        old_sig_scores_dict = model["signatureScores"];

        pbar = ProgressBar(len(old_sig_scores_dict));
        for sig in sigs:
            if(sig.name in old_sig_scores_dict):
                try:
                    sig_scores_dict[sig.name] = data.eval_signature(sig);
                except ValueError:  #Only thrown when the signature has no genes in the data
                    pass #Just discard the signature then
                pbar.update();
        pbar.complete();

        if(args["precomputed"]):
            for precomputed_file in args["precomputed"]:
                precomputed_sig_scores = Signatures.load_precomputed(precomputed_file, data.col_labels);
                sig_scores_dict.update(precomputed_sig_scores);

        #Adds in quality score as a pre-computed signature
        if(not args["nomodel"]):
            sig_scores_dict["FP_Quality"] = SigScoreMethods.SignatureScores(sample_qc_scores,"FP_Quality",data.col_labels,isFactor=False, isPrecomputed=True);

        model["signatureScores"] = sig_scores_dict;



def merge_projections(sub_projections, sub_data, holdout_data):
    """
    Finds approximate locations for held-out data in 2d projections

    :param sub_coordinates:
    :param sub_data:
    :param holdout_data:
    :return:
    """

    all_projections = dict();

    for key,p in sub_projections.items():
        K=10;
        subsample_dist = cdist(holdout_data.T, sub_data.T, metric='euclidean');
        #subsample_dist shape is not_sub_data N_Samples x sub_data N_Samples

        subsample_dist_arg = subsample_dist.argsort(axis=1);
        subsample_dist_arg = subsample_dist_arg[:, np.arange(K)]; #Take closest K

        x_coords = p[0,:][subsample_dist_arg];
        y_coords = p[1,:][subsample_dist_arg];

        #New X,Y coordinate is average of nearest neighbors
        holdout_x = np.mean(x_coords, axis=1);
        holdout_y = np.mean(y_coords, axis=1);

        NUM_SAMPLES = sub_data.shape[1] + holdout_data.shape[1];
        all_data = np.zeros((2,NUM_SAMPLES));
        all_data[0,:sub_data.shape[1]] = p[0,:];
        all_data[1,:sub_data.shape[1]] = p[1,:];
        all_data[0,sub_data.shape[1]:] = holdout_x;
        all_data[1,sub_data.shape[1]:] = holdout_y;
        all_projections[key] = all_data;

    return all_projections;


