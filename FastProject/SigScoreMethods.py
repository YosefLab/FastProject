# -*- coding: utf-8 -*-
""" Different ways to evlauate the signature score

Each method should have the same signature so they can be swapped

Right now, each takes in a wrapped data object and a signature object

Specified by the command-line argument: --sig_score_method

"""

from .Signatures import SignatureScores;
import numpy as np;


def naive_eval_signature(self, signature, zeros, min_signature_genes):
    """
    Naive eval signature just sums the columns * sign
    Equivalent to all weights = 1
    """
    sig_vector = signature.sig_indices(self.row_labels);
    ii = np.nonzero(sig_vector)[0];

    num_matched_genes = ii.size;

    if(num_matched_genes == 0):
        raise ValueError("No genes match signature");

    if(num_matched_genes < min_signature_genes):
        raise ValueError("Too few genes match signature");

    data = self.base[ii, :];
    weights = np.ones(data.shape);
    sig_vector = sig_vector[ii, :];

    pdata = data * sig_vector * weights;

    sig_scores = pdata.sum(axis=0);
    sig_scores /= np.sum(np.abs(sig_vector) * weights, axis=0);  # Only normalize by weights in the signature

    sig_obj = SignatureScores(sig_scores, signature.name, self.col_labels,
                              isFactor=False, isPrecomputed=False, numGenes=num_matched_genes);

    return sig_obj;


def weighted_eval_signature(self, signature, zeros, min_signature_genes):
    """
    Weighted average using weights in self.weights
    """
    sig_vector = signature.sig_indices(self.row_labels);
    ii = np.nonzero(sig_vector)[0];

    num_matched_genes = ii.size;

    if(num_matched_genes == 0):
        raise ValueError("No genes match signature");

    if(num_matched_genes < min_signature_genes):
        raise ValueError("Too few genes match signature");

    weights = self.weights[ii, :];
    data = self.base[ii, :];
    sig_vector = sig_vector[ii, :];

    pdata = data * sig_vector * weights;

    sig_scores = pdata.sum(axis=0);
    norm_factor = np.sum(np.abs(sig_vector) * weights, axis=0);  # Only normalize by weights in the signature
    norm_factor[norm_factor == 0] = 1.0;
    sig_scores /= norm_factor;

    sig_obj = SignatureScores(sig_scores, signature.name, self.col_labels,
                              isFactor=False, isPrecomputed=False, numGenes=num_matched_genes);

    return sig_obj;


def imputed_eval_signature(self, signature, zeros, min_signature_genes):
    """
    Imputes likely values for zeros using the weight
    Weights represent p(not expressed | not detected)
    """
    sig_vector = signature.sig_indices(self.row_labels);
    ii = np.nonzero(sig_vector)[0];

    num_matched_genes = ii.size;

    if(num_matched_genes == 0):
        raise ValueError("No genes match signature");

    if(num_matched_genes < min_signature_genes):
        raise ValueError("Too few genes match signature");

    weights = self.weights[ii, :];
    data = self.base[ii, :];
    sig_vector = sig_vector[ii, :];
    zeros = zeros[ii, :];

    # impute values
    # weight represents p(not e | not d)
    mu = (data * (~zeros)).sum(axis=1, keepdims=True) / (~zeros).sum(axis=1, keepdims=True);
    mu[np.isnan(mu)] = 0.0;  # All zeros, mu is zero

    pe_nd = 1 - weights;  # Probability expressed | not detected
    pne_nd = weights;

    sig_scores_real = ((data * ~zeros) * sig_vector).sum(axis=0);
    sig_scores_imputed = np.sum(pe_nd * mu * sig_vector * zeros + pne_nd * data * sig_vector * zeros, axis=0);
    sig_scores = sig_scores_real + sig_scores_imputed;

    sig_obj = SignatureScores(sig_scores, signature.name, self.col_labels,
                              isFactor=False, isPrecomputed=False, numGenes=num_matched_genes);

    return sig_obj;


def nonzero_eval_signature(self, signature, zeros, min_signature_genes):
    """
    Evalute the signature scores using only nonzero entries
    Sum columns and normalize by # of nonzero items
    """
    sig_vector = signature.sig_indices(self.row_labels);
    ii = np.nonzero(sig_vector)[0];

    num_matched_genes = ii.size;

    if(num_matched_genes == 0):
        raise ValueError("No genes match signature");

    if(num_matched_genes < min_signature_genes):
        raise ValueError("Too few genes match signature");

    data = self.base[ii, :];
    sig_vector = sig_vector[ii, :];
    zeros = zeros[ii, :];

    sig_scores_real = ((data * ~zeros) * sig_vector).sum(axis=0);
    norm_factor = (~zeros).sum(axis=0);
    norm_factor[norm_factor == 0] = 1.0
    sig_scores = sig_scores_real / norm_factor;

    sig_obj = SignatureScores(sig_scores, signature.name, self.col_labels,
                              isFactor=False, isPrecomputed=False, numGenes=num_matched_genes);

    return sig_obj;
