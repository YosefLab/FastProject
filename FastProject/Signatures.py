# -*- coding: utf-8 -*-
"""
Created on Wed Jan 21 16:34:49 2015

@author: David
"""
from __future__ import division, print_function;

import numpy as np;
from sklearn.metrics.pairwise import pairwise_distances;
from .Utils import ProgressBar;

def read_signatures(filename='', match_terms=[]):
    """Calls either read_signatures_txt or read_signatures_gmt when appropriate
    """

    if(filename.lower().endswith('.gmt')):
        return read_signatures_gmt(filename, match_terms);
    else:
        return read_signatures_txt(filename, match_terms);


def read_signatures_txt(filename='', match_terms=[]):
    """Reads in a signature file.  Returns a list of Signature objects
    
    Parameters
    ----------
    filename : string
        Name (and path if not in working directory) of signature file to read
        If no file is entered, opens a file dialog
    match_terms : list(String)
        List of terms to be matched against signature names.
        If empty or omitted, reads all signatures
        Signature is retained if any of the terms match.
    
    Returns
    -------
    signatures : list(FastProject.Signature)
        The filtered signature list
    """
    
    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();
    
    ff = open(filename, 'r');
    found_signatures = dict();    ## dictionary of signatures
    
    if(type(match_terms) is str):
        match_terms = [match_terms];
    
    match_terms = [term.lower() for term in match_terms];
    
    try:    

        for i, line in enumerate(ff):
            row_data = line.strip().split('\t');
            if(len(row_data) == 2):
                signed = False;
            elif(len(row_data) == 3):
                signed = True;
            else:
                raise ValueError("Line " + str(i) + " Signature file should contain 2 (unsigned) or 3 (signed) columns");

            name = row_data[0];
            if(not found_signatures.has_key(name)):  ## Add the signature if we haven't seen it yet

                ## Only add signatures if they match one of the supplied terms
                matched = False if len(match_terms) > 0 else True
                lname = name.lower();                
                for term in match_terms:
                    if(lname.find(term) >= 0):
                        matched = True;
                        break;
                
                if(matched):
                    sig = Signature(dict(), signed, filename, name);
                    found_signatures.update({name: sig});
                else:
                    continue;
            
            sig_dict = found_signatures[name].sig_dict;
            
            
            if signed:
                sig_sign = row_data[1].lower();
                if(sig_sign == 'plus'):
                    sig_val = 1;
                elif(sig_sign == 'minus'):
                    sig_val = -1;
                elif(sig_sign == 'both'):
                    sig_val = 0;
                elif(sig_sign == 'mius'):  #LOL, close enough
                    sig_val = -1;
                ## Other spelling variants go here
                else:
                    print("Error on line ", str(i), " Couldn't read signature value.");
                    print("   :", line);
                    continue;
                    
                sig_dict[row_data[2]] = sig_val;
                
            else:  #unsigned case
                sig_val = 0;
                sig_dict[row_data[1]] = sig_val;
     
    except:
        raise;
    finally:
        ff.close();
    
    return [found_signatures[key] for key in found_signatures.keys()]  #dict to list

def read_signatures_gmt(filename='', match_terms=[]):
    """Reads in a signature file in GMT format.

    Each row is tab delimited with gene set name, description, and list of genes
    Signed signatures should be on adjacent rows with identical set names except for
    a _plus or _minus suffix.

    Example:
    Memory_plus description gene1 gene2 gene3
    Memory_minus description gene1 gene2 gene3

    Returns a list of Signature objects

    Parameters
    ----------
    filename : string
        Name (and path if not in working directory) of signature file to read
        If no file is entered, opens a file dialog
    match_terms : list(String)
        List of terms to be matched against signature names.
        If empty or omitted, reads all signatures
        Signature is retained if any of the terms match.

    Returns
    -------
    signatures : list(FastProject.Signature)
        The filtered signature list
    """

    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();

    if(type(match_terms) is str):
        match_terms = [match_terms];

    match_terms = [term.lower() for term in match_terms];

    with open(filename, 'r') as ff:
        found_signatures = dict();    ## dictionary of signatures

        end_of_file = False;
        while(not end_of_file):
            line = ff.readline();
            if(line == ""):
                end_of_file = True;
                continue;

            row_data = line.strip().split('\t');
            name = row_data[0];

            ## Only add signatures if they match one of the supplied terms
            matched = False if len(match_terms) > 0 else True
            lname = name.lower().rstrip("_minus").rstrip("_plus");
            for term in match_terms:
                if(term in lname):
                    matched = True;
                    break;
            if(not matched): continue;

            if(name.lower().endswith("_plus") or name.lower().endswith("_minus")):
                #Signature is first of a signed pair
                row2_data = ff.readline().strip().split('\t');
                name2 = row2_data[0];
                if(name.lower().endswith("_plus") and name2.lower().endswith("_minus")):
                    plus_row = row_data;
                    minus_row = row2_data;
                elif(name2.lower().endswith("_plus") and name.lower().endswith("_minus")):
                    minus_row = row_data;
                    plus_row = row2_data;
                else:
                    raise ValueError("Missing pair for signed signature: " + name);

                sig_name = plus_row[0];
                ii = sig_name.lower().rfind("_plus");
                sig_name = sig_name[:ii]
                sig = Signature(dict(), True, filename,sig_name)

                for gene_name in plus_row[2:]:
                    sig.sig_dict.update({gene_name: 1});
                for gene_name in minus_row[2:]:
                    sig.sig_dict.update({gene_name: -1});

                found_signatures.update({name: sig});

            else:
                if(found_signatures.has_key(name)):
                    raise ValueError("Duplicate Signature name in Signature file: " + name);

                sig = Signature(dict(), False, filename, name);
                for gene_name in row_data[2:]:
                    sig.sig_dict.update({gene_name: 0});
                found_signatures.update({name: sig});


    return found_signatures.values();  #dict to list

def filter_sig_list(signatures, match_terms):
    """
    Filters the list of signature objects to retain only those that contain
    one of the specified terms
    
    Parameters
    ----------
    signatures : list(Signature)
        List of signatures to be filtered
    match_terms : list(String)
        List of terms to be matched against signature names.
        Signature is retained if any of the terms match.
    
    Returns
    -------
    filtered_signatures : list(FastProject.Signature)
        The filtered signature list
        
    
    """
    
    if(type(match_terms) is str):       #So you dont need to wrap a single string in []
        match_terms = [match_terms];    
    
    match_terms = map(lambda term: term.lower(), match_terms);    
    
    filtered_signatures = list();
    for sig in signatures:
        name = sig.name.lower();
        
        for term in match_terms:            
            if(name.find(term) >= 0):
                filtered_signatures.append(sig);
                break;
    
    
    return filtered_signatures;


def conformity_with_p(data_loc, sig_values, NEIGHBORHOOD_SIZE = 0.1):
    """
    Score each sample based on how similar its signature score is to its neighborhood
    Then compare similarity to shuffled data and derive a p-value
    
    Parameters
    ---------
    data_loc : array-like, shape (Num_Dimensions, Num_Samples)
        locations of points used to calculate neareset neighbors
    sig_values : array-like, 1D, shape (Num_Samples)
        Signature value for each sample.  Get using Signature.eval_data
    NEIGHBORHOOD_SIZE : float
        Approximate size of the neighborhood to consider around points

    Returns
    -------
    dissimilarity : array-like, 1D, shape (Num_Samples)
        Score for each sample as to how dissimilar it is from neighboring points
    p_value : float
        Probability that the scores are significant
        
    """
    
    distance_matrix = pairwise_distances(data_loc.T, metric='euclidean');

    weights = np.exp(-1 * distance_matrix**2 / NEIGHBORHOOD_SIZE);
    np.fill_diagonal(weights,0); #Don't count self

    neighborhood = sig_values.reshape((1,len(sig_values)));

    sum_weights = np.sum(weights, axis=1);
    neighborhood_prediction = np.sum(neighborhood * weights, axis=1) \
                / sum_weights;
    
    
    ##Neighborhood dissimilarity score = |actual - predicted|
    dissimilarity = np.abs(sig_values - neighborhood_prediction);

    #Minimum of 100 iterations.  Add in more if needed.
    NUM_RAND_TRIALS_L1 = 100
    random_dissimilarity = np.zeros(NUM_RAND_TRIALS_L1);
    
    for i in range(NUM_RAND_TRIALS_L1):
        random_sig_values = np.random.permutation(sig_values);
        random_neighborhood = random_sig_values.reshape((1,len(random_sig_values)));

        neighborhood_prediction = np.sum(random_neighborhood * weights, axis=1) \
                / sum_weights;

        random_dissimilarity[i] = np.median(
                                    np.abs(
        random_sig_values - neighborhood_prediction
        ));
        
    #Compare number of times the random permutation has a better (median) sig score
    #than the signature score of non-shuffled.  Ratio is p value.
    count_random_wins = np.count_nonzero(random_dissimilarity <= np.median(dissimilarity));
    p_value = (1 + count_random_wins) / (1 + NUM_RAND_TRIALS_L1);

    if(p_value < 0.05):

        NUM_RAND_TRIALS_L2 = 100
        random_dissimilarity = np.zeros(NUM_RAND_TRIALS_L2);

        for i in range(NUM_RAND_TRIALS_L2):
            random_sig_values = np.random.permutation(sig_values);
            random_neighborhood = random_sig_values.reshape((1,len(random_sig_values)));

            neighborhood_prediction = np.sum(random_neighborhood * weights, axis=1) \
                                      / sum_weights;

            random_dissimilarity[i] = np.median(
                np.abs(
                    random_sig_values - neighborhood_prediction
                ));

        #Compare number of times the random permutation has a better (median) sig score
        #than the signature score of non-shuffled.  Ratio is p value.
        count_random_wins += np.count_nonzero(random_dissimilarity <= np.median(dissimilarity));
        p_value = (1 + count_random_wins) / (1 + NUM_RAND_TRIALS_L1 + NUM_RAND_TRIALS_L2);


    if(p_value < 0.01):

        NUM_RAND_TRIALS_L3 = 500
        random_dissimilarity = np.zeros(NUM_RAND_TRIALS_L3);

        for i in range(NUM_RAND_TRIALS_L3):
            random_sig_values = np.random.permutation(sig_values);
            random_neighborhood = random_sig_values.reshape((1,len(random_sig_values)));

            neighborhood_prediction = np.sum(random_neighborhood * weights, axis=1) \
                                      / sum_weights;

            random_dissimilarity[i] = np.median(
                np.abs(
                    random_sig_values - neighborhood_prediction
                ));

        #Compare number of times the random permutation has a better (median) sig score
        #than the signature score of non-shuffled.  Ratio is p value.
        count_random_wins += np.count_nonzero(random_dissimilarity <= np.median(dissimilarity));
        p_value = (1 + count_random_wins) / (1 + NUM_RAND_TRIALS_L1 + NUM_RAND_TRIALS_L2 + NUM_RAND_TRIALS_L3);

    #Scale dissimilarity by the MAD
    #  this allows dissimilarity values to better be compared between signatures of varying magnitudes
    MAD = np.median(np.abs(sig_values - np.median(sig_values)));
    if(MAD != 0):
        dissimilarity /= MAD;
    
    return dissimilarity, p_value

def sigs_vs_projections(projections, sig_scores, NEIGHBORHOOD_SIZE = 0.1):
    """
    Evaluates the significance of each signature vs each projection

    :param projections: dict of (string) => (numpy.ndarray of shape 2xNum_Samples)
        Maps projections to their spatial coordinates for each sample
    :param sig_scores: dict of (string) => (numpy.ndarray of shape Num_Samples)
        Maps signature names to their value at each coordinate
    :return:
    """
    sp_row_labels = sig_scores.keys();
    sp_col_labels = projections.keys();

    N_SAMPLES = sig_scores[sp_row_labels[0]].shape[0];
    N_SIGNATURES = len(sp_row_labels);
    N_PROJECTIONS = len(sp_col_labels);


    sig_proj_matrix   = np.zeros((N_SIGNATURES,N_PROJECTIONS));
    sig_proj_matrix_p = np.zeros((N_SIGNATURES,N_PROJECTIONS));


    #Build a matrix of all signatures
    sig_score_matrix = np.zeros((N_SAMPLES, N_SIGNATURES));

    for j, sig in enumerate(sp_row_labels):
        sig_score_matrix[:,j] = sig_scores[sig];


    NUM_RAND_TRIALS_L1 = 500;
    print();
    print("Evaluating Signatures against Projections");
    pp = ProgressBar(N_PROJECTIONS*NUM_RAND_TRIALS_L1);
    for i, proj in enumerate(sp_col_labels):
        data_loc = projections[proj];

        distance_matrix = pairwise_distances(data_loc.T, metric='euclidean');

        weights = np.exp(-1 * distance_matrix**2 / NEIGHBORHOOD_SIZE);
        np.fill_diagonal(weights,0); #Don't count self
        weights /= np.sum(weights, axis=1, keepdims=True);

        neighborhood_prediction = np.dot(weights, sig_score_matrix);


        ##Neighborhood dissimilarity score = |actual - predicted|
        dissimilarity = np.abs(sig_score_matrix - neighborhood_prediction);
        med_dissimilarity = np.median(dissimilarity, axis=0);


        #Minimum of 100 iterations.  Add in more if needed.
        random_win_count = np.zeros(N_SIGNATURES);

        for j in range(NUM_RAND_TRIALS_L1):
            random_sig_values = np.random.permutation(sig_score_matrix); #Shuffles on first index, samples

            neighborhood_prediction = np.dot(weights, random_sig_values);

            med_random_dissimilarity = np.median(
                np.abs(
                    random_sig_values - neighborhood_prediction
                ), axis=0);

            random_win_count[med_random_dissimilarity <= med_dissimilarity] += 1;
            pp.update();

        p_values = (1 + random_win_count) / (1 + NUM_RAND_TRIALS_L1);

        sig_proj_matrix[:,i] = med_dissimilarity;
        sig_proj_matrix_p[:,i] = p_values;


    pp.complete();

    MAD = np.median(np.abs(sig_score_matrix - np.median(sig_score_matrix, axis=0)), axis=0);
    MAD[MAD==0] = 1;
    MAD = MAD.reshape((N_SIGNATURES, 1));
    sig_proj_matrix /= MAD;

    return (sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p);

def sigs_vs_projections_v2(projections, sig_scores, NEIGHBORHOOD_SIZE = 0.1):
    """
    Evaluates the significance of each signature vs each projection

    :param projections: dict of (string) => (numpy.ndarray of shape 2xNum_Samples)
        Maps projections to their spatial coordinates for each sample
    :param sig_scores: dict of (string) => (numpy.ndarray of shape Num_Samples)
        Maps signature names to their value at each coordinate
    :return:
    """
    sp_row_labels = sig_scores.keys();
    sp_col_labels = projections.keys();
    sp_col_labels.sort();

    N_SAMPLES = sig_scores[sp_row_labels[0]].shape[0];
    N_SIGNATURES = len(sp_row_labels);
    N_PROJECTIONS = len(sp_col_labels);


    sig_proj_matrix   = np.zeros((N_SIGNATURES,N_PROJECTIONS));
    sig_proj_matrix_p = np.zeros((N_SIGNATURES,N_PROJECTIONS));


    #Build a matrix of all signatures
    sig_score_matrix = np.zeros((N_SAMPLES, N_SIGNATURES));

    for j, sig in enumerate(sp_row_labels):
        sig_score_matrix[:,j] = sig_scores[sig];


    print();
    print("Evaluating Signatures against Projections");
    pp = ProgressBar(N_PROJECTIONS);
    for i, proj in enumerate(sp_col_labels):
        data_loc = projections[proj];

        distance_matrix = pairwise_distances(data_loc.T, metric='euclidean');

        weights = np.exp(-1 * distance_matrix**2 / NEIGHBORHOOD_SIZE);
        np.fill_diagonal(weights,0); #Don't count self
        weights /= np.sum(weights, axis=1, keepdims=True);

        neighborhood_prediction = np.dot(weights, sig_score_matrix);


        ##Neighborhood dissimilarity score = |actual - predicted|
        dissimilarity = np.abs(sig_score_matrix - neighborhood_prediction);
        med_dissimilarity = np.median(dissimilarity, axis=0);

        NUM_REPLICATES = 10000;
        ordered_sig_vector = np.arange(N_SAMPLES).reshape(N_SAMPLES,1);
        random_sig_values = np.repeat(ordered_sig_vector, NUM_REPLICATES, axis=1)
        for j in xrange(random_sig_values.shape[1]):
            np.random.shuffle(random_sig_values[:,j]);

        random_predictions = np.dot(weights, random_sig_values);
        random_scores = np.median(np.abs(ordered_sig_vector - random_predictions), axis=0);

        mu = np.mean(random_scores);
        sigma = np.std(random_scores);

        from scipy.stats import norm;
        p_values = norm.cdf((med_dissimilarity - mu)/sigma);


        sig_proj_matrix[:,i] = med_dissimilarity;
        sig_proj_matrix_p[:,i] = p_values;
        pp.update();

    sig_proj_matrix_p = p_to_q(sig_proj_matrix_p);
    sig_proj_matrix_p[sig_proj_matrix_p == 0] = 1e-300; #Correct for -inf
    sig_proj_matrix_p = np.log10(sig_proj_matrix_p);

    pp.complete();

    return (sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p);

def p_to_q(p_values):
    """
    Uses the Benjamini-Hochberg procedure to convert p_values to q_values

    :param p_values:  numpy.ndarray of p_values
    :return q_values: numpy.ndarray of q_values, same shape as p_values
    """
    original_shape = p_values.shape;
    p_vals_flat = p_values.flatten();
    rank = p_vals_flat.argsort().argsort()+1;
    num_tests = p_values.size;
    q_vals = p_vals_flat * num_tests / rank;
    q_vals.shape = original_shape;

    return q_vals;

def load_precomputed(filename, sample_labels):
    """
    Reads precomputed signature values form a tab-delimited text file
    First row of the file contains sample labels that the signatures correspond with
    Each subsequent row contains a signature name in the first column, followed by the signature values,
        one for each sample label in the file

    :param filename: signature score file name
    :param sample_labels: labels for which we want the signature scores
    :return: a dictionary mapping signature names (string) to signature scores (np.ndarray)
        Each signature score consists of an array with signatures corresponding, by position,
        with the sample labels in the sample_labels argument.
    """
    with open(filename, 'r') as fin:
        line1 = fin.readline().strip().split('\t');
        if(line1[0] == ""): #Remove empty upper-left cell if present
            line1 = line1[1:];

        #match indices between signatures in file and sample_labels
        #want x such that file_cols[x] == sample labels
        target_l = [sl.lower() for sl in sample_labels];
        source_l = [sl.lower() for sl in line1];
        translation_indices = np.zeros(len(target_l), dtype=np.int32);
        for i in xrange(translation_indices.size):
            try:
                translation_indices[i] = source_l.index(target_l[i]);
            except ValueError:
                raise ValueError("Error: Missing value in precomputed signatures for sample " + target_l[i]);

        #Gather signatures
        sig_scores = dict();
        for line in fin:
            line = line.strip();
            if(line == ""): continue;
            s_line = line.split("\t");
            sig_name = s_line[0];
            sig_vals = np.array([float(x) for x in s_line[1:]]);
            sig_vals = sig_vals[translation_indices];
            sig_scores[sig_name] = sig_vals;

        return sig_scores;

class Signature:
    
    def __init__(self, sig_dict, signed, filename, name):
        self.sig_dict = sig_dict;
        self.signed = signed;
        self.source = filename;
        self.name = name;

    def sig_indices(self, genes):
        """Helper method
        
        Returns an array with length = len(genes)
        Entries in the array are 0 if the gene is not in the signature
        Otherwise, value is determined by the signature type"""
        
        out = np.zeros((len(genes), 1), dtype=np.float64);
        
        for i, gene in enumerate(genes):
            if(self.sig_dict.has_key(gene)):
                val = self.sig_dict[gene];
                if(val == 1 or val == 0):
                    out[i] = 1;
                if(val == -1):
                    out[i] = -1;

        #If signed, weight the indices such that the sum of positive signatures
        #counts as much as the sum of negative signatures
        #Weights result in mean(data[pos_sig])/2 - mean(data[neg_sig])
        #                   =  mean(data*sig)
        #Results in mean(out) = 0
        #           mean(|out|) = 1 
        if(self.signed):
            num_pos = np.count_nonzero(out ==  1);
            num_neg = np.count_nonzero(out == -1);
            if(num_pos > 0 and num_neg > 0):
                num_total = num_pos + num_neg;
                
                pos_weight = num_total/num_pos/2;
                neg_weight = num_total/num_neg/2;
    
                out[out==1]  = pos_weight;
                out[out==-1] = neg_weight*-1;
         
        return out;
    