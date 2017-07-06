# -*- coding: utf-8 -*-
"""Functions for loading and manipulating signatures
"""
from __future__ import absolute_import, print_function, division;

import random
import numpy as np;
from sklearn.metrics.pairwise import pairwise_distances;
from scipy.spatial.distance import cdist;
from scipy.stats import norm
from .Utils import ProgressBar;
from .Global import RANDOM_SEED;
from . import SigScoreMethods;
from .SigScoreMethods import SignatureScores

#This is used to cache the background distribution used when evaluating
#Signatures vs projections.  No need to regenerate the random indices
#when the size has not changed.  Saves significant time for large N_SAMPLES
_bg_dist = np.zeros((0,0));
def get_bg_dist(N_SAMPLES, NUM_REPLICATES):
    global _bg_dist;
    if(_bg_dist.shape[0] != N_SAMPLES or _bg_dist.shape[1] != NUM_REPLICATES):
        np.random.seed(RANDOM_SEED);
        _bg_dist = np.random.rand(N_SAMPLES, NUM_REPLICATES);
        _bg_dist = np.argsort(_bg_dist, axis=0);

    return _bg_dist;

def read_signatures(filename='', match_terms=[]):
    """Calls either read_signatures_txt or read_signatures_gmt when appropriate
    """

    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();

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
    
    ff = open(filename, 'rU');
    found_signatures = dict();    ## dictionary of signatures
    
    if(type(match_terms) is str):
        match_terms = [match_terms];
    
    match_terms = [term.lower() for term in match_terms];
    
    try:    

        for i, line in enumerate(ff):
            sline = line.strip();
            if(sline[0] == '#' or sline[0:2] == '//'):
                continue;

            row_data = sline.split('\t');
            if(len(row_data) == 2):
                signed = False;
            elif(len(row_data) == 3):
                signed = True;
            else:
                raise ValueError("Line " + str(i) + " Signature file should contain 2 (unsigned) or 3 (signed) columns");

            name = row_data[0];
            if(name not in found_signatures):  ## Add the signature if we haven't seen it yet

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
                gene_name = row_data[2].upper();
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

            else:  #unsigned case
                sig_val = 0;
                gene_name = row_data[1].upper();

            sig_dict[gene_name] = sig_val;
    except:
        raise;
    finally:
        ff.close();
    
    return [found_signatures[key] for key in found_signatures]  #dict to list

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

    plus_terms = ['plus', 'up'];
    minus_terms = ['minus', 'dn', 'down'];

    def sig_parts(name):
        """
        Returns: The name of the signature with the plus or minus suffix removed.
                and whether or not the signature is positive or negative
                    : as a string, either 'plus' or 'minus'
        """
        lname = name.lower();
        sign = 'unsigned';
        ii = len(lname);
        for term in plus_terms:
            if(lname.endswith('_' + term)):
                ii = lname.rfind('_' + term);
                sign = 'plus';

        if(sign != 'plus'):
            for term in minus_terms:
                if(lname.endswith('_' + term)):
                    ii = lname.rfind('_' + term);
                    sign = 'minus';

        return name[0:ii], sign;

    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();

    if(type(match_terms) is str):
        match_terms = [match_terms];

    match_terms = [term.lower() for term in match_terms];

    with open(filename, 'rU') as ff:
        found_signatures = dict();    ## dictionary of signatures
        end_of_file = False;

        while(not end_of_file):
            line = ff.readline().strip();
            if(line == ""):
                end_of_file = True;
                continue;
            if(line[0] == '#' or line[0:2] == '//'):
                continue;

            row_data = line.split('\t');
            name = row_data[0];
            root_name, sign = sig_parts(name);

            ## Only add signatures if they match one of the supplied terms
            matched = False if len(match_terms) > 0 else True
            lname = root_name.lower();
            for term in match_terms:
                if(term in lname):
                    matched = True;
                    break;
            if(not matched): continue;

            isSigned = sign == 'plus' or sign == 'minus';

            #Check if sig exists, create/insert if not
            if(root_name in found_signatures):
                sig = found_signatures[root_name];
            else:
                sig = Signature(dict(), isSigned, filename, root_name);
                found_signatures.update({root_name: sig});

            if(sign == 'plus'):
                sig_val = 1;
            elif(sign == 'minus'):
                sig_val = -1;
            elif(sign == 'unsigned'):
                sig_val = 0;
            else:
                raise Exception("This should not happen");

            for gene_name in row_data[2:]:
                gene_name_split = gene_name.split(',')
                if(len(gene_name_split) == 1):
                    sig.sig_dict.update({gene_name.upper(): sig_val});
                elif(len(gene_name_split) == 2):
                    sig.sig_dict.update({gene_name_split[0].upper(): float(gene_name_split[1])});
                else:
                    raise Exception("Error in signature: " + name)


    return list(found_signatures.values());  #dict to list

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
    
    match_terms = [term.lower() for term in match_terms];    
    
    filtered_signatures = list();
    for sig in signatures:
        name = sig.name.lower();
        
        for term in match_terms:            
            if(name.find(term) >= 0):
                filtered_signatures.append(sig);
                break;
    
    
    return filtered_signatures;

def sigs_vs_projections(projections, sig_scores_dict, random_sig_scores_dict, NEIGHBORHOOD_SIZE = 0.33):
    """
    Evaluates the significance of each signature vs each projection

    :param projections: dict of (string) => (numpy.ndarray of shape 2xNum_Samples)
        Maps projections to their spatial coordinates for each sample
    :param sig_scores_dict: dict of (string) => (numpy.ndarray of shape Num_Samples)
        Maps signature names to their value at each coordinate
    :return: tuple

        sp_row_labels: List of Strings. Labels for rows of the output matrices
        sp_col_labels: List of Strings.  Labels for columns of the output matrices
        sig_proj_matrix: numpy.ndarray, NUM_SIGNATURES x NUM_PROJECTIONS
                        sig_proj dissimilarity score
        sig_proj_matrix_p: numpy.ndarray, NUM_SIGNATURES x NUM_PROJECTIONS
                        sig_proj dissimilarity p value
    """
    np.random.seed(RANDOM_SEED);
    sp_row_labels = [];
    sp_row_labels_factors = [];
    sp_row_labels_pnum = [];

    #Remove signatures that are factor signatures and precomputed-numerical signatures
    for name, sig_scores in sig_scores_dict.items():
        if(sig_scores.isPrecomputed):
            if(sig_scores.isFactor):
                sp_row_labels_factors.append(name);
            else:
                sp_row_labels_pnum.append(name);
        else:
            sp_row_labels.append(name);



    sp_col_labels = list(projections.keys());
    sp_col_labels.sort();

    N_SAMPLES = len(sig_scores_dict[
        list(sig_scores_dict.keys())[0]
    ].sample_labels)
    N_SIGNATURES = len(sp_row_labels);
    N_SIGNATURES_FACTORS = len(sp_row_labels_factors);
    N_SIGNATURES_PNUM = len(sp_row_labels_pnum);
    N_PROJECTIONS = len(sp_col_labels);


    sig_proj_matrix   = np.zeros((N_SIGNATURES,N_PROJECTIONS));
    sig_proj_matrix_p = np.zeros((N_SIGNATURES,N_PROJECTIONS));

    factor_sig_proj_matrix   = np.zeros((N_SIGNATURES_FACTORS,N_PROJECTIONS));
    factor_sig_proj_matrix_p = np.zeros((N_SIGNATURES_FACTORS,N_PROJECTIONS));

    pnum_sig_proj_matrix = np.zeros((N_SIGNATURES_PNUM, N_PROJECTIONS));
    pnum_sig_proj_matrix_p = np.zeros((N_SIGNATURES_PNUM, N_PROJECTIONS));

    #Build a matrix of all signatures
    sig_score_matrix = np.zeros((N_SAMPLES, N_SIGNATURES));

    for j, sig in enumerate(sp_row_labels):
        sig_score_matrix[:,j] = sig_scores_dict[sig].ranks;

    # Build a matrix of random signatures
    random_sig_score_matrix = np.zeros((N_SAMPLES, len(random_sig_scores_dict)));
    random_sig_score_keys = list(random_sig_scores_dict.keys());

    for j, sig in enumerate(random_sig_score_keys):
        random_sig_score_matrix[:, j] = random_sig_scores_dict[sig].ranks;

    #Build one-hot matrices for each factor
    factor_dict = dict();
    for sig in sp_row_labels_factors:
        factor_values = sig_scores_dict[sig].scores;
        factor_levels = list(set(factor_values)); #Makes unique
        factor_frequencies = np.zeros(len(factor_levels));
        factor_matrix = np.zeros((N_SAMPLES, 0));
        for j, fval in enumerate(factor_levels):
            factor_matrix_row = np.zeros((N_SAMPLES, 1));
            equal_ii = [i for i,val in enumerate(factor_values) if val == fval];
            factor_matrix_row[equal_ii] = 1;
            factor_frequencies[j] = len(equal_ii) / len(factor_values);
            factor_matrix = np.concatenate((factor_matrix, factor_matrix_row), axis=1);
        factor_dict[sig] = (factor_levels, factor_frequencies, factor_matrix);


    print();
    print("Evaluating Signatures against Projections");
    pp = ProgressBar(N_PROJECTIONS);
    for i, proj in enumerate(sp_col_labels):
        data_loc = projections[proj];

        distance_matrix = cdist(data_loc.T, data_loc.T, metric='euclidean');

        weights = np.exp(-1 * distance_matrix**2 / NEIGHBORHOOD_SIZE**2);
        weights[np.arange(N_SAMPLES), np.arange(N_SAMPLES)] = 0; #Don't count self
        weights_norm_factor = np.sum(weights, axis=1, keepdims=True);
        weights_norm_factor[weights_norm_factor == 0] = 1.0;
        weights /= weights_norm_factor;

        neighborhood_prediction = np.dot(weights, sig_score_matrix);


        ##Neighborhood dissimilarity score = |actual - predicted|
        dissimilarity = np.abs(sig_score_matrix - neighborhood_prediction);
        med_dissimilarity = np.median(dissimilarity, axis=0);

        # Calculate scores for random signatures
        random_neighborhood_prediction = np.dot(weights, random_sig_score_matrix);
        random_dissimilarity = np.abs(random_sig_score_matrix - random_neighborhood_prediction);
        random_med_dissimilarity = np.median(random_dissimilarity, axis=0);

        # Group by number of genes
        backgrounds = dict();
        for k, rsig in enumerate(random_sig_score_keys):
            numGenes = random_sig_scores_dict[rsig].numGenes;
            if(numGenes not in backgrounds):
                backgrounds.update({numGenes: []});
            backgrounds[numGenes].append(random_med_dissimilarity[k]);

        bg_stat = np.zeros((len(backgrounds), 3));
        for k, numGenes in enumerate(backgrounds.keys()):
            mu_x = np.mean(backgrounds[numGenes]);
            std_x = np.std(backgrounds[numGenes]);
            bg_stat[k, 0] = numGenes;
            bg_stat[k, 1] = mu_x;
            bg_stat[k, 2] = std_x;

        mu = np.zeros(med_dissimilarity.shape);
        sigma = np.zeros(med_dissimilarity.shape);
        for k in range(med_dissimilarity.size):
            # find background with closest number of genes
            numGenes = sig_scores_dict[sp_row_labels[k]].numGenes;
            row_i = np.argmin(np.abs(numGenes - bg_stat[:, 0]));
            mu[k] = bg_stat[row_i, 1];
            sigma[k] = bg_stat[row_i, 2];


        p_values = norm.cdf((med_dissimilarity - mu)/sigma);

        sig_proj_matrix[:,i] = 1 - med_dissimilarity / N_SAMPLES;
        sig_proj_matrix_p[:,i] = p_values;


        # Calculate significance for precomputed numerical signatures
        # This is done separately because there are likely to be many repeats (e.g. for a time coordinate)
        for j,sig in enumerate(sp_row_labels_pnum):

            sig_scores = sig_scores_dict[sig].ranks;

            if((sig_scores == sig_scores[0]).all()): # If they are all the same, p-value is 1.0
                pnum_sig_proj_matrix[j, i] = 0.0;
                pnum_sig_proj_matrix_p[j, i] = 1.0;
                continue;

            sig_scores = sig_scores.reshape((sig_scores.size,1));
            sig_predictions = np.dot(weights, sig_scores);
            dissimilarity = np.abs(sig_scores - sig_predictions);
            med_dissimilarity = np.median(dissimilarity);

            # Now compute a background
            NUM_REPLICATES = 10000;
            random_sig_values = get_bg_dist(N_SAMPLES, NUM_REPLICATES);
            bg_values = sig_scores.ravel()[random_sig_values];
            random_predictions = np.dot(weights, bg_values);
            random_scores = np.median(np.abs(bg_values - random_predictions), axis=0);

            mu = np.mean(random_scores);
            sigma = np.std(random_scores);
            if(sigma != 0):
                p_value = norm.cdf((med_dissimilarity - mu)/sigma);
            else:
                p_value = 1.0;


            pnum_sig_proj_matrix[j,i] = 1 - med_dissimilarity / N_SAMPLES;
            pnum_sig_proj_matrix_p[j,i] = p_value;

        #Calculate significance for Factor signatures
        for j,sig in enumerate(sp_row_labels_factors):
            factor_levels, factor_frequencies, factor_matrix = factor_dict[sig];

            if((factor_frequencies == 1).any()):  # Act accordingly if all values are the same
                factor_sig_proj_matrix[j, i] = 0.0;
                factor_sig_proj_matrix_p[j, i] = 1.0;
                continue;

            N_LEVELS = len(factor_levels);
            factor_predictions = np.dot(weights, factor_matrix);

            dissimilarity = 1 - np.sum(factor_matrix * factor_predictions, axis=1);
            med_dissimilarity = np.median(dissimilarity);

            #Now...compute a background
            NUM_REPLICATES = 1000;
            column_assignments = np.random.choice(N_LEVELS, NUM_REPLICATES, p = factor_frequencies);
            column_assignments = factor_frequencies[column_assignments];
            column_assignments = column_assignments.reshape((1,NUM_REPLICATES));
            rand_factors = np.random.rand(N_SAMPLES, NUM_REPLICATES);
            rand_factors = (rand_factors < column_assignments).astype('int');
            random_predictions = np.dot(weights, rand_factors);
            for ii in range(random_predictions.shape[0]):
                random_predictions[ii] = np.random.permutation(random_predictions[ii]);
            rand_med_dissimilarity = np.median(1-random_predictions, axis=0);

            mu = np.mean(rand_med_dissimilarity);
            sigma = np.std(rand_med_dissimilarity);
            if(sigma == 0):
                p_value = 1;
            else:
                p_value = norm.cdf((med_dissimilarity - mu)/sigma);

            factor_sig_proj_matrix[j,i] = 1 - med_dissimilarity;
            factor_sig_proj_matrix_p[j,i] = p_value;


        pp.update();

    #Concatenate the Factor sig-proj entires back in
    sig_proj_matrix = np.concatenate((sig_proj_matrix, factor_sig_proj_matrix, pnum_sig_proj_matrix), axis=0);
    sig_proj_matrix_p = np.concatenate((sig_proj_matrix_p, factor_sig_proj_matrix_p, pnum_sig_proj_matrix_p), axis=0);
    sp_row_labels = sp_row_labels + sp_row_labels_factors + sp_row_labels_pnum;

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
    q_vals[q_vals > 1.0] = 1.0;

    return q_vals;

def load_precomputed(filename, sample_labels):
    """
    Reads precomputed signature values form a tab-delimited text file
    First row of the file contains sample labels that the signatures correspond with
    Each subsequent row contains a signature name in the first column,
         followed by the signature type (either 'numerical' or 'factor')
         followed by the signature values, one for each sample label in the file

    :param filename: signature score file name
    :param sample_labels: labels for which we want the signature scores
    :return: Dictionary of Signature Name (String) -> Signature (SignatureScores)
    """

    with open(filename, 'rU') as fin:
        #Determine column labels to apply
        line1 = fin.readline().rstrip().split('\t');
        line2 = fin.readline().rstrip().split('\t');

        if(len(line1) == len(line2)): #First two entries must be empty or column headers
            line1 = line1[2:];
        elif(len(line1) == len(line2)-2):
            line1 = line1;
        else: #Other arrangements signify some unusual formatting
            raise ValueError("Error in header line of precomputed signature file.\n"
                 + "First row should contain tab-separated list of samples");


        #match indices between signatures in file and sample_labels
        #want x such that file_cols[x] == sample labels
        target_l = [sl.lower() for sl in sample_labels];
        source_l = [sl.lower() for sl in line1];
        translation_indices = np.zeros(len(target_l), dtype=np.int32);
        for i in range(translation_indices.size):
            try:
                translation_indices[i] = source_l.index(target_l[i]);
            except ValueError:
                raise ValueError("Error: Missing value in precomputed signatures for sample " + target_l[i]);

        fin.seek(0);
        xx = fin.readline();

        #Gather signatures
        sig_scores = dict();
        for line in fin:
            line = line.strip();
            if(line == ""): continue;
            s_line = line.split("\t");
            sig_name = s_line[0];
            sig_type = s_line[1].strip().lower();
            sig_val_cells = s_line[2:];

            if(sig_type == 'numerical'):
                sig_isFactor = False;
                try:
                    sig_vals = np.array([float(x) for x in sig_val_cells]);
                    sig_vals = sig_vals[translation_indices];
                except ValueError as e:
                    print(e.message);
                    print('Error in precomputed signature:', sig_name);
                    for i,x in enumerate(sig_val_cells):
                        try:
                            y = float(x);
                        except ValueError:
                            print("Error in column", i);
                            print("Bad value:", x);
                    raise Exception('Failed to load precomputed signature. Correct file format and re-run.');
            elif(sig_type == 'factor'):
                sig_isFactor = True;
                sig_vals = [sig_val_cells[i] for i in translation_indices];
            else:
                raise ValueError('Column 2 of precomputed signature file should specify either "numerical" or "factor"');

            sig_scores[sig_name] = SignatureScores(sig_vals, sig_name, sample_labels, sig_isFactor, isPrecomputed=True, numGenes=0);

        return sig_scores;


def calculate_sig_scores(data, signatures, method='weighted_avg',
                         zero_locations=None, min_signature_genes=1e99):
    """
    Generate random signatures to be used as a background distribution

    Parameters
    ----------
    data: ExpressionData
        Data matrix to use to calculate signature scores
    signatures: list of Signature
        Signatures of which scores are computed
    method: str
        Which scoring method to use
    zero_locations: boolean numpy.ndarray
        Same size as data
        True where 'data' was originally zero
        Used for normalization methods that transform zeros to some other value
    min_signature_genes:  numeric
        Signatures that match less that `min_signature_genes` are discarded

    Returns
    -------
    dict of str -> SignatureScores
        key is name of Signature
        value is SignatureScores object
    """

    sig_scores_dict = dict()

    # Determine signature score evaluation method
    if method == "naive":
        sig_score_method = SigScoreMethods.naive_eval_signature
    elif method == "weighted_avg":
        sig_score_method = SigScoreMethods.weighted_eval_signature
    elif method == "imputed":
        sig_score_method = SigScoreMethods.imputed_eval_signature
    elif method == "only_nonzero":
        sig_score_method = SigScoreMethods.nonzero_eval_signature

    pbar = ProgressBar(len(signatures))
    for sig in signatures:

        try:
            sig_scores_dict[sig.name] = sig_score_method(
                data, sig, zero_locations, min_signature_genes)
        except ValueError:  # Thrown when the signature has too few genes in the data
            pass  # Just discard the signature then

        pbar.update()
    pbar.complete()

    return sig_scores_dict


def generate_random_sigs(features, signed, sizes=None, num_per_size=3000):
    """
    Generate random signatures to be used as a background distribution

    Parameters
    ----------
    features: list of str
        List of features (e.g. gene symbols) to use in signatures
    signed: bool
        Whether or not the signature should be signed
    sizes: list of int
        How many features to draw for background signatures
    num_per_size: int
        How many signatures to create at each size in `sizes`

    Returns
    -------
    list of Signature
        Randomly generate background signatures
    """

    if sizes is None:
        sizes = [5, 10, 20, 50, 100, 200]

    random_sigs = []

    if signed:
        signs = [-1, 1]
    else:
        signs = [1]

    for size in sizes:
        for j in range(num_per_size):
            new_sig_dict = dict()
            new_sig_genes = random.sample(features, size)
            new_sig_signs = np.random.choice(signs, size)

            for gene, sign in zip(new_sig_genes, new_sig_signs):
                new_sig_dict.update({gene: int(sign)})

            new_sig = Signature(new_sig_dict, signed, 'x',
                                name="RANDOM_BG_" + str(size) + "_" + str(j))
            random_sigs.append(new_sig)

    return random_sigs

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
            if(gene in self.sig_dict):
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
        if(False):
            num_pos = np.count_nonzero(out ==  1);
            num_neg = np.count_nonzero(out == -1);
            if(num_pos > 0 and num_neg > 0):
                num_total = num_pos + num_neg;
                
                pos_weight = num_total/num_pos/2;
                neg_weight = num_total/num_neg/2;
    
                out[out==1]  = pos_weight;
                out[out==-1] = neg_weight*-1;
         
        return out;
