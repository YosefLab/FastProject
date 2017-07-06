# -*- coding: utf-8 -*-
""" Main FastProject Pipeline

Here is the main execution path for FastProject
While this defines the types of processining and the
order in which it is done, the goal should be to organize
as much of that off into other modules.
"""
from __future__ import division, print_function;
import numpy as np;
import pandas as pd;
import time;
from . import Filters;
from . import Transforms;
from . import Signatures;
from . import SigScoreMethods;
from . import Projections;
from . import SubSample;
from .DataTypes import ExpressionData;
from .Utils import ProgressBar;
from .Global import FP_Output;
from . import NormalizationMethods;


def Analysis(expressionMatrix, signatures, precomputed_signatures, housekeeping_genes, input_projections, input_weights, kwargs):
    """Run the full FastProject analysis pipeline

    Parameters
    ----------
    expressionMatrix : ExpressionData
    signatures : list of Signatures.Signature
    precomputed_signatures : dict
        Keys are precomputed signature names (str)
        Values are signature levels/scores (SigScoreMethods.SignatureScores)
    housekeeping_genes : list of str
    input_projections : dict
        Keys are of type str, representing projection names
        Values are of type 2xN pandas.DataFrame with column names matching
        sample names in `expressionMatrix`
    input_weights : pandas.DataFrame or None
        Same size as expressionMatrix
        Values are floats from 0.0 to 1.0
    **kwargs : dict
        Additional options

    Returns
    -------
    Models : dict
        An object containing all analysis results
        modelName(str) -> modelData(dict)
            "Data" -> ExpressionData
            "signatureScores" -> list of SigScoreMethods.SignatureScores
            "sampleLabels" -> list of str, labels for each sample
            "projectionData" -> list of dict
                "filter" -> str, name of filter used
                "genes" -> list of str, list of genes that passed filter
                "pca" -> bool, whether or not PCA performed before project
                "projections" -> dict
                    projection name (str) -> 2xN numpy.ndarray projection coords
                "sigProjMatrix" -> np.ndarray, (sig, proj) consistency scores
                "sigProjMatrix_p" -> np.ndarray, (sig, proj) p-values
                "projectionKeys" -> list of str, column labels for sigProjMatrix
                "signatureKeys" -> list of str, row labels for sigProjMatrix
                "clusters" -> dict
                    projection name (str) -> projection cluster assignments (dict)
                        cluster method (str) -> cluster assignments (1D np.ndarray)
                "loadings" -> np.ndarray
                    np.ndarray, shape [num_genes x 3], first 3 PCA loadings
                    only present if 'pca' is true

    qc_info : pandas.DataFrame
        Columns are "Scores" and "Passes"
        "Score": quality score
        "Passes": bool, whether or not the sample passed
        Index: sample label

    """

    start_time = time.time();

    #Wrap data in ExpressionData object, add as a Model
    all_data = expressionMatrix;
    edata = expressionMatrix;

    # if(kwargs["subsample_size"] > edata.shape[1]):
    #     kwargs["subsample_size"] = None;
    kwargs["subsample_size"] = None;

    #holdouts = None
    #if(kwargs["subsample_size"])
    #    holdouts, edata = SubSample.split_samples(edata, kwargs["subsample_size"])

    if(kwargs["threshold"] is None):
        # Default threshold is 20% of samples - post sub-sampling
        kwargs["threshold"] = int(0.2 * edata.shape[1]);

    #Hold on to originals so we don't lose data after filtering in case it's needed later
    original_data = pd.DataFrame(edata.base,
                                 index=edata.row_labels,
                                 columns=edata.col_labels)
    original_data.index = original_data.index.str.upper()

    #Filtering
    edata = Filters.apply_filters(edata, kwargs["threshold"], kwargs["nofilter"], kwargs["lean"]);

    Models = dict();
    emodel = dict({"Data": edata});
    Models.update({"Expression": emodel});

    #%% Probability transform
    if(not kwargs["nomodel"]):

        FP_Output('\nEstimating non-detect events');
        p_nd = Transforms.estimate_non_detection(original_data)

        FP_Output('\nCorrecting for false-negatives using housekeeping gene levels');
        (fit_func, params) = Transforms.create_false_neg_map(original_data, p_nd, housekeeping_genes);

        if(input_weights is None):
            weights = Transforms.compute_weights(fit_func, params, edata, p_nd)
        else:
            weights = input_weights.loc[edata.row_labels, edata.col_labels].values;

        edata.weights = weights;

        sample_passes_qc, sample_qc_scores = Transforms.quality_check(params);
        qc_info = pd.DataFrame({"Score": sample_qc_scores, "Passes": sample_passes_qc}, index=edata.col_labels);

        #If specified, remove items that did not pass qc check
        if(kwargs["qc"]):
            for name, model in Models.items():
                dataMatrix = model["Data"];
                model["Data"] = dataMatrix.subset_samples(sample_passes_qc);

            sample_qc_scores = sample_qc_scores[sample_passes_qc];
            p_nd = p_nd.loc[:, sample_passes_qc]
    else:
        sample_qc_scores = None
        qc_info = pd.DataFrame(0.0, columns=["Score", "Passes"],
                               index=edata.col_labels)
        p_nd = pd.DataFrame(1.0, index=original_data.index,
                            columns=original_data.columns)

    # Necessary because the matrix might be modified when data failing qc is removed
    edata = Models["Expression"]["Data"];

    # Make an extra matrix just storing the location of zeros in the original matrix
    zero_locations = p_nd.loc[edata.row_labels].values

    # Make an extra quality score of just the zero proportion in each sample
    zeros_qscore = zero_locations.mean(axis=0)

    # Transforms.z_normalize(edata);

    # Generate some random signatures for testing purposes
    # for size in [5, 10, 20, 50, 100, 200]:
    #     for j in range(100):
    #         new_sig_dict = dict();
    #         new_sig_genes = random.sample(edata.row_labels, size);
    #         new_sig_signs = np.random.choice([1, -1], size);
    #         for gene, sign in zip(new_sig_genes, new_sig_signs):
    #             new_sig_dict.update({gene: int(sign)});
    #         new_sig = Signatures.Signature(new_sig_dict, True, 'x', "RANDOM_SIGNED_" + str(size) + "_" + str(j));
    #         sigs.append(new_sig);

    # # Generate some positive random signatures too (only + sign)
    # for size in [5, 10, 20, 50, 100, 200]:
    #     for j in range(100):
    #         new_sig_dict = dict();
    #         new_sig_genes = random.sample(edata.row_labels, size);
    #         new_sig_signs = np.random.choice([1], size);
    #         for gene, sign in zip(new_sig_genes, new_sig_signs):
    #             new_sig_dict.update({gene: int(sign)});
    #         new_sig = Signatures.Signature(new_sig_dict, True, 'x', "RANDOM_UNSIGNED_" + str(size) + "_" + str(j));
    #         sigs.append(new_sig);

    # Generate random signatures for background significance
    random_sigs = Signatures.generate_random_sigs(edata.row_labels, signed=False)

    for name, model in Models.items():

        data = model["Data"];

        FP_Output('\nModel: ', name)

        #Evaluate Signatures
        FP_Output("\nEvaluating signature scores on samples...");

        # Determine normalization method
        if(kwargs["sig_norm_method"] == "none"):
            sig_norm_method = NormalizationMethods.no_normalization;
        elif(kwargs["sig_norm_method"] == "znorm_columns"):
            sig_norm_method = NormalizationMethods.col_normalization;
        elif(kwargs["sig_norm_method"] == "znorm_rows"):
            sig_norm_method = NormalizationMethods.row_normalization;
        elif(kwargs["sig_norm_method"] == "znorm_rows_then_columns"):
            sig_norm_method = NormalizationMethods.row_and_col_normalization;
        elif(kwargs["sig_norm_method"] == "rank_norm_columns"):
            sig_norm_method = NormalizationMethods.col_rank_normalization;

        # Normalize data with it
        sig_data = data.get_normalized_copy(sig_norm_method);

        sig_scores_dict = Signatures.calculate_sig_scores(sig_data, signatures,
                              method=kwargs["sig_score_method"],
                              zero_locations=zero_locations,
                              min_signature_genes=kwargs["min_signature_genes"])

        # Add in precomputed signature scores 'meta-data'
        # Need to filter these as the order/number of samples might have changed
        for name, sigscores in precomputed_signatures.items():
            ii = [sigscores.sample_labels.index(x) for x in sig_data.col_labels]
            new_labels = [sigscores.sample_labels[i] for i in ii]
            new_scores = [sigscores.scores[i] for i in ii]
            sigscores.sample_labels = new_labels
            sigscores.scores = new_scores

        sig_scores_dict.update(precomputed_signatures)

        #Adds in quality score as a pre-computed signature
        if(sample_qc_scores is not None): #Might be None if --nomodel option is selected
            sig_scores_dict["FP_Quality"] = SigScoreMethods.SignatureScores(sample_qc_scores,"FP_Quality",data.col_labels,isFactor=False, isPrecomputed=True, numGenes=0);

        #Adds in zero proportion as another pre-computed signature
        sig_scores_dict["Zero_Proportion"] = SigScoreMethods.SignatureScores(zeros_qscore,"Zero_Proportion",data.col_labels,isFactor=False, isPrecomputed=True, numGenes=0);

        FP_Output("\nEvaluating null signature scores on samples...");

        random_sig_scores_dict = Signatures.calculate_sig_scores(sig_data, random_sigs,
                              method=kwargs["sig_score_method"],
                              zero_locations=zero_locations,
                              min_signature_genes=kwargs["min_signature_genes"])

        model["signatureScores"] = sig_scores_dict;
        model["projectionData"] = [];
        model["sampleLabels"] = data.col_labels;

        for filter_name in data.filters:
            projData = dict();

            FP_Output("\nFilter-Level:", filter_name);
            #%% Dimensional Reduction procedures
            FP_Output("\nProjecting data into 2 dimensions");

            projections, pcdata = Projections.generate_projections(data, filter_name,
                                                                   input_projections, kwargs["lean"]);

            #Evaluate Clusters
            FP_Output("Evaluating Clusters...");
            clusters = Projections.define_clusters(projections);

            #%% Evaluating signatures against projections
            # sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections(projections, sig_scores_dict);
            sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections(projections, sig_scores_dict, random_sig_scores_dict);

            #Store in projData
            projData["filter"] = filter_name;
            projData["genes"] = data.filtered_genes(filter_name);
            projData["pca"] = False;
            projData["projections"] = projections;
            projData["sigProjMatrix"] = sig_proj_matrix;
            projData["sigProjMatrix_p"] = sig_proj_matrix_p;
            projData["projectionKeys"] = sp_col_labels;
            projData["signatureKeys"] = sp_row_labels;
            projData["clusters"] = clusters;
            model["projectionData"].append(projData);


            #Now do it all again using the principal component data
            # if(kwargs["pca_filter"]):
            #     pcdata = Projections.filter_PCA(pcdata, scores=sample_qc_scores, variance_proportion=0.25);
            # else:
            #     pcdata = Projections.filter_PCA(pcdata, variance_proportion=0.25, min_components = 30);

            #%% Dimensional Reduction procedures
            FP_Output("Projecting PC data into 2 dimensions");

            projections, pcdata2 = Projections.generate_projections(pcdata, filter_name, lean=kwargs["lean"]);

            #Evaluate Clusters
            FP_Output("Evaluating Clusters...");
            clusters = Projections.define_clusters(projections);

            #%% Evaluating signatures against projections
            #sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections(projections, sig_scores_dict);
            sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections(projections, sig_scores_dict, random_sig_scores_dict);

            projData = dict();
            projData["filter"] = filter_name;
            projData["genes"] = data.filtered_genes(filter_name);
            projData["pca"] = True;
            projData["projections"] = projections;
            projData["sigProjMatrix"] = sig_proj_matrix;
            projData["sigProjMatrix_p"] = sig_proj_matrix_p;
            projData["projectionKeys"] = sp_col_labels;
            projData["signatureKeys"] = sp_row_labels;
            projData["clusters"] = clusters;
            projData["loadings"] = pcdata.loadings[:,0:3];
            model["projectionData"].append(projData);

    #Filter output signatures
    for name, model in Models.items():
        signatureScores = model["signatureScores"];

        signature_significance = pd.DataFrame();
        for projData in model['projectionData']:
            sp = pd.DataFrame(projData['sigProjMatrix_p'], index=projData['signatureKeys']).min(axis=1);
            signature_significance = pd.concat((signature_significance, sp), axis=1);

        signature_significance = signature_significance.min(axis=1);

        #Determine a threshold of significance
        #If too many samples, hard limit the number of output signatures to conserve file size
        if(kwargs["all_sigs"]): # Keep all sigs if this flag is set
            threshold = 1e99; 
        else:
            OUTPUT_SIGNATURE_LIMIT = min(200, signature_significance.size);
            if(all_data.shape[1] > 2000):
                aa = np.argsort(signature_significance);
                threshold = signature_significance[aa[OUTPUT_SIGNATURE_LIMIT - 1]];
            else:
                threshold = -1.3;
                aa = np.argsort(signature_significance);
                # Ensure at least N = Output_Signature_Limit signatures
                # If the Nth best signature is over the threshold, raise the threshold
                if(signature_significance[aa[OUTPUT_SIGNATURE_LIMIT-1]] > threshold):
                    threshold = signature_significance[aa[OUTPUT_SIGNATURE_LIMIT-1]];

        #Iterate back through and prune signatures worse than threshold
        #Create a dictionary of sigs to keep
        keep_sig = dict();
        for name, sig_score in signatureScores.items():
            if(not sig_score.isPrecomputed and signature_significance[name] > threshold):
                keep_sig.update({name: False});
            else:
                keep_sig.update({name: True});

        #Remove values in the model's signatureScores dict
        new_signatureScores = {}
        for name, sig_score in signatureScores.items():
            if(keep_sig[name] == True):
                new_signatureScores[name] = sig_score

        signatureScores = new_signatureScores

        #Remove values in each filters sigProjMatrix and the sigProjMatrix keys
        for projData in model['projectionData']:
            projData_keep_sig = np.array([keep_sig[sig] for sig in projData["signatureKeys"]]);
            projData["sigProjMatrix"] = projData["sigProjMatrix"][projData_keep_sig,:];
            projData["sigProjMatrix_p"] = projData["sigProjMatrix_p"][projData_keep_sig,:];
            projData["signatureKeys"] = [x for x in projData["signatureKeys"] if keep_sig[x]];


    #Merge all the holdouts back into the model
    #if(kwargs["subsample_size"] is not None):
    #    FP_Output("Merging held-out samples back in")
    #    SubSample.merge_samples(all_data, Models, signatures, prob_params, kwargs);

    # Cluster genes
    from scipy.cluster.hierarchy import leaves_list, linkage;
    edata = Models["Expression"]["Data"];
    edata_norm = edata.get_normalized_copy(NormalizationMethods.row_normalization)
    linkage_matrix = linkage(edata_norm);
    leaves_i = leaves_list(linkage_matrix);
    edata_norm = edata_norm.subset_genes(leaves_i);
    Models["Expression"]["Data"] = edata_norm;

    FP_Output("\nFastProject Analysis Complete")
    elapsed_time = time.time() - start_time;
    FP_Output("Elapsed Time {:.2f} seconds".format(elapsed_time));

    return Models, qc_info;
