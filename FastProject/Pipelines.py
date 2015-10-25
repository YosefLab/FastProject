from __future__ import division, print_function;
import os;
import numpy as np;
import pandas as pd;
import time;
import scipy.stats;
import logging;
from FastProject import Filters;
from FastProject import FileIO;
from FastProject import Transforms;
from FastProject import Signatures;
from FastProject import Projections;
from FastProject import SubSample;
from FastProject.DataTypes import ExpressionData, ProbabilityData, PCData;
from FastProject.Utils import ProgressBar;
from FastProject import HtmlViewer;
from .Global import options, args, FP_Output;


def FullOutput():

    #Create directory for all outputs
    if(options.output):
        dir_name = options.output;
    else:
        default_dir_name = 'FastProject_Output';
        if(os.path.isdir(default_dir_name)):
            i = 1;
            while(True):
                dir_name = default_dir_name + str(i);
                if(not os.path.isdir(dir_name)):
                    break;
                else:
                    i = i+1;
        else:
            dir_name = default_dir_name;

    FileIO.make_dirs(dir_name);
    logging.basicConfig(format='%(asctime)s %(message)s', filename=os.path.join(dir_name, 'fastproject.log'), level=logging.INFO);
    logging.info("Running FastProject Analysis");
    logging.info("Data File: " + args[0]);
    for key in options.__dict__:
        logging.info(key + ": " + str(options.__dict__[key]));

    start_time = time.time();
    housekeeping_filename = options.housekeeping;


    #%% Read expression data from file
    if(len(args) > 0):
        filename = args[0];
    else:
        raise ValueError("Argument Error:  data_file not specified.\nExiting...");

    if(not os.path.isfile(filename)):
        raise ValueError("Argument Error: data file not found.\nExiting...");

    (edata, genes, cells) = FileIO.read_matrix(filename);

    FP_Output("Imported ", edata.shape[0], " genes across ", edata.shape[1], " samples");

    #%% Load Signature file
    sigs = [];
    if(options.signatures):
        sig_file = options.signatures;
        if(not os.path.isfile(sig_file)):
            raise ValueError("Option Error: signature file " + sig_file + " not found.\nExiting...");

        sigs = Signatures.read_signatures(sig_file);
    if(options.precomputed):
        precomputed_sig_file = options.precomputed;
        if(not os.path.isfile(precomputed_sig_file)):
            raise ValueError("Option Error: precomputed signature file " + precomputed_sig_file + " not found.\nExiting...");
        _throwaway = Signatures.load_precomputed(options.precomputed, cells);

    if(not options.signatures and not options.precomputed): #Need one or the other here
        raise ValueError("Option Error: Must specify either a signature file or a pre-computed signature file.\nExiting...");


    #Wrap data in ExpressionData object, add as a Model
    all_data = ExpressionData(edata, genes, cells); #TODO: Does this need to be copied?
    edata = ExpressionData(edata, genes, cells);

    if(options.subsample_size > edata.shape[1]):
        options.subsample_size = None;

    holdouts = None;
    if(options.subsample_size):
        holdouts, edata = SubSample.split_samples(edata, options.subsample_size);


    #Hold on to originals so we don't lose data after filtering in case it's needed later
    original_data = edata.copy();

    #Filtering
    filter_dict = {};
    if(options.nofilter):
        edata = Filters.filter_genes_novar(edata);

        filter_dict.update({'No_Filter': set(edata.row_labels)});
    else:
        edata = Filters.filter_genes_threshold(edata, 0.2);

        #HDT Filtering
        FP_Output("Removing genes with unimodal distribution across samples using Hartigans DT...");
        hdt_mask = Filters.filter_genes_hdt(edata, 0.05);
        #Fano Filtering
        FP_Output("Applying Fano-Filtering...");
        fano_mask = Filters.filter_genes_fano(edata, 2);

        filter_dict.update({
            'Threshold': set(edata.row_labels), #None means 'use all genes'. This set only used when outputting filter
            'HDT': set([edata.row_labels[i] for i,x in enumerate(hdt_mask) if x]),
            'Fano': set([edata.row_labels[i] for i,x in enumerate(fano_mask) if x])
        });

    edata.filters = filter_dict;

    Models = dict();
    emodel = dict({"Data": edata});
    Models.update({"Expression": emodel});

    #%% Probability transform
    if(not options.nomodel):

        FP_Output('\nFitting expression data to exp/norm mixture model');
        (pdata, mu_h, mu_l, st_h, Pi) = Transforms.probability_of_expression(edata);
        prob_params = (mu_h, mu_l, st_h, Pi);


        FP_Output('\nCorrecting for false-negatives using housekeeping gene levels');
        (fit_func, params) = Transforms.create_false_neg_map(original_data, housekeeping_filename);
        (pdata, fn_prob) = Transforms.correct_for_fn(pdata, mu_h, fit_func, params);
        fn_prob[edata > 0] = 0;

        pdata = ProbabilityData(pdata, edata);
        pmodel = dict({"Data": pdata});
        Models.update({"Probability": pmodel});

        edata.weights = 1-fn_prob;
        pdata.weights = 1-fn_prob;

        sample_passes_qc, sample_qc_scores = Transforms.quality_check(params);
        FileIO.write_qc_file(dir_name, sample_passes_qc, sample_qc_scores, edata.col_labels);

        #If specified, remove items that did not pass qc check
        if(options.qc):
            for name, dataMatrix in Models.items():
                Models[name] = dataMatrix.subset_samples(sample_passes_qc);

            sample_qc_scores = sample_qc_scores[sample_passes_qc];
    else:
        sample_qc_scores = None;
        prob_params = None;


    Transforms.z_normalize(edata);

    fout_js = HtmlViewer.get_output_js_handle(dir_name);

    for name, model in Models.items():
        model_dir = os.path.join(dir_name, name);
        data = model["Data"];
        try:
            os.makedirs(model_dir);
        except OSError:
            pass;

        FP_Output('\nModel: ', name)

        #Evaluate Signatures
        FP_Output("\nEvaluating signature scores on samples...");

        sig_scores_dict = dict();

        pbar = ProgressBar(len(sigs));
        for sig in sigs:
            try:
                sig_scores_dict[sig.name] = data.eval_signature(sig);
            except ValueError:  #Only thrown when the signature has no genes in the data
                pass #Just discard the signature then
            pbar.update();
        pbar.complete();

        if(options.precomputed):
            precomputed_sig_scores = Signatures.load_precomputed(options.precomputed, data.col_labels);
            sig_scores_dict.update(precomputed_sig_scores);

        #Adds in quality score as a pre-computed signature
        if(sample_qc_scores is not None): #Might be None if --nomodel option is selected
            sig_scores_dict["FP_Quality"] = Signatures.SignatureScores(sample_qc_scores,"FP_Quality",data.col_labels,isFactor=False, isPrecomputed=True);

        model["signatureScores"] = sig_scores_dict;
        model["projectionData"] = [];
        model["sampleLabels"] = data.col_labels;

        for filter_name in filter_dict.keys():
            projData = dict();

            FP_Output("\nFilter-Level:", filter_name);
            #%% Dimensional Reduction procedures
            FP_Output("\nProjecting data into 2 dimensions");

            projections, pcdata = Projections.generate_projections(data, filter_name);

            #Evaluate Clusters
            FP_Output("Evaluating Clusters...");
            clusters = Projections.define_clusters(projections);

            #%% Evaluating signatures against projections
            sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections(projections, sig_scores_dict);

            #Store in projData
            projData["filter"] = filter_name;
            projData["genes"] = data.filters[filter_name];
            projData["pca"] = False;
            projData["projections"] = projections;
            projData["sigProjMatrix"] = sig_proj_matrix;
            projData["sigProjMatrix_p"] = sig_proj_matrix_p;
            projData["projectionKeys"] = sp_col_labels;
            projData["signatureKeys"] = sp_row_labels;
            projData["clusters"] = clusters;
            model["projectionData"].append(projData);


            #Now do it all again using the principal component data
            if(options.pca_filter):
                pcdata = Projections.filter_PCA(pcdata, scores=sample_qc_scores, variance_proportion=0.25);
            else:
                pcdata = Projections.filter_PCA(pcdata, variance_proportion=0.25, min_components = 30);

            #%% Dimensional Reduction procedures
            FP_Output("Projecting PC data into 2 dimensions");

            projections, pcdata2 = Projections.generate_projections(pcdata, filter_name);

            #Evaluate Clusters
            FP_Output("Evaluating Clusters...");
            clusters = Projections.define_clusters(projections);

            #%% Evaluating signatures against projections
            sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections(projections, sig_scores_dict);

            projData = dict();
            projData["filter"] = filter_name;
            projData["genes"] = data.filters[filter_name];
            projData["pca"] = True;
            projData["projections"] = projections;
            projData["sigProjMatrix"] = sig_proj_matrix;
            projData["sigProjMatrix_p"] = sig_proj_matrix_p;
            projData["projectionKeys"] = sp_col_labels;
            projData["signatureKeys"] = sp_row_labels;
            projData["clusters"] = clusters;
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
        OUTPUT_SIGNATURE_LIMIT = 200;
        if(all_data.shape[1] > 2000 and len(signatureScores) > OUTPUT_SIGNATURE_LIMIT):
            aa = np.argsort(signature_significance);
            threshold = signature_significance[aa[OUTPUT_SIGNATURE_LIMIT]];
        else:
            threshold = -1.301;

        #Iterate back through and prune signatures worse than threshold
        #Create a dictionary of sigs to keep
        keep_sig = dict();
        for name, sig_score in signatureScores.items():
            if(not sig_score.isPrecomputed and signature_significance[name] > threshold):
                keep_sig.update({name: False});
            else:
                keep_sig.update({name: True});

        #Remove values in the model's signatureScores dict
        for name, sig_score in signatureScores.items():
            if(keep_sig[name] == False):
                signatureScores.pop(name);

        #Remove values in each filters sigProjMatrix and the sigProjMatrix keys
        for projData in model['projectionData']:
            projData_keep_sig = np.array([keep_sig[sig] for sig in projData["signatureKeys"]]);
            projData["sigProjMatrix"] = projData["sigProjMatrix"][projData_keep_sig,:];
            projData["sigProjMatrix_p"] = projData["sigProjMatrix_p"][projData_keep_sig,:];
            projData["signatureKeys"] = [x for x in projData["signatureKeys"] if keep_sig[x]];



    #Write signatures to file
    #Assemble signatures into an object, then convert to JSON variable and write
    sig_dict = {};
    for sig in sigs:
        sig_genes = sig.sig_dict.keys();
        sig_values = sig.sig_dict.values();
        sort_i = np.array(sig_values).argsort()[::-1];#Put positive signatures first
        sig_genes = [sig_genes[i] for i in sort_i];
        sig_values = [sig_values[i] for i in sort_i];
        sig_dict.update({sig.name: {'Genes':sig_genes, 'Signs':sig_values}});
    fout_js.write(HtmlViewer.toJS_variable("FP_Signatures", sig_dict));

    #Merge all the holdouts back into the model
    if(options.subsample_size is not None):
        FP_Output("Merging held-out samples back in")
        SubSample.merge_samples(all_data, Models, sigs, prob_params);

    #Write the original data matrix to the javascript file.
    #First, cluster genes
    from scipy.cluster.hierarchy import leaves_list, linkage;
    edata = Models["Expression"]["Data"];
    linkage_matrix = linkage(edata);
    leaves_i = leaves_list(linkage_matrix);
    edata_clustered = edata[leaves_i, :];
    edata_clustered.row_labels = [edata.row_labels[i] for i in leaves_i];

    data_json = dict({
        'data': edata_clustered,
        'gene_labels': edata_clustered.row_labels,
        'sample_labels': edata_clustered.col_labels,
    });
    fout_js.write(HtmlViewer.toJS_variable("FP_ExpressionMatrix", data_json));

    FileIO.write_models(dir_name, Models);

    #Remove model["Data"] since only the expression data is written for JS
    #  and it's already written above in FP_ExpressionMatrix
    for name, model in Models.items():
        model.pop("Data");

    fout_js.write(HtmlViewer.toJS_variable("FP_Models", Models));

    fout_js.close();
    HtmlViewer.copy_html_files(dir_name);
    FP_Output("\nFastProject Analysis Complete")
    elapsed_time = time.time() - start_time;
    FP_Output("Elapsed Time {:.2f} seconds".format(elapsed_time));
