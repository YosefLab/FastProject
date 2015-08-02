from __future__ import division, print_function;
import os;
import numpy as np;
import time;
import scipy.stats;
from FastProject import Filters;
from FastProject import FileIO;
from FastProject import Transforms;
from FastProject import Signatures;
from FastProject import Projections;
from FastProject.DataTypes import ExpressionData, ProbabilityData, PCData;
from FastProject.Utils import ProgressBar;
from FastProject import HtmlViewer;

def FullOutput(options, args):

    start_time = time.time();
    if(options.housekeeping):
        housekeeping_filename = options.housekeeping;
    else:
        housekeeping_filename = '';  #If not in interactive mode and no -h specified, just use the default list

    def get_housekeeping_file():
        fn = housekeeping_filename;
        return fn;

    #%% Read expression data from file
    if(len(args) > 0):
        filename = args[0];
    else:
        raise ValueError("Argument Error:  data_file not specified.\nExiting...");

    if(not os.path.isfile(filename)):
        raise ValueError("Argument Error: data file not found.\nExiting...");

    (edata, genes, cells) = FileIO.read_matrix(filename);

    print("Imported ", edata.shape[0], " genes across ", edata.shape[1], " samples");

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


    #Wrap data in ExpressionData object
    edata = ExpressionData(edata, genes, cells);

    if(options.subsample_size > edata.shape[1]):
        options.subsample_size = None;

    #Hold on to originals so we don't lose data after filtering in case it's needed later
    original_data = edata.copy();

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

    if(not os.path.isdir(dir_name)):
        os.makedirs(dir_name);

    #Filtering
    filter_dict = {};
    if(options.nofilter):
        edata = Filters.filter_genes_novar(edata);

        filter_dict.update({'None': set()});
    else:
        edata = Filters.filter_genes_threshold(edata, 0.2);

        if(options.subsample_size):
            sub_ii = np.random.choice(edata.shape[1], options.subsample_size, replace=False);
            fdata = edata.subset_samples(sub_ii);
        else:
            fdata = edata;

        #HDT Filtering
        print("Removing genes with unimodal distribution across samples using Hartigans DT...");
        hdt_mask = Filters.filter_genes_hdt(fdata, 0.05);
        #Fano Filtering
        print("Applying Fano-Filtering...");
        fano_mask = Filters.filter_genes_fano(fdata, 2);

        filter_dict.update({
            'None': set(edata.row_labels), #None means 'use all genes'. This set only used when outputting filter
            'HDT': set([edata.row_labels[i] for i,x in enumerate(hdt_mask) if x]),
            'Fano': set([edata.row_labels[i] for i,x in enumerate(fano_mask) if x])
        });

    edata.filters = filter_dict;

    #%% Probability transform
    housekeeping_filename = get_housekeeping_file();

    print()
    print('Fitting expression data to exp/norm mixture model');
    (pdata, mu_h) = Transforms.probability_of_expression(edata, options.subsample_size);


    print();
    print('Correcting for false-negatives using housekeeping gene levels');
    (fit_func, params) = Transforms.create_false_neg_map(original_data, housekeeping_filename);
    (pdata, fn_prob) = Transforms.correct_for_fn(pdata, mu_h, fit_func, params);
    fn_prob[edata > 0] = 0;

    pdata = ProbabilityData(pdata, edata);

    edata.weights = 1-fn_prob;
    pdata.weights = 1-fn_prob;

    sample_passes, sample_scores = Transforms.quality_check(params);
    FileIO.write_qc_file(os.path.join(dir_name, "QC_Report.html"), sample_passes, sample_scores, edata.col_labels);

    if(options.qc):
        pdata = pdata.subset_samples(sample_passes);
        edata = edata.subset_samples(sample_passes);
        sample_scores = sample_scores[sample_passes];

    if(options.subsample_size > edata.shape[1]):
        options.subsample_size = None;

    Transforms.z_normalize(edata);

    model_names = ['Expression', 'Probability'];
    model_data = [edata, pdata];

    fout_js = open(dir_name + os.sep + "FP_data.jsdata", 'w');
    js_models = [];

    for name, data in zip(model_names, model_data):
        model_dir = os.path.join(dir_name, name);
        try:
            os.makedirs(model_dir);
        except OSError:
            pass;

        print();
        print('Model: ', name)

        #Evaluate Signatures
        print();
        print("Evaluating signature scores on samples...");

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

        #Prompt to save data
        out_file = 'SignatureScores.txt';
        FileIO.write_signature_scores(os.path.join(model_dir, out_file), sig_scores_dict, data.col_labels);

        #Save data to js model as well
        js_model_dict = {'model': name};
        js_model_dict.update({'signatureScores': sig_scores_dict})
        js_model_dict.update({'sampleLabels': data.col_labels});
        js_model_dict.update({'projectionData': []})
        js_models.append(js_model_dict);

        for filter_name in filter_dict.keys():
            if(filter_name == "None"):
                filter_dir = os.path.join(model_dir, "No_Filter");
            else:
                filter_dir = os.path.join(model_dir, filter_name + "_Filter");
            try:
                os.makedirs(filter_dir);
            except OSError:
                pass;

            print();
            print("Filter-Level:", filter_name);
            #%% Dimensional Reduction procedures
            print();
            print("Projecting data into 2 dimensions");

            projections, pcdata = Projections.generate_projections(data, filter_name, options.subsample_size);

            #Evaluate Clusters
            print("Evaluating Clusters...");
            clusters = Projections.define_clusters(projections);

            #%% Evaluating signatures against projections
            sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections(projections, sig_scores_dict ,subsample_size=options.subsample_size);

            #Save Projections
            FileIO.write_projection_file(os.path.join(filter_dir, 'Projections.txt'), data.col_labels, projections);

            #Save Clusters
            FileIO.write_cluster_file(os.path.join(filter_dir, 'Clusters.txt'), data.col_labels, clusters)

            #Output matrix of p-values for conformity scores
            FileIO.write_matrix(os.path.join(filter_dir, "DissimilarityMatrix.txt"),sig_proj_matrix, sp_row_labels, sp_col_labels);
            FileIO.write_matrix(os.path.join(filter_dir, "PMatrix.txt"),sig_proj_matrix_p, sp_row_labels, sp_col_labels);

            #Output genes used in filter
            FileIO.write_filter_file(os.path.join(filter_dir, 'ProjectedGenes.txt'), data.filters[filter_name]);

            #Output JS
            js_filt_dict = dict();
            js_model_dict['projectionData'].append(js_filt_dict);
            js_filt_dict.update({'filter': filter_name});
            js_filt_dict.update({'genes': data.filters[filter_name]});
            js_filt_dict.update({'pca': False});
            js_filt_dict.update({'projections': projections});
            js_filt_dict.update({'sigProjMatrix': sig_proj_matrix});
            js_filt_dict.update({'sigProjMatrix_p': sig_proj_matrix_p});
            js_filt_dict.update({'projectionKeys': sp_col_labels});
            js_filt_dict.update({'signatureKeys': sp_row_labels});
            js_filt_dict.update({'clusters': clusters});



            #Now do it all again using the principal component data
            if(options.pca_filter):
                pcdata = Projections.filter_PCA(pcdata, scores=sample_scores, variance_proportion=0.25);
            else:
                pcdata = Projections.filter_PCA(pcdata, variance_proportion=0.25, min_components = 30);

            #%% Dimensional Reduction procedures
            print();
            print("Projecting data into 2 dimensions");

            projections, pcdata2 = Projections.generate_projections(pcdata, filter_name, options.subsample_size);

            #Evaluate Clusters
            print("Evaluating Clusters...");
            clusters = Projections.define_clusters(projections);

            #%% Evaluating signatures against projections
            sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections(projections, sig_scores_dict, subsample_size = options.subsample_size);

            #Save Projections
            FileIO.write_projection_file(os.path.join(filter_dir, 'Projections-PC.txt'), pcdata.col_labels, projections);

            #Save Clusters
            FileIO.write_cluster_file(os.path.join(filter_dir, 'Clusters-PC.txt'), pcdata.col_labels, clusters)

            #Output matrix of p-values for conformity scores
            FileIO.write_matrix(os.path.join(filter_dir, "DissimilarityMatrix-PC.txt"),sig_proj_matrix, sp_row_labels, sp_col_labels);
            FileIO.write_matrix(os.path.join(filter_dir, "PMatrix-PC.txt"),sig_proj_matrix_p, sp_row_labels, sp_col_labels);

            #Output JS
            js_filt_dict = dict();
            js_model_dict['projectionData'].append(js_filt_dict);
            js_filt_dict.update({'filter': filter_name});
            js_filt_dict.update({'genes': data.filters[filter_name]});
            js_filt_dict.update({'pca': True});
            js_filt_dict.update({'projections': projections});
            js_filt_dict.update({'sigProjMatrix': sig_proj_matrix});
            js_filt_dict.update({'sigProjMatrix_p': sig_proj_matrix_p});
            js_filt_dict.update({'projectionKeys': sp_col_labels});
            js_filt_dict.update({'signatureKeys': sp_row_labels});
            js_filt_dict.update({'clusters': clusters});


    fout_js.write(HtmlViewer.toJS_variable("FP_Models", js_models));

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

    #Write the original data matrix to the javascript file.
    #First, cluster genes
    from scipy.cluster.hierarchy import leaves_list, linkage;
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

    fout_js.close();
    HtmlViewer.copy_html_files(dir_name);
    print();
    print("FastProject Analysis Complete")
    elapsed_time = time.time() - start_time;
    print("Elapsed Time {:.2f} seconds".format(elapsed_time));
