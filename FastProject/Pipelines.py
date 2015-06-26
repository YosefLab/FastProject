from __future__ import division, print_function;
import sys;
import os;
import numpy as np;
import time;
from FastProject import Filters;
from FastProject import FileIO;
from FastProject import Transforms;
from FastProject import Signatures;
from FastProject import Projections;
from FastProject.DataTypes import ExpressionData, ProbabilityData, PCData;
from FastProject.Utils import ProgressBar;
from FastProject import HtmlViewer;

def SingleOutput(options, args):

    if(options.housekeeping):
        housekeeping_filename = options.housekeeping;
    else:
        if(options.interactive):
            housekeeping_filename = "NOT_SPECIFIED";
        else:
            housekeeping_filename = '';  #If not in interactive mode and no -h specified, just use the default list

    def get_housekeeping_file():
        fn = housekeeping_filename;
        if(fn == "NOT_SPECIFIED"):
            while(True):
                fn = raw_input("Enter name of housekeeping file, or press ENTER to use default list: ");

                if(fn == '' or os.path.isfile(fn)):
                    break;
                else:
                    print("Error : File not found");
                    continue;

        return fn;

    #%% Read expression data from file
    while(True):
        if(options.interactive):
            filename = raw_input("Enter name of data file: ");
        else:
            if(len(args) > 0):
                filename = args[0];
                if(not os.path.isfile(filename)):
                    raise ValueError("Argument Error: data file not found.\nExiting...");
            else:
                raise ValueError("Argument Error:  data_file not specified.\nExiting...");

        if(os.path.isfile(filename)):
            (edata, genes, cells) = FileIO.read_matrix(filename);
            break;
        else:
            print("Error : File not found");
            continue

    print("Imported ", edata.shape[0], " genes across ", edata.shape[1], " samples");

    #Wrap data in ExpressionData object
    edata = ExpressionData(edata, genes, cells);

    #Hold on to originals so we don't lose data after filtering in case it's needed later
    original_data = edata;

    #Create directory for all outputs
    tt = time.localtime();
    dot_index = filename.rfind('.');
    if(dot_index > -1):
        dir_name = filename[0:dot_index];
    else:
        dir_name = filename;

    dir_name = dir_name + '_' + str(tt.tm_year) + '{:0>2d}'.format(tt.tm_mon) + '{:0>2d}'.format(tt.tm_mday);
    dir_name = dir_name + '_' + '{:0>2d}'.format(tt.tm_hour) + '{:0>2d}'.format(tt.tm_min) + '{:0>2d}'.format(tt.tm_sec);

    os.makedirs(dir_name);

    #%% Filtering genes
    while(True):  #Loop exited with 'break', see below

        original = edata.shape;
        if(options.interactive):
            print();
            print("Filtering Options Available");
            print();
            print("\t0.  Continue")
            print("\t1.  Remove housekeeping genes");
            print("\t2.  Remove inactive genes");
            print("\t3.  Filter genes for biomadility using HDT");
            print("\t4.  Apply Fano-factor filtering.");
            print();

            choice = raw_input("Make a selection: ");

            try:
                choice = int(choice);
            except ValueError:
                print("Error : Bad value.  Enter a number")
                continue;
        else:
            if(len(options.filters) == 0):
                choice = 0;
            else:
                choice = options.filters[0];
                options.filters = options.filters[1:];
                try:
                    choice = int(choice);
                    if(choice > 4):
                        raise ValueError;
                except ValueError:
                    print("Error:  Bad option for -f, --filters.  Exiting");
                    sys.exit();
        print();
        if(choice==0):
            break;
        elif(choice==1): #Housekeeping
            print("Removing housekeeping genes...");
            housekeeping_filename = get_housekeeping_file();
            edata = Filters.filter_housekeeping(edata, housekeeping_filename);
        elif(choice==2): #Threshold of activation
            print("Removing genes inactive in > 80% samples...");
            edata = Filters.filter_genes_threshold(edata, 0.2);
        elif(choice==3): #HDT test
            print("Removing genes with unimodal distribution across samples using Hartigans DT...");
            gene_passes = Filters.filter_genes_hdt(edata, 0.05, return_mask=True);
            edata.projection_mask = np.logical_and(edata.projection_mask, gene_passes);
        elif(choice==4): #Save to file
            print("Applying Fano-Filtering...");
            gene_passes = Filters.filter_genes_fano(edata, 2);
            edata.projection_mask = np.logical_and(edata.projection_mask, gene_passes);
        else:
            print("Error : Invalid Choice\n");
            continue;

        if(0 < choice < 3):
            print("Removed ", original[0]-edata.shape[0], " Genes");
            print(edata.shape[0], " Genes retained");
        if(2 < choice < 5):
            print(np.sum(edata.projection_mask), "Genes saved for projection.")

    data = edata;  #'data' is used for projections/signatures.  Can be overwritten with the probability data object

    #%% Probability transform
    housekeeping_filename = get_housekeeping_file();

    print()
    print('Fitting expression data to exp/norm mixture model');
    (prob, mu_h) = Transforms.probability_of_expression(data);


    print();
    print('Correcting for false-negatives using housekeeping gene levels');
    (fit_func, params) = Transforms.create_false_neg_map(original_data, housekeeping_filename);
    (prob, fn_prob) = Transforms.correct_for_fn(prob, mu_h, fit_func, params);
    fn_prob[data > 0] = 0;

    prob = ProbabilityData(prob, data);

    data.weights = 1-fn_prob;
    prob.weights = 1-fn_prob;

    sample_passes, sample_scores = Transforms.quality_check(params);

    if(options.qc):
        prob = prob.subset_samples(sample_passes);
        data = data.subset_samples(sample_passes);
        sample_scores = sample_scores[sample_passes];

    Transforms.z_normalize(data);

    while(True):

        if(options.interactive and (not options.probability)):
            choice = raw_input("Transform data to probability values? [y/n]: ");
        else:
            choice = 'y' if options.probability else 'n';

        if(choice.lower()[0] == 'y'):
            data = prob;
            break;

        elif(choice.lower()[0] == 'n'):
            break;

        else:
            print("Error : invalid input.  Enter Y or N");
            continue;

    #%% PCA transform?
    while(True):

        if(options.interactive and (not options.pca)):
            choice = raw_input("Perform PCA before projection? [y/n]: ");
        else:
            choice = 'y' if options.pca else 'n';

        if(choice.lower()[0] == 'y'):
            #Transform into top N principal components
            pc_data = Projections.perform_weighted_PCA(data);
            pc_data = Projections.filter_PCA(pc_data, scores=sample_scores, variance_proportion=0.25);
            data = pc_data;
            break;
        elif(choice.lower()[0] == 'n'):
            break;
        else:
            print("Error : invalid input.  Enter Y or N");
            continue;

    #%% Prompt: Save transformed data to file
    while(True):

        if(options.interactive and (not options.output)):
            print();
            choice = raw_input("Save transformed data to file? [y/n]: ");
        else:
            choice = 'y' if options.output else 'n';

        if(choice.lower()[0] == 'y'):
            #Save data
            lastdot = filename.rfind('.');
            if(lastdot > -1):
                out_file = filename[0:lastdot] + '_xfrm' + filename[lastdot:];
            else:
                out_file = filename + '_xfrm';

            FileIO.write_data(dir_name + os.sep + out_file, data)
            break;
        elif(choice.lower()[0] == 'n'):
            break;
        else:
            print("Error : invalid input.  Enter Y or N");
            continue;

    #%% Dimensional Reduction procedures
    print();
    print("Projecting data into 2 dimensions");

    projections = Projections.generate_projections(data);

    #Save data
    out_file = 'projections.txt'
    Projections.write_projection_file(dir_name + os.sep + out_file, data.col_labels, projections);


    #%% Output Projection Plots
    print();
    print('Outputting Projection Plots to File');
    pp = ProgressBar(len(projections));
    for proj_name in projections:
        FileIO.write_scatter_plot(
            filename = dir_name + os.sep + proj_name,
            x_coords = projections[proj_name][0,:],
            y_coords = projections[proj_name][1,:],
            xlabel = 'Dim 1',
            ylabel = 'Dim 2',
            title = proj_name);
        pp.update();
    pp.complete();

    #%% Signature file
    sigs = [];
    USE_SIGNATURES = False;
    #Use signature?
    while(True):

        if(options.interactive and (not options.signatures)):
            print();
            print("Enter name of signature file or press enter to omit signature analysis.")
            choice = raw_input("File name : ");
        else:
            if(options.signatures):
                choice = options.signatures;
            else:
                choice = '';

        if(len(choice) == 0):
            break;

        elif(os.path.isfile(choice)): #Load signature
            USE_SIGNATURES = True;
            sigs = Signatures.read_signatures(choice);
            break;
        else:
            print("Error : File not found");
            continue;

    #Filter Signatures to use
    if(USE_SIGNATURES and options.interactive):
        while(True):
            ##List signatures
            print("\n" + str(len(sigs)) + " Signatures Loaded: " + "\n");
            for i,sig in enumerate(sigs):
                print(sig.name);
                if(i == 50):
                    print("...rest of output suppressed, " + str(len(sigs)-50) + " signatures remaining")
                    break;

            print();
            print("Filter signatures using keywords.");
            print("Signatures containing any word in the list will be retained.");
            print("Enter a list of keywords (comma seperated), or hit enter to continue");
            keywords = raw_input("Keywords: ");

            if(len(keywords)==0):
                break;
            else:
                #Filter based on keywords
                keywords = [word.strip() for word in keywords.split(',')];
                sigs = Signatures.filter_sig_list(sigs, keywords);

    #Evaluate Signatures
    if(USE_SIGNATURES):
        print();
        print("Evaluating signature scores on samples...");

        sig_scores = dict();

        pbar = ProgressBar(len(sigs));
        for sig in sigs:
            try:
                sig_scores[sig.name] = data.eval_signature(sig);
            except ValueError:  #Only thrown when the signature has no genes in the data
                pass #Just discard the signature then
            pbar.update();
        pbar.complete();

        if(options.precomputed):
            precomputed_sig_scores = Signatures.load_precomputed(options.precomputed, data.col_labels);
            sig_scores.update(precomputed_sig_scores);

        #Prompt to save data
        out_file = 'sig_scores.txt';
        FileIO.write_signature_scores(dir_name + os.sep + out_file, sig_scores, data.col_labels);



    #%% Evaluating signatures against projections
    if(USE_SIGNATURES):
        sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections_v2(projections, sig_scores);


        #Output matrix of p-values for conformity scores
        FileIO.write_matrix(dir_name + os.sep + "p_matrix.txt",sig_proj_matrix_p, sp_row_labels, sp_col_labels);


        #%% Output top N plots
        N_PLOTS = 30;
        flat_indices = np.argsort(sig_proj_matrix_p, axis=None);
        row_indices, col_indices = np.unravel_index(flat_indices, sig_proj_matrix_p.shape);

        print();
        print('Outputting Projection-Signature Plots to File');
        pp = ProgressBar(N_PLOTS);
        for i in range(N_PLOTS):
            r = row_indices[i];
            c = col_indices[i];

            sig_name = sp_row_labels[r];
            sig_name = sig_name.replace("\\","_"); #Make valid file name
            sig_name = sig_name.replace("/", "_");
            proj_name = sp_col_labels[c];

            FileIO.write_scatter_plot(
                filename = dir_name + os.sep + sig_name+'_'+proj_name,
                x_coords = projections[proj_name][0,:],
                y_coords = projections[proj_name][1,:],
                colors   = sig_scores[sig_name],
                title    = proj_name + '\n' + sig_name + '\np = ' + '{:.3f}'.format(sig_proj_matrix_p[r,c]));

            pp.update();

        pp.complete();

        #Output to the viewer
        from FastProject import HtmlViewer;

        if(type(data) is PCData and type(data.parent_data) is ExpressionData):
            label = "ExpressionPD";
        elif(type(data) is PCData and type(data.parent_data) is ProbabilityData):
            label = "ProbabilityPC";
        elif(type(data) is ExpressionData):
            label = "Expression";
        elif(type(data) is ProbabilityData):
            label = "Probability";
        else:
            raise Exception("Unrecognized Data Type");

        HtmlViewer.copy_html_files(dir_name);
        #Wrap data into an object
        js_out = dict();
        js_out.update({'Projections': projections});
        js_out.update({'SigScores': sig_scores});
        js_out.update({'SigProjMatrix': sig_proj_matrix});
        js_out.update({'SigProjMatrix_p': sig_proj_matrix_p});
        js_out.update({'ProjectionKeys': sp_col_labels});
        js_out.update({'SignatureKeys': sp_row_labels});

        #Ouput object to js file
        fout_js = open(dir_name + os.sep + "FP_data.js", 'w');
        fout_js.write(HtmlViewer.toJS_variable("FP_" + label, js_out));
        fout_js.close();


    print();
    print("FastProject Analysis Complete")

def FullOutput(options, args):
    start_time = time.time();
    if(options.housekeeping):
        housekeeping_filename = options.housekeeping;
    else:
        if(options.interactive):
            housekeeping_filename = "NOT_SPECIFIED";
        else:
            housekeeping_filename = '';  #If not in interactive mode and no -h specified, just use the default list

    def get_housekeeping_file():
        fn = housekeeping_filename;
        if(fn == "NOT_SPECIFIED"):
            while(True):
                fn = raw_input("Enter name of housekeeping file, or press ENTER to use default list: ");

                if(fn == '' or os.path.isfile(fn)):
                    break;
                else:
                    print("Error : File not found");
                    continue;
        return fn;

    #%% Read expression data from file
    while(True):
        if(options.interactive):
            filename = raw_input("Enter name of data file: ");
        else:
            if(len(args) > 0):
                filename = args[0];
                if(not os.path.isfile(filename)):
                    raise ValueError("Argument Error: data file not found.\nExiting...");
            else:
                raise ValueError("Argument Error:  data_file not specified.\nExiting...");

        if(os.path.isfile(filename)):
            (edata, genes, cells) = FileIO.read_matrix(filename);
            break;
        else:
            print("Error : File not found");
            continue

    print("Imported ", edata.shape[0], " genes across ", edata.shape[1], " samples");

    #Wrap data in ExpressionData object
    edata = ExpressionData(edata, genes, cells);

    #Hold on to originals so we don't lose data after filtering in case it's needed later
    original_data = edata;

    #Create directory for all outputs
    tt = time.localtime();
    dot_index = filename.rfind('.');
    if(dot_index > -1):
        dir_name = filename[0:dot_index];
    else:
        dir_name = filename;

    dir_name = dir_name + '_' + str(tt.tm_year) + '{:0>2d}'.format(tt.tm_mon) + '{:0>2d}'.format(tt.tm_mday);
    dir_name = dir_name + '_' + '{:0>2d}'.format(tt.tm_hour) + '{:0>2d}'.format(tt.tm_min) + '{:0>2d}'.format(tt.tm_sec);

    os.makedirs(dir_name);

    #%% Filtering genes
    while(True):  #Loop exited with 'break', see below

        original = edata.shape;
        if(options.interactive):
            print();
            print("Filtering Options Available");
            print();
            print("\t0.  Continue")
            print("\t1.  Remove housekeeping genes");
            print("\t2.  Remove inactive genes");
            print("\t3.  Filter genes for biomadility using HDT");
            print("\t4.  Apply Fano-factor filtering.");
            print();

            choice = raw_input("Make a selection: ");

            try:
                choice = int(choice);
            except ValueError:
                print("Error : Bad value.  Enter a number")
                continue;
        else:
            if(len(options.filters) == 0):
                choice = 0;
            else:
                choice = options.filters[0];
                options.filters = options.filters[1:];
                try:
                    choice = int(choice);
                    if(choice > 4):
                        raise ValueError;
                except ValueError:
                    print("Error:  Bad option for -f, --filters.  Exiting");
                    sys.exit();
        print();
        if(choice==0):
            break;
        elif(choice==1): #Housekeeping
            print("Removing housekeeping genes...");
            housekeeping_filename = get_housekeeping_file();
            edata = Filters.filter_housekeeping(edata, housekeeping_filename);
        elif(choice==2): #Threshold of activation
            print("Removing genes inactive in > 80% samples...");
            edata = Filters.filter_genes_threshold(edata, 0.2);
        elif(choice==3): #HDT test
            print("Removing genes with unimodal distribution across samples using Hartigans DT...");
            gene_passes = Filters.filter_genes_hdt(edata, 0.05, return_mask=True);
            edata.projection_mask = np.logical_and(edata.projection_mask, gene_passes);
        elif(choice==4): #Save to file
            print("Applying Fano-Filtering...");
            gene_passes = Filters.filter_genes_fano(edata, 2);
            edata.projection_mask = np.logical_and(edata.projection_mask, gene_passes);
        else:
            print("Error : Invalid Choice\n");
            continue;

        if(0 < choice < 3):
            print("Removed ", original[0]-edata.shape[0], " Genes");
            print(edata.shape[0], " Genes retained");
        if(2 < choice < 5):
            print(np.sum(edata.projection_mask), "Genes saved for projection.")

    data = edata;  #'data' is used for projections/signatures.  Can be overwritten with the probability data object

    #%% Probability transform
    housekeeping_filename = get_housekeeping_file();

    print()
    print('Fitting expression data to exp/norm mixture model');
    (prob, mu_h) = Transforms.probability_of_expression(data);


    print();
    print('Correcting for false-negatives using housekeeping gene levels');
    (fit_func, params) = Transforms.create_false_neg_map(original_data, housekeeping_filename);
    (prob, fn_prob) = Transforms.correct_for_fn(prob, mu_h, fit_func, params);
    fn_prob[data > 0] = 0;

    prob = ProbabilityData(prob, data);

    data.weights = 1-fn_prob;
    prob.weights = 1-fn_prob;

    sample_passes, sample_scores = Transforms.quality_check(params);

    if(options.qc):
        prob = prob.subset_samples(sample_passes);
        data = data.subset_samples(sample_passes);
        sample_scores = sample_scores[sample_passes];

    Transforms.z_normalize(data);

    #Perform PCA on the data, with and without the probability xfrm
    pc_data = Projections.perform_weighted_PCA(data);
    pc_prob = Projections.perform_weighted_PCA(prob);

    if(options.pca_filter):
        pc_data = Projections.filter_PCA(pc_data, scores=sample_scores, variance_proportion=0.25);
        pc_prob = Projections.filter_PCA(pc_prob, scores=sample_scores, variance_proportion=0.25);
    else:
        pc_data = Projections.filter_PCA(pc_data, variance_proportion=0.25, min_components = 30);
        pc_prob = Projections.filter_PCA(pc_prob, variance_proportion=0.25, min_components = 30);



    #%% Signature file
    sigs = [];
    #Use signature?
    while(True):

        if(options.interactive and (not options.signatures)):
            print();
            print("Enter name of signature file...")
            choice = raw_input("File name : ");
        else:
            if(options.signatures):
                choice = options.signatures;
            else:
                choice = '';

        if(len(choice) == 0):
            break;

        elif(os.path.isfile(choice)): #Load signature
            sigs = Signatures.read_signatures(choice);
            break;
        else:
            print("Error : File not found");
            continue;

    #Filter Signatures to use
    if(options.interactive):
        while(True):
            ##List signatures
            print("\n" + str(len(sigs)) + " Signatures Loaded: " + "\n");
            for i,sig in enumerate(sigs):
                print(sig.name);
                if(i == 50):
                    print("...rest of output suppressed, " + str(len(sigs)-50) + " signatures remaining")
                    break;

            print();
            print("Filter signatures using keywords.");
            print("Signatures containing any word in the list will be retained.");
            print("Enter a list of keywords (comma separated), or hit enter to continue");
            keywords = raw_input("Keywords: ");

            if(len(keywords)==0):
                break;
            else:
                #Filter based on keywords
                keywords = [word.strip() for word in keywords.split(',')];
                sigs = Signatures.filter_sig_list(sigs, keywords);


    data_matrices = [data, prob, pc_data, pc_prob];
    data_labels = ["Expression", "Probability", "ExpressionPC", "ProbabilityPC"];

    fout_js = open(dir_name + os.sep + "FP_data.jsdata", 'w');

    for label, data in zip(data_labels, data_matrices):
        print();
        print("Evaluating: ", label);
        #%% Dimensional Reduction procedures
        print();
        print("Projecting data into 2 dimensions");

        projections = Projections.generate_projections(data);

        #Save data
        out_file = label + '_Projections.txt'
        Projections.write_projection_file(dir_name + os.sep + out_file, data.col_labels, projections);

        #Evaluate Signatures
        print();
        print("Evaluating signature scores on samples...");

        sig_scores = dict();

        pbar = ProgressBar(len(sigs));
        for sig in sigs:
            try:
                sig_scores[sig.name] = data.eval_signature(sig);
            except ValueError:  #Only thrown when the signature has no genes in the data
                pass #Just discard the signature then
            pbar.update();
        pbar.complete();

        if(options.precomputed):
            precomputed_sig_scores = Signatures.load_precomputed(options.precomputed, data.col_labels);
            sig_scores.update(precomputed_sig_scores);

        #Prompt to save data
        out_file = label + '_SigScores.txt';
        FileIO.write_signature_scores(dir_name + os.sep + out_file, sig_scores, data.col_labels);



        #%% Evaluating signatures against projections
        sp_row_labels, sp_col_labels, sig_proj_matrix, sig_proj_matrix_p = Signatures.sigs_vs_projections_v2(projections, sig_scores);


        #Evaluate Clusters
        print("Evaluating Clusters...");
        clusters = Projections.define_clusters(projections);

        #Output matrix of p-values for conformity scores
        FileIO.write_matrix(dir_name + os.sep + label + "_DissimilarityMatrix.txt",sig_proj_matrix, sp_row_labels, sp_col_labels);
        FileIO.write_matrix(dir_name + os.sep + label + "_PMatrix.txt",sig_proj_matrix_p, sp_row_labels, sp_col_labels);

        #Wrap data into an object
        js_out = dict();
        js_out.update({'Projections': projections});
        js_out.update({'SigScores': sig_scores});
        js_out.update({'SigProjMatrix': sig_proj_matrix});
        js_out.update({'SigProjMatrix_p': sig_proj_matrix_p});
        js_out.update({'ProjectionKeys': sp_col_labels});
        js_out.update({'SignatureKeys': sp_row_labels});
        js_out.update({'SampleLabels': data.col_labels});
        js_out.update({'Clusters': clusters});

        fout_js.write(HtmlViewer.toJS_variable("FP_" + label, js_out));

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
    data_json = dict({
        'data': edata,
        'gene_labels': edata.row_labels,
        'sample_labels': edata.col_labels,
    });
    fout_js.write(HtmlViewer.toJS_variable("FP_ExpressionMatrix", data_json));

    fout_js.close();
    HtmlViewer.copy_html_files(dir_name);
    print();
    print("FastProject Analysis Complete")
    elapsed_time = time.time() - start_time;
    print("Elapsed Time {:.2f} seconds".format(elapsed_time));
