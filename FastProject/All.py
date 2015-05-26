#For compatibility with Python 3.x
from __future__ import division, print_function;
import sys;
if(sys.version_info.major > 2):
    raw_input = input;

from FastProject import Filters;
from FastProject import FileIO;
from FastProject import Transforms;
from FastProject import Signatures;
from FastProject import Projections;
from FastProject.DataTypes import ExpressionData, ProbabilityData, PCData;
from FastProject.Utils import ProgressBar;
from FastProject import HtmlViewer;
import os;
import numpy as np;
import time;

from optparse import OptionParser

parser = OptionParser('usage: %prog [options] data_file');
parser.add_option("-k", "--housekeeping", metavar="FILE",
                  help="Read list of housekeeping genes from FILE.  Uses default list if not specified");
parser.add_option("-s", "--signatures", metavar="FILE",
                  help="Loads signatures from FILE.  Otherwise, signature analysis is skipped in unless in interactive mode.");
parser.add_option("-f","--filters",default="0",help="""Specifies filters to be used on genes\n\n1. Remove Housekeeping
2. Threshold (Active in at least 20% of samples)
3. Bimodal (Using Hartigans Dip Test p<0.05

e.g. -f 1, or -f 123""");
parser.add_option("-o", "--output", action="store_true", default=False, help="Outputs data after filtering and any transforms");
parser.add_option("-q", "--qc", action="store_true", default=False, help="Performs a quality check on samples, filtering samples that do not pass");
parser.add_option("-i", "--interactive", action="store_true", default=False, help="Prompts options via command line instead");


(options, args) = parser.parse_args();

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
        print("\t4.  Save result of filtering to file");
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
                if(choice > 3):
                    raise ValueError;
            except ValueError:
                print("Error:  Bad option for -f, --filters.  Exiting");
                sys.exit();

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
        edata = Filters.filter_genes_hdt(edata, 0.05);
    elif(choice==4): #Save to file
        out_file = raw_input("Enter name of file to create : ");
        FileIO.write_data(dir_name + os.sep + out_file, edata);
        print("Data saved to " + out_file);
    else:
        print("Error : Invalid Choice\n");
        continue;

    if(0 < choice < 4):
        print();
        print("Removed ", original[0]-edata.shape[0], " Genes");
        print(edata.shape[0], " Genes retained");

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

prob = ProbabilityData(prob, data);

sample_passes, sample_scores = Transforms.quality_check(params);

if(options.qc):
    prob = prob.subset_samples(sample_passes);
    data = data.subset_samples(sample_passes);

Transforms.z_normalize(data);

#Perform PCA on the data, with and without the probability xfrm
pc_data = Projections.perform_weighted_PCA(data, sample_scores**-1);
pc_data = PCData(pc_data, data);
pc_data = Projections.filter_PCA(pc_data, sample_scores);


pc_prob = Projections.perform_weightd_PCA(prob, sample_scores**-1);
pc_prob = PCData(pc_prob, prob);
pc_prob = Projections.filter_PCA(pc_prob, sample_scores);

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

fout_js = open(dir_name + os.sep + "FP_data.js", 'w');

for label, data in zip(data_labels, data_matrices):
    print("Evaluating: ", label);
    #%% Dimensional Reduction procedures
    print();
    print("Projecting data into 2 dimensions");

    projections = Projections.generate_projections(data);

    #Save data
    out_file = label + '_Projections.txt'
    Projections.write_projection_file(dir_name + os.sep + out_file, cells, projections);

    #Evaluate Signatures
    print();
    print("Evaluating signature scores on samples...");

    sig_scores = dict();

    pbar = ProgressBar(len(sigs));
    for sig in sigs:
        try:
            sig_scores[sig.name] = data.eval_signature(sig, fn_prob);
        except ValueError:  #Only thrown when the signature has no genes in the data
            pass #Just discard the signature then
        pbar.update();
    pbar.complete();

    #Prompt to save data
    out_file = label + '_SigScores.txt';
    FileIO.write_signature_scores(dir_name + os.sep + out_file, sig_scores, cells);



    #%% Evaluating signatures against projections
    N_NEIGHBORS = 5;
    sig_proj_matrix   = np.zeros((len(sig_scores),len(projections)));
    sig_proj_matrix_p = np.zeros((len(sig_scores),len(projections)));

    sp_row_labels = sig_scores.keys();
    sp_col_labels = projections.keys();

    print();
    print("Evaluating Signatures against Projections");
    pp = ProgressBar(len(sp_row_labels) * len(sp_col_labels));
    for i, sig in enumerate(sp_row_labels):
        for j, proj in enumerate(sp_col_labels):
            dissimilarity, p = Signatures.conformity_with_p(projections[proj],sig_scores[sig],N_NEIGHBORS);
            sig_proj_matrix[i,j] = np.median(dissimilarity);
            sig_proj_matrix_p[i,j] = p;
            pp.update();

    pp.complete();

    #Output matrix of p-values for conformity scores
    FileIO.write_matrix(dir_name + os.sep + label + "_PMatrix.txt",sig_proj_matrix_p, sp_row_labels, sp_col_labels);

    #Wrap data into an object
    js_out = dict();
    js_out.update({'Projections': projections});
    js_out.update({'SigScores': sig_scores});
    js_out.update({'SigProjMatrix': sig_proj_matrix});
    js_out.update({'SigProjMatrix_p': sig_proj_matrix_p});
    js_out.update({'ProjectionKeys': sp_col_labels});
    js_out.update({'SignatureKeys': sp_row_labels});

    fout_js.write(HtmlViewer.toJS_variable("FP_" + label, js_out));

fout_js.close();
HtmlViewer.copy_html_files(dir_name);
print("FastProject Analysis Complete")
