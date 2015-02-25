# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 15:16:22 2015

@author: David
"""

#For compatibility with Python 3.x
from __future__ import division, print_function;
import sys;    
if(sys.version_info.major > 2):
    raw_input = input;

from . import Filters;
from . import FileIO;
from . import Transforms;
from . import Signatures;
from . import Projections;
from .Utils import ProgressBar;
import os;
import numpy as np;
import time;

from optparse import OptionParser

PCA_TRANSFORM = False;

HAS_NUMBA = False;
try:
    import numba
    HAS_NUMBA = True;
except ImportError:
    HAS_NUMBA = False;
    

parser = OptionParser('usage: %prog [options] data_file');
parser.add_option("-k", "--housekeeping", metavar="FILE", 
                  help="Read list of housekeeping genes from FILE.  Uses default list if not specified");
parser.add_option("-s", "--signatures", metavar="FILE",
                  help="Loads signatures from FILE.  Otherwise, signature analysis is skipped in unless in interactive mode.");
parser.add_option("-f","--filters",default="0",help="""Specifies filters to be used on genes\n\n1. Remove Housekeeping
2. Threshold (Active in at least 20% of samples)
3. Bimodal (Using Hartigans Dip Test p<0.05

e.g. -f 1, or -f 123""");
parser.add_option("-p","--probability", action="store_true", default=False, help="Projects using probability of expression rather than log expression level");
parser.add_option("-c","--pca", action="store_true", default=False, help="Reduces to a smaller number of principal components before projection");
parser.add_option("-o", "--output", action="store_true", default=False, help="Outputs data after filtering and any transforms");
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
        else:
            print();
            print("Argument Error:  data_file not specified.\nExiting...");
            sys.exit();
    
    if(os.path.isfile(filename)):
        (data, genes, cells) = FileIO.read_matrix(filename);
        break;
    else:
        print("Error : File not found");
        continue

print("Imported ", data.shape[0], " genes across ", data.shape[1], " samples");

#Hold on to originals so we don't lose data after filtering in case it's needed later
original_data = data;
original_genes = genes;
original_cells = cells;

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
    
    original = data.shape;    
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
        (data,genes) = Filters.filter_housekeeping(data, genes, housekeeping_filename);
    elif(choice==2): #Threshold of activation
        print("Removing genes inactive in > 80% samples...");
        (data, genes) = Filters.filter_genes_threshold(data, genes, 0.2);
    elif(choice==3): #HDT test
        if(HAS_NUMBA):
            print("Removing genes with unimodal distribution across samples using Hartigans DT...");
            (data, genes) = Filters.filter_genes_hdt(data, genes, 0.05);
        else:
            print("Unavailable: This filtering method requires the python package 'numba'");
    elif(choice==4): #Save to file
        out_file = raw_input("Enter name of file to create : ");
        FileIO.write_matrix(dir_name + os.sep + out_file, data, genes, cells);
        print("Data saved to " + out_file);
    else:
        print("Error : Invalid Choice\n");
        continue;

    if(choice > 0 and choice < 4):
        print();
        print("Removed ", original[0]-data.shape[0], " Genes");
        print(data.shape[0], " Genes retained");
    
#%% Probability transform
housekeeping_filename = get_housekeeping_file();
prob, fn_prob = Transforms.probability_transform(data, original_data, original_genes, housekeeping_filename);    
        
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
        choice = raw_input("Reduce data into top 30 pincipal components before projection? [y/n]: ");
    else:
        choice = 'y' if options.pca else 'n';
    
    if(choice.lower()[0] == 'y'):
        #Transform into top N principal components
        PCA_TRANSFORM = True;
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
            
        FileIO.write_matrix(dir_name + os.sep + out_file, data, genes, cells)
        break;
    elif(choice.lower()[0] == 'n'):
        break;
    else:
        print("Error : invalid input.  Enter Y or N");
        continue;

#%% Dimensional Reduction procedures
print();
print("Projecting data into 2 dimensions");

if(PCA_TRANSFORM):
    pc_data, pc_rows = Projections.perform_PCA(data, genes, 30);
    projections = Projections.generate_projections(pc_data);
else:
    projections = Projections.generate_projections(data);

#Save data
out_file = 'projections.txt'
Projections.write_projection_file(dir_name + os.sep + out_file, cells, projections);
        

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

if(USE_SIGNATURES):
    print();
    print("Evaluating signature scores on samples...");
    
    sig_scores = dict();
    
    pbar = ProgressBar(len(sigs));
    for sig in sigs:
        sig_scores[sig.name] = sig.eval_data(data, genes, fn_prob);
        pbar.update();
    pbar.complete();
        
    #Prompt to save data
    out_file = 'sig_scores.txt';
    FileIO.write_signature_scores(dir_name + os.sep + out_file, sig_scores, cells);
      
        

#%% Evaluating signatures against projections
if(USE_SIGNATURES):
    N_NEIGHBORS = 5;
    sig_proj_matrix   = np.zeros((len(sigs),len(projections)));
    sig_proj_matrix_p = np.zeros((len(sigs),len(projections)));
    
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
        proj_name = sp_col_labels[c];
        
        FileIO.write_scatter_plot(
            filename = dir_name + os.sep + sig_name+'_'+proj_name,
            x_coords = projections[proj_name][0,:], 
            y_coords = projections[proj_name][1,:],
            colors   = sig_scores[sig_name],
            title    = proj_name + '\n' + sig_name + '\np = ' + '{:.3f}'.format(sig_proj_matrix_p[r,c]));
        
        pp.update();
        
    pp.complete();

