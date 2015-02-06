# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 15:16:22 2015

@author: David
"""

#For compatibility with Python 3.x
from __future__ import division, print_function;
if(not raw_input):
    raw_input = input;

from DimReduce import Filters;
from DimReduce import FileIO;
from DimReduce import Transforms;
from DimReduce import Signatures;
from DimReduce import Projections;
from DimReduce.Utils import ProgressBar;
import os;
import numpy as np;

housekeeping_filename = '';
def get_housekeeping_file():
    fn = housekeeping_filename;
    if(len(fn) == 0):
        while(True):
            fn = raw_input("Enter name of housekeeping file: ");
            
            if(os.path.isfile(fn)):
                break;
            else:
                print("Error : File not found");
                continue;
    
    return fn;

#%% Read expression data from file
while(True):
    filename = raw_input("Enter name of data file: ");
    
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

#%% Filtering genes
while(True):  #Loop exited with 'break', see below
    
    original = data.shape;    
    
    print();    
    print("Filtering Options Available");
    print();
    print("\t0.  Continue")
    print("\t1.  Remove housekeeping genes");
    print("\t2.  Remove inactive genes");
    print("\t3.  Filter genes for biomadility using HDT");
    print("\t4.  Save result of filtering to file");
    
    choice = raw_input("Make a selection: ");
    
    try:
        choice = int(choice);
    except ValueError:
        print("Error : Bad value.  Enter a number")
        continue;
    
    if(choice==0):
        break;
    elif(choice==1): #Housekeeping
        housekeeping_filename = get_housekeeping_file();
        (data,genes) = Filters.filter_housekeeping(data, genes, housekeeping_filename);
    elif(choice==2): #Threshold of activation
        (data, genes) = Filters.filter_genes_threshold(data, genes, 0.2);
    elif(choice==3): #HDT test
        (data, genes) = Filters.filter_genes_hdt(data, genes, 0.05);
    elif(choice==4): #Save to file
        out_file = raw_input("Enter name of file to create : ");
        FileIO.write_matrix(out_file, data, genes, cells);
        print("Data saved to " + out_file);
    else:
        print("Error : Invalid Choice\n");
        continue;

    if(choice > 0 and choice < 4):
        print();
        print("Removed ", original[0]-data.shape[0], " Genes");
        print(data.shape[0], " Genes retained");
    
#%% Probability transform?
while(True):
    
    choice = raw_input("Transform data to probability values? [y/n]: ");
    if(choice.lower()[0] == 'y'):
        #Transform based on probabilities
        housekeeping_filename = get_housekeeping_file();
        data = Transforms.probability_transform(data, original_data, original_genes, housekeeping_filename);
        
        #Prompt to save data
        while(True):
            print();
            choice = raw_input("Save result to file? [y/n]: ");
            if(choice.lower()[0] == 'y'):
                #Save data
                out_file = raw_input("Enter name of file to create : ");
                FileIO.write_matrix(out_file, data, genes, cells)
                break;
            elif(choice.lower()[0] == 'n'):
                break;
            else:
                print("Error : invalid input.  Enter Y or N");
                continue;
        break;
    elif(choice.lower()[0] == 'n'):
        break;
    else:
        print("Error : invalid input.  Enter Y or N");
        continue;

#%% Signature file
sigs = [];
USE_SIGNATURES = False;
while(True):
    print();
    print("Enter name of signature file or press enter to omit signature analysis.")
    choice = raw_input("File name : ");
    
    if(len(choice) == 0):
        break;

    elif(os.path.isfile(choice)): #Load signature
        USE_SIGNATURES = True;
        sigs = Signatures.read_signatures(choice);
        break;
    else:
        print("Error : File not found");
        continue;
         

if(USE_SIGNATURES):
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
        sig_scores[sig.name] = sig.eval_data(data, genes);
        pbar.update();
    pbar.complete();
        
    #Prompt to save data
    while(True):
        print();
        choice = raw_input("Save signature scores to file? [y/n]: ");
        if(choice.lower()[0] == 'y'):
            #Save data
            out_file = raw_input("Enter name of file to create : ");
            FileIO.write_signature_scores(out_file, sig_scores, cells);
            break;
        elif(choice.lower()[0] == 'n'):
            break;
        else:
            print("Error : invalid input.  Enter Y or N");
            continue;
        

#%% Dimensional Reduction procedures
print();
print("Projecting data into 2 dimensions");
print();

projections = Projections.generate_projections(data);

#Prompt to save data
while(True):
    print();
    choice = raw_input("Save projections to file? [y/n]: ");
    if(choice.lower()[0] == 'y'):
        #Save data
        out_file = raw_input("Enter name of file to create : ");
        Projections.write_projection_file(out_file, cells, projections);
        break;
    elif(choice.lower()[0] == 'n'):
        break;
    else:
        print("Error : invalid input.  Enter Y or N");
        continue;

#%% Output Projection Plots
print();
print('Outputting Projection Plots to File');
pp = ProgressBar(len(projections));
for proj_name in projections:
    FileIO.write_scatter_plot(
        filename = proj_name,
        x_coords = projections[proj_name][0,:],
        y_coords = projections[proj_name][1,:],
        xlabel = 'Dim 1',
        ylabel = 'Dim 2',
        title = proj_name);
    pp.update();
pp.complete();

#%% Evaluating signatures against projections
if(USE_SIGNATURES):
    N_NEIGHBORS = 5;
    sig_proj_matrix   = np.zeros((len(sigs),len(projections)));
    sig_proj_matrix_p = np.zeros((len(sigs),len(projections)));
    
    sp_row_labels = sig_scores.keys();
    sp_col_labels = projections.keys();
    
    for i, sig in enumerate(sp_row_labels):
        for j, proj in enumerate(sp_col_labels):
            print(sig, ' : ', proj);
            dissimilarity, p = Signatures.conformity_with_p(projections[proj],sig_scores[sig],N_NEIGHBORS);
            sig_proj_matrix[i,j] = np.median(dissimilarity);
            sig_proj_matrix_p[i,j] = p;
    
    
    #Output matrix of p-values for conformity scores
    FileIO.write_matrix("p_matrix.txt",sig_proj_matrix_p, sp_row_labels, sp_col_labels);
    
    
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
            filename = sig_name+'_'+proj_name,
            x_coords = projections[proj_name][0,:], 
            y_coords = projections[proj_name][1,:],
            colors   = sig_scores[sig_name],
            title    = proj_name + '\n' + sig_name + '\np = ' + '{:.3f}'.format(sig_proj_matrix_p[r,c]));
        
        pp.update();
        
    pp.complete();
