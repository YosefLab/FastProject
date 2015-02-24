# -*- coding: utf-8 -*-
"""
Created on Sat Jan 03 21:29:39 2015

@author: daved_000
"""
import numpy as np;
import matplotlib.pyplot as plt;

def read_matrix(filename='', delimiter = '\t'):
    """Reads data in matrix format from <filename>  Returns data_matrix, 
    row_labels, and column_labels"""
    
    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();
    
    print("Loading data from " + filename + " ...");
    #load matrix data
    #First read in list of cell names
    
    ff = open(filename,'r');
    col_labels = ff.readline().strip().split(delimiter);
    col_labels = [col_label.strip() for col_label in col_labels];
    
    #col_labels in first line may or may not be shifted by one to line up with 
    #data in following lines.  Determine # of datapoints from second line
    line2 = ff.readline().strip().split(delimiter);
    if(len(line2) == len(col_labels)):
        col_labels = col_labels[1:];
    
    ff.close();
    
    #read all columns, start at 1 to skip first column of row labels
    cols_to_read = np.arange(1,len(col_labels)+1);
    
    data = np.loadtxt(filename,delimiter=delimiter,skiprows=1,usecols=cols_to_read);
    
    if(data.std() > 10 or data.mean() > 10):
        print "Log-Transforming data..."
        data = np.log(data + 1);
    
    #read in gene names
    ff = open(filename,'r');
    xx = ff.readline();
    
    row_labels = list();
    for line in ff:
      entries = line.split(delimiter);
      row_labels.append(entries[0].split('.')[0]);  #remove anything after the decimal
      #or
      #row_labels.append(entries[0]);
    
    ff.close();
    
    return (data, row_labels, col_labels);

def write_matrix(filename, data, row_labels, col_labels):
    ff = open(filename, 'w');
    
    ff.write('\t' + '\t'.join(col_labels) + '\n'); #Extra tab so col labels line up with data
    for i in range(data.shape[0]):
        ff.write(row_labels[i] + '\t');
        ff.write('\t'.join(["{:.5f}".format(num).rstrip('0').rstrip('.') for num in data[i,:]]));
        ff.write('\n');
    
    ff.close();  
    
def read_matrix_nolabels(filename = '', delimiter = '\t'):
    
    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();

    #read data    
    data = np.loadtxt(filename, delimiter=delimiter);
    
    #generate labels
    num_rows = data.shape[0];
    num_cols = data.shape[1];
    
    row_labels = [str(i) for i in xrange(num_rows)];
    col_labels = [str(i) for i in xrange(num_cols)];
    
    return (data, row_labels, col_labels);
    
def write_signature_scores(filename, sig_scores, col_labels):
    row_labels = sig_scores.keys();
    data_matrix = np.zeros((len(row_labels), len(col_labels)));
    for i, row in enumerate(row_labels):
        data_matrix[i,:] = sig_scores[row];
    
    write_matrix(filename, data_matrix, row_labels, col_labels);

def write_scatter_plot(filename, x_coords, y_coords, colors=[], xlabel='', ylabel='', title=''):
    #Add .png extension if filename lacks it
    if(not filename.endswith(".png")):
        filename = filename + ".png";
        
    plt.ioff();
    ff = plt.figure();
    if(len(colors) == len(x_coords)):
        plt.scatter(x_coords, y_coords, c=colors);
        plt.colorbar();
    else:
        plt.scatter(x_coords, y_coords);
    
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);
    plt.title(title);
    
    plt.savefig(filename, bbox_inches='tight');

    plt.close(ff);

