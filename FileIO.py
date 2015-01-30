# -*- coding: utf-8 -*-
"""
Created on Sat Jan 03 21:29:39 2015

@author: daved_000
"""
import numpy as np;

def read_matrix(filename='', delimiter = '\t'):
    """Reads data in matrix format from <filename>  Returns data_matrix, 
    row_labels, and column_labels"""
    
    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();
    
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