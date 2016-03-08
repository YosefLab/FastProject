# -*- coding: utf-8 -*-
"""Code to handle reading and writing various files

Methods that interact with files are found here.
Except for cases when it made more sense to put them
in another module
  For example: reading Signatures is handled in the Signatures module

"""
from __future__ import absolute_import, print_function, division;
import os;
import numpy as np;
import shutil;
from . import HtmlViewer;
import pandas as pd;


def make_dirs(root_directory):
    """Creates relevant directories that will be needed for output files.

    Parameters
    ----------
    root_directory : str
        Root directory for FastProject output

    """

    try:
        os.makedirs(root_directory);
    except OSError:
        pass;

    html_sub_directory = os.path.join(root_directory, HtmlViewer.OUTPUT_RESOURCE_DIRECTORY);
    try:
        os.makedirs(html_sub_directory);
    except OSError:
        pass;


def read_matrix(filename='', delimiter = '\t'):
    """Reads gene expression matrix from `filename`

    Parameters
    ----------
    filename : str
        Name (with path) of file to read from
    delimiter : str, optional
        Delimiter used in file to separate columns
        Default: '\t' (tab-delimited)

    Returns
    -------
    data_matrix : numpy.ndarray
        Contains expression levels from `filename`
        Size is NUM_GENES x NUM_SAMPLES
    row_labels: [str]
        Label for each row - names for each gene
    col_labels: [str]
        Label for each column - names for each sample
    """

    if(filename == ''):
        from Tkinter import Tk
        from tkFileDialog import askopenfilename
        Tk().withdraw();
        filename = askopenfilename();
    
    print("Loading data from " + filename + " ...");
    #load matrix data
    #First read in list of cell names
    
    with open(filename,'rU') as ff:
        col_labels = ff.readline().strip().split(delimiter);
        col_labels = [col_label.strip().strip("'\"") for col_label in col_labels];

        #col_labels in first line may or may not be shifted by one to line up with
        #data in following lines.  Determine # of datapoints from second line
        line2 = ff.readline().strip().split(delimiter);
        if(len(line2) == len(col_labels)):
            col_labels = col_labels[1:];


    #read all columns, start at 1 to skip first column of row labels
    cols_to_read = np.arange(1,len(col_labels)+1);
    
    data = np.loadtxt(filename,delimiter=delimiter,skiprows=1,usecols=cols_to_read);
    
    if(data.std() > 10 or data.mean() > 10):
        print("Log-Transforming data...");
        data = np.log(data + 1);
    
    #read in gene names
    with open(filename,'rU') as ff:
        xx = ff.readline();

        row_labels = list();
        for line in ff:
          entries = line.split(delimiter);
          row_labels.append(entries[0].strip("'\"").upper());
    

    #Test for unique row and column labels
    if(len(set(row_labels)) != len(row_labels)):
        #Find duplicates
        row_labels_copy = [r for r in row_labels];
        for row_label in set(row_labels_copy):
            row_labels_copy.remove(row_label);

        #Warning for duplicate genes
        print("WARNING: Row labels (gene identifiers) are not unique.\n" + ", ".join(row_labels_copy));
        print("Removing these genes...");

        row_labels_copy = set(row_labels_copy);
        to_remove_i = [i for i,x in enumerate(row_labels) if x in row_labels_copy];
        to_keep_i = np.setxor1d(np.arange(len(row_labels)), to_remove_i);
        data = data[to_keep_i,:];
        row_labels = [row_labels[i] for i in to_keep_i];

    if(len(set(col_labels)) != len(col_labels)):
        #Find duplicates
        print("WARNING: Column labels (sample names) are not unique.\n Made unique by adding a suffix.")
        _uniquify(col_labels);

    return (data, row_labels, col_labels);

def _uniquify(labels):
    """
    Takes a list of strings, and makes it unique by appending _1, _2, etc to duplicates
    In-Place

    :param labels: Input list of strings
    :return: None, modification is in-place
    """

    count_dict = dict();
    for label in labels:
        if(label not in count_dict):
            count_dict.update({label: 0});

        count_dict[label] = count_dict[label] + 1;

    for i in range(len(labels))[::-1]:
        if(count_dict[label] > 1):
            labels[i] = labels[i] + "_" + str(count_dict[label]);
            count_dict[label] = count_dict[label] - 1;


def write_matrix(filename, data, row_labels, col_labels):
    ff = open(filename, 'w');
    
    ff.write('\t' + '\t'.join(col_labels) + '\n'); #Extra tab so col labels line up with data
    for i in range(data.shape[0]):
        ff.write(row_labels[i] + '\t');
        ff.write('\t'.join(["{:.5f}".format(num).rstrip('0').rstrip('.') for num in data[i,:]]));
        ff.write('\n');
    
    ff.close();  

def write_data(filename, data):
    write_matrix(filename, data, data.row_labels, data.col_labels); 
    
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
    
    row_labels = [str(i) for i in range(num_rows)];
    col_labels = [str(i) for i in range(num_cols)];
    
    return (data, row_labels, col_labels);
    
def write_signature_scores(filename, sig_scores_dict, col_labels):
    with open(filename, 'w') as ff:
        #Write column headers to first row
        ff.write('\t' + '\t'.join(col_labels) + '\n');

        #Iterate over signatures and write scores
        for name in sig_scores_dict:
            sig_scores = sig_scores_dict[name]
            if(sig_scores.isFactor):
                row = [name] + sig_scores.scores;
            else:
                row = [name] + ["{:.5f}".format(num).rstrip('0').rstrip('.') for num in sig_scores.scores];

            ff.write('\t'.join(row) + '\n');

def write_projection_file(filename, sample_labels, projections):
    """
    Outputs the coordinates for each projection to a file.

    Parameters
    ----------
    filename : String
        File to write
    sample_labels : list(String) len=Num_Samples
        Labels for each sample
    projections : dict(string, (2 x Num_Samples) numpy.ndarray)
        dictionary mapping the projection type (e.g. "tSNE") to an array containing
        the two-dimensional coordinates for each sample in the projection.

    """

    ff = open(filename,'w');

    for proj in projections:
        coordinates = projections[proj];
        for i in range(coordinates.shape[1]):
            ff.write(proj + '\t');
            ff.write(sample_labels[i] + '\t');
            ff.write('{:.5f}'.format(coordinates[0,i]) + '\t');
            ff.write('{:.5f}'.format(coordinates[1,i]));
            ff.write('\n');

    ff.close();

def write_cluster_file(filename, sample_labels, clusters):
    """
    Outputs cluster assignments to a file

    :param filename: String
        File to be written
    :param sample_labels: list(String)
        Labels for each sample
    :param clusters: Dict(projection_name (String) -> Projection Clusters)
        Projection Clusters: Dict(cluster_name (String) -> Cluster Assignments)
        Cluster Assignments: numpy.ndarray (len = num_samples)
            Cluster assignments for each sample

    """

    with open(filename, 'w') as fout:
        fout.write('\t'.join(sample_labels) + '\n');

        for proj in clusters:
            proj_clusters = clusters[proj];
            for cluster in proj_clusters:
                assignments = proj_clusters[cluster];
                fout.write(proj + " " + cluster + "\t");
                fout.write('\t'.join([str(x) for x in assignments]));
                fout.write('\n');

def write_filter_file(filename, filter):
    """
    Outputs cluster assignments to a file

    :param filename: String
        File to be written
    :param filter: set(String)
        Set of gene names included in the filter

    """

    with open(filename, 'w') as fout:
        fout.write('\t'.join(filter) + '\n');

def write_models(directory, Models):
    """
    Writes the resulting Models to the output directory tree

    :param directory: Root directory for FP output
    :param Models: Nested dictionary structure containing all output data
    :return:
    """

    for name in Models:
        model = Models[name];
        model_dir = os.path.join(directory, name);

        #Sig Scores
        out_file = 'SignatureScores.txt';
        write_signature_scores(os.path.join(model_dir, out_file), model["signatureScores"], model["sampleLabels"]);

        for projData in model["projectionData"]:
            if(projData["filter"] == "No_Filter"):
                filter_dir = os.path.join(model_dir, "No_Filter");
            else:
                filter_dir = os.path.join(model_dir, projData["filter"] + "_Filter");

            if(projData["pca"]):
                filter_dir = filter_dir + "_PCA";

            try:
                os.makedirs(filter_dir);
            except OSError:
                pass;

            #Save Projections
            write_projection_file(os.path.join(filter_dir, 'Projections.txt'),
                                 model["sampleLabels"],
                                  projData["projections"]);

            #Save Clusters
            write_cluster_file(os.path.join(filter_dir, 'Clusters.txt'),
                               model["sampleLabels"],
                               projData["clusters"])

            #Output matrix of p-values for conformity scores
            write_matrix(os.path.join(filter_dir, "DissimilarityMatrix.txt"),
                         projData["sigProjMatrix"],
                         projData["signatureKeys"],
                         projData["projectionKeys"]);

            write_matrix(os.path.join(filter_dir, "PMatrix.txt"),
                        projData["sigProjMatrix_p"],
                        projData["signatureKeys"],
                        projData["projectionKeys"]);

            #Output genes used in filter
            write_filter_file(os.path.join(filter_dir, 'ProjectedGenes.txt'),
                              projData["genes"]);

            #If PCA, output the Loadings
            if(projData["pca"]):
               write_matrix(os.path.join(filter_dir, "PCA_Loadings.txt"),
                            projData["loadings"],
                            projData["genes"],
                            ["PC "+str(x+1) for x in range(projData["loadings"].shape[1])]);

def write_weights(directory, data):
    """
    Writes the matrix of weights to file.

    :param directory: Folder to which file is written
    :param data: ExpressionData or ProbabilityData object
    :return: None
    """

    weight_frame = pd.DataFrame(data.weights, index=data.row_labels, columns=data.col_labels);
    out_file = os.path.join(directory, "weights.txt");
    weight_frame.to_csv(out_file, sep="\t");



def write_qc_file(directory, sample_passes, sample_scores, sample_labels):
    """
    Outputs a file with Quality Report information for the samples

    :param directory: String
        Output directory to use
    :param sample_passes: np.ndarray, dtype=bool, size=NUM_SAMPLES
        Whether or not the sample passes the quality check
    :param sample_scores: np.ndarray, dtype=float, size=NUM_SAMPLES
        A samples quality score
    :param sample_labels: list(String)
        Labels for each sample
    :return: None
    """

    # Compute MAD
    med = np.median(sample_scores);
    mad = np.median(np.abs(sample_scores - med));
    cutoff = med - 1.6 * mad;

    # Build data object
    values = HtmlViewer.toJS_variable("values", sample_scores);

    # Build html Table
    ii = np.argsort(sample_scores);

    sample_scores = sample_scores[ii];
    sample_labels = [sample_labels[i] for i in ii];

    out_table = [];
    out_table.append("<table>");
    out_table.append("<tr>");
    out_table.append("<th>Sample</th>");
    out_table.append("<th>Score</th>");
    out_table.append("</tr>");

    for score, label in zip(sample_scores, sample_labels):
        if(score < cutoff):
            out_table.append('<tr class="failed_row">');
        else:
            out_table.append("<tr>");
        out_table.append("<td>" + label + "</td>");
        out_table.append("<td>" + "{:.3f}".format(score) + "</td>");
        out_table.append("</tr>");

    out_table.append("</table>");


    #Output Html File
    out_html = os.path.join(directory,"QC_Report.html");

    shutil.copy(HtmlViewer.RESOURCE_DIR + os.sep + "QC_Report.html",out_html);

    with open(out_html, "r") as fin:
        contents = fin.read();

    contents = contents.replace("<% table %>", "\n".join(out_table));
    contents = contents.replace("<% values %>", values);
    contents = contents.replace("<% MAD %>", HtmlViewer.toJS_variable("MAD", cutoff));

    with open(out_html, "w") as fout:
        fout.write(contents);





