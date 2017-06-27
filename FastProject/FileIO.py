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
from .Global import FP_Output, get_housekeeping_dir
import pandas as pd

try:
    from pandas.errors import ParserError
except ImportError:
    from pandas.parser import CParserError as ParserError


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

        try:
            os.makedirs(model_dir);
        except OSError:
            pass;

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
            write_matrix(os.path.join(filter_dir, "ConsistencyMatrix.txt"),
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
    :param data: ExpressionData
    :return: None
    """

    weight_frame = pd.DataFrame(data.weights, index=data.row_labels, columns=data.col_labels);
    out_file = os.path.join(directory, "weights.txt");
    weight_frame.to_csv(out_file, sep="\t");


def write_qc_file(directory, qc_info):
    """ Outputs a file with Quality Report information for the samples

    Parameters
    ----------
    directory : str
        Output directory to use
    qc_info : pandas.DataFrame
        Columns are "Scores" and "Passes"
        "Score": quality score
        "Passes": bool, whether or not the sample passed
        Index: sample label

    """

    sample_scores = qc_info["Score"].values

    # Compute MAD
    med = np.median(sample_scores);
    mad = np.median(np.abs(sample_scores - med));
    cutoff = med - 1.6 * mad;

    # Build data object
    values = HtmlViewer.toJS_variable("values", sample_scores);

    # Build html Table
    qc_info = qc_info.sort_values("Score");

    out_table = [];
    out_table.append("<table>");
    out_table.append("<tr>");
    out_table.append("<th>Sample</th>");
    out_table.append("<th>Score</th>");
    out_table.append("</tr>");

    for label in qc_info.index:
        score = qc_info.loc[label, "Score"];
        if(score < cutoff):
            out_table.append('<tr class="failed_row">');
        else:
            out_table.append("<tr>");
        out_table.append("<td>" + label + "</td>");
        out_table.append("<td>" + "{:.3f}".format(score) + "</td>");
        out_table.append("</tr>");

    out_table.append("</table>");

    # Output Html File
    out_html = os.path.join(directory, "QC_Report.html");

    shutil.copy(HtmlViewer.RESOURCE_DIR + os.sep + "QC_Report.html", out_html);

    with open(out_html, "r") as fin:
        contents = fin.read();

    contents = contents.replace("<% table %>", "\n".join(out_table));
    contents = contents.replace("<% values %>", values);
    contents = contents.replace("<% MAD %>", HtmlViewer.toJS_variable("MAD", cutoff));

    with open(out_html, "w") as fout:
        fout.write(contents);


def load_input_projections(projection_files, sample_names):
    """
    Loads input projection coordinates into module variable _input_projections

    :param projection_files: List of file locations (str) that contain projection coordinates
    :param sample_names: List of sample names.  Used to validate input

    Returns
    -------
    input_projections : dict or None
        Keys are of type str, representing projection names
        Values are of type 2xN pandas.DataFrame with column names matching
        sample names in `expressionMatrix`
    """
    input_projections = {};

    sample_names_norm = [x.upper() for x in sample_names];

    for projection_file in projection_files:
        if(not os.path.isfile(projection_file)):
            raise ValueError("Option Error: projection file " + projection_file + " not found.");

        try:
            loaded_coords = pd.read_table(projection_file, header=None);
        except ParserError:
            # if Header row has 2 entries and rest has 3, get this error
            # Try to load data anyway and rework it back into the correct initial format
            loaded_coords = pd.read_table(projection_file);
            loaded_coords = pd.DataFrame({0: loaded_coords.index,
                                          1: loaded_coords.iloc[:, 0].values,
                                          2: loaded_coords.iloc[:, 1].values},
                                         index=np.arange(loaded_coords.shape[0]));

        # Verify the file
        # Check for 3 columns
        if(loaded_coords.shape[1] != 3):
            raise ValueError("Format Error: projection file " + projection_file +
                             " should have 3 columns (sample name, x coordinate, y coordinate).");

        # Fix in case header was included with first cell empty
        if(pd.isnull(loaded_coords.iloc[0, 0])):
            loaded_coords.iloc[0, 0] = "_____";  # Dummy value, won't match gene

        # Normalize names to upper-case temporarily
        loaded_coords[0] = [x.upper() for x in loaded_coords[0]];

        # Did they accidentally include a heading?
        if(loaded_coords.iloc[0, 0] not in sample_names_norm):
            loaded_coords = loaded_coords.drop(loaded_coords.index[0], axis='index');  # drop first row

        # Check for duplicate names
        if(loaded_coords[0].duplicated().any()):
            duplicated_names = loaded_coords.loc[loaded_coords[0].duplicated(), 0].values;
            raise ValueError("Format Error: projection file " + projection_file +
                             " contains duplicate sample names:\n" + "\n".join(duplicated_names));

        # Set as the index
        loaded_coords = loaded_coords.set_index(0);

        # Check that all data matrix sample names are in the projection file
        not_matched = np.array([x not in loaded_coords.index for x in sample_names_norm]);
        if(not_matched.any()):
            missing_samples = [sample_names[i] for i in np.nonzero(not_matched)[0]];
            raise ValueError("Format Error: projection file " + projection_file +
                             " missing entries for:\n" + "\n".join(missing_samples));

        # If we got here, it's safe to re-order the file
        loaded_coords = loaded_coords.loc[sample_names_norm, :];
        loaded_coords.index = sample_names;  # Restore original case
        loaded_coords["X"] = loaded_coords[1];
        loaded_coords["Y"] = loaded_coords[2];
        loaded_coords = loaded_coords.drop(1, axis="columns").drop(2, axis="columns");

        # Check for null values and throw error
        if(loaded_coords.isnull().any().any()):
            bad_rows = loaded_coords.index[loaded_coords.isnull().any(axis=1)];
            raise ValueError("Format Error: projection file " + projection_file +
                             " contains NaNs for samples:\n" + "\n".join(bad_rows));

        proj_name = os.path.splitext(projection_file)[0];
        input_projections.update({proj_name: loaded_coords});

    return input_projections;


def load_input_weights(weight_file, genes, cells):
    """
    Used to pre-load an input weight matrix

    Parameters
    ----------
    weight_file : str
        File to load weights from
    genes : list of str
        Names of genes in data matrix (for validation)
    cells : list of str
        Names of cells in data matrix (for validation)

    Returns
    -------
    input_weights : pandas.DataFrame
        Shape is len(`genes`) x len(`cells`)
        Values are floats from 0.0 to 1.0

    """

    if(not os.path.isfile(weight_file)):
        raise ValueError("Option Error: weights file " + weight_file + " not found.");

    weights = pd.read_table(weight_file, header=0, index_col=0);
    weights.index = [x.upper() for x in weights.index];  # Normalize gene names to upper-case

    # Correct any NaNs
    if(pd.isnull(weights).any().any()):
        FP_Output("NaN's found in input weight matrix.  Replacing with zeros...");
        weights[pd.isnull(weights)] = 0.0;

    # Check that we have a gene/sample pair for everyone
    missing_genes = np.array([x not in weights.index for x in genes]);
    missing_samples = np.array([x not in weights.columns for x in cells]);

    if(missing_genes.any()):
        missing_gene_labels = [x for i, x in enumerate(genes) if missing_genes[i]];
        raise ValueError("Format Error: weights file " + weight_file +
                         " missing entries for:\n" + "\n".join(missing_gene_labels));

    if(missing_samples.any()):
        missing_sample_labels = [x for i, x in enumerate(cells) if missing_samples[i]];
        raise ValueError("Format Error: weights file " + weight_file +
                         " missing entries for:\n" + "\n".join(missing_sample_labels));

    return weights;


def saveResultstoDisk(models, signatures, qc_info, dir_name):
    """
    Save the results of a FastProject analysis to a directory structure in `dir_name`
    """

    fout_js = HtmlViewer.get_output_js_handle(dir_name);

    # Write signatures to file
    # Assemble signatures into an object, then convert to JSON variable and write
    sig_dict = {};
    for sig in signatures:
        sig_genes = list(sig.sig_dict.keys());
        sig_values = list(sig.sig_dict.values());
        sort_i = np.array(sig_values).argsort()[::-1];  # Put positive signatures first
        sig_genes = [sig_genes[i] for i in sort_i];
        sig_values = [sig_values[i] for i in sort_i];
        sig_dict.update({sig.name: {'Genes': sig_genes, 'Signs': sig_values}});
    fout_js.write(HtmlViewer.toJS_variable("FP_Signatures", sig_dict));

    edata = models["Expression"]["Data"];
    data_json = dict({
        'data': edata,
        'gene_labels': edata.row_labels,
        'sample_labels': edata.col_labels,
    });
    fout_js.write(HtmlViewer.toJS_variable("FP_ExpressionMatrix", data_json));

    write_models(dir_name, models);
    write_weights(dir_name, models["Expression"]["Data"]);

    # Remove model["Data"] since only the expression data is written for JS
    #  and it's already written above in FP_ExpressionMatrix
    for name, model in models.items():
        model.pop("Data");

    fout_js.write(HtmlViewer.toJS_variable("FP_Models", models));

    fout_js.close();
    HtmlViewer.copy_html_files(dir_name);

    write_qc_file(dir_name, qc_info);


def load_housekeeping_genes(filename=""):
    """
    Loads housekeeping genes from `filename`

    If `filename` is blank, then loads all the default list

    Returns: list of str
        Housekeeping gene identifiers
    """

    housekeeping_files = list()

    if(filename != ""):  # If file specified, use that file
        housekeeping_files.append(filename)
    else:  # Otherwise, use all the files in housekeeping directory
        housekeeping_dir = get_housekeeping_dir()
        files = os.listdir(housekeeping_dir)
        for ff in files:
            housekeeping_files.append(os.path.join(housekeeping_dir, ff))

    housekeeping_genes = list()
    for hkf in housekeeping_files:
        with open(hkf, 'rU') as fin:
            for line in fin.readlines():
                housekeeping_genes.append(line.strip().lower())

    return housekeeping_genes
