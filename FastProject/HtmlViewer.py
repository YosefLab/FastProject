"""
This module contains functions for generating html and js files so that
plots can be easily viewed.
"""
import shutil;
import os;

this_directory = os.path.dirname(os.path.abspath(__file__));
RESOURCE_DIR = this_directory + os.sep + "Viewer Resources";

def copy_html_file(destination_dir):
    shutil.copy(RESOURCE_DIR + os.sep + "viewer.html", destination_dir + os.sep + "viewer.html");



def ndarray_to_JS(np_array, format_str = "{:.3f}"):
    """
    Returns a string representation of np_array in JS list of lists syntax
    :param np_array: Array to be converted (1 or 2 dimensional)
    :return: String representation.  Not a complete JS statement

    For example, output might look like '[[1,2],[3,5]]'

    """
    if(np_array.ndim == 1):
        value = "[" + ", ".join([format_str.format(x) for x in np_array]) + "]";
    else:
        rows = list();
        for row in np_array:
            rows.append("[" + ", ".join([format_str.format(x) for x in row]) + "]");

        value = "[" + ','.join(rows) + "]";

    return value;

def dict_to_JS(py_dict, object_name, transpose_arrays = False):
    """
    Formats the dictionary of numpy ndarrays into text that can be written to a .js file for html access

    The result is a javascript object with keys equal to the projection names

    :param py_dict: dictionary with key=String, value = numpy.ndarray
    :param object_name: string representing the name of the JS object to create
    :param transpose_arrays: bool to indicate whether or not the arrays in dict should be
                             transposed before being written.
    :return: string representation to be written to JS file.
    """

    #Efficient way to do this is to build a list of strings and join it at the end
    lines = list();
    lines.append("var " + object_name + " = {};"); #Initialize empty JS object

    for key in py_dict.keys():
        line_start = object_name + "['" + key + "'] = ";
        np_array = py_dict[key];
        if(transpose_arrays):
            np_array = np_array.T;

        value = ndarray_to_JS(np_array);

        whole_line = line_start + value + ';';
        lines.append(whole_line);

    output = '\n'.join(lines) + '\n';

    return output;

def string_list_to_JS(stringlist, object_name):
    """
    Converts a list of strings into a javascript list of strings;
    :param stringlist: list of strings
    :param object_name: name of the JS object to be created
    :return: text that can be inserted into a js file
    """

    line_start = object_name + " = ";
    js_array = "[" + ", ".join(["'"+ss+"'" for ss in stringlist]) + "]";
    output = line_start + js_array + ";\n";

    return output;

def generate_js_data_file(filename, projections, sig_scores, sig_proj_matrix, sig_proj_matrix_p, projection_keys, signature_keys):
    """

    :param filename: (String) name of file to write.  Should end in .js
    :param projections: dictionary with key = (String) name of projection, value = (2xN_samples) numpy.ndarray
    :param sig_scores:  dictionary with key= (String) name of signature, value = (N_samples) numpy.ndarray
    :param sig_proj_matrix: (N_Signatures x N_Projections) numpy.ndarray
    :param sig_proj_matrix_p: (N_Signatures x N_Projections) numpy.ndarray
    :param projection_keys: List of Strings.  Projection names in order of sig_proj_matrix columns
    :param signature_keys: List of Strings.  Signatures names in order of sig_prog_matrix rows
    :return: None
    """

    with open(filename, 'w') as fout:
        fout.write(dict_to_JS(projections, "FP_Projections", transpose_arrays=True));
        fout.write(dict_to_JS(sig_scores, "FP_SigScores"))
        fout.write("FP_SigProjMatrix = " + ndarray_to_JS(sig_proj_matrix) + ";\n");
        fout.write("FP_SigProjMatrix_p = " + ndarray_to_JS(sig_proj_matrix_p) + ";\n");
        fout.write(string_list_to_JS(projection_keys,"FP_ProjectionKeys"));
        fout.write(string_list_to_JS(signature_keys,"FP_SignatureKeys"));




