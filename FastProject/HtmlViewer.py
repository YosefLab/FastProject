# -*- coding: utf-8 -*-
"""Contains functions for generating result html and js files 

"""
import shutil;
import os;
import sys
import numpy as np;
from .Global import get_viewer_resource_dir;

if sys.version_info.major == 2:
    STR_TYPES = (str,unicode)
else:
    STR_TYPES = (str,)


RESOURCE_DIR = get_viewer_resource_dir();
OUTPUT_RESOURCE_DIRECTORY = "_resources";

def get_output_js_handle(destination_dir):
    fout_js = open(os.path.join(destination_dir,OUTPUT_RESOURCE_DIRECTORY,"FP_data.jsdata"), 'w');
    return fout_js;

def copy_html_files(destination_dir):
    shutil.copy(RESOURCE_DIR + os.sep + "Results.html", destination_dir + os.sep + "Results.html");

    html_sub_directory = os.path.join(destination_dir,OUTPUT_RESOURCE_DIRECTORY)

    shutil.copy(RESOURCE_DIR + os.sep + "jquery-2.1.4.min.js", html_sub_directory + os.sep + "jquery-2.1.4.min.js");
    shutil.copy(RESOURCE_DIR + os.sep + "bootstrap.css", html_sub_directory + os.sep + "bootstrap.css");
    shutil.copy(RESOURCE_DIR + os.sep + "bootstrap.min.js", html_sub_directory + os.sep + "bootstrap.min.js");
    shutil.copy(RESOURCE_DIR + os.sep + "d3.min.js", html_sub_directory + os.sep + "d3.min.js");
    shutil.copy(RESOURCE_DIR + os.sep + "d3tip.js", html_sub_directory + os.sep + "d3tip.js"); #Extension For tooltips
    shutil.copy(RESOURCE_DIR + os.sep + "ColorScatter.js", html_sub_directory + os.sep + "ColorScatter.js");
    shutil.copy(RESOURCE_DIR + os.sep + "HeatMap.js", html_sub_directory + os.sep + "HeatMap.js");
    shutil.copy(RESOURCE_DIR + os.sep + "Clustered_HeatMap.js", html_sub_directory + os.sep + "Clustered_HeatMap.js");
    shutil.copy(RESOURCE_DIR + os.sep + "jszip.min.js", html_sub_directory + os.sep + "jszip.min.js");
    shutil.copy(RESOURCE_DIR + os.sep + "styles.css", html_sub_directory + os.sep + "styles.css");
    shutil.copy(RESOURCE_DIR + os.sep + "Utilities.js", html_sub_directory + os.sep + "Utilities.js");

def toJS_variable(variable_name, obj):
    return "var " + variable_name + " = " + toJS(obj) + ";\n";

def toJS(obj, level=0):
    """
    Convert a python object into a javascript format.
    Only certain object types are supported:
        List >> Javascript Array
        Numpy.ndarray >> Javascript Array
        Dictionary >> Javascript Object if keys are strings
    :param obj: object to be converted
    :param level: determines indentation level
    :return: string representation of the object in javascript
    """

    try:
        obj.to_JSON;
    except AttributeError:
        pass;
    else:
        return obj.to_JSON();

    if(isinstance(obj, STR_TYPES)):
        return '"' + obj + '"';
    if(type(obj) is int or type(obj) is float or isinstance(obj, np.number)):
        return str(obj);
    if(issubclass(type(obj),np.ndarray) and obj.ndim < 3):
        return ndarray_to_JS(obj);
    if(type(obj) is list or type(obj) is set):
        return '[' + ','.join([toJS(x) for x in obj]) + ']';
    if(type(obj) is bool):
        if(obj): return 'true';
        else: return 'false';
    if(type(obj) is dict):
        pairs = list();
        indentation = "".join(["\t"]*(level+1));
        for key in obj:
            if(isinstance(key, STR_TYPES) or type(key) is int or type(key) is float):
                pairs.append("\n"+indentation+toJS(key, level+1)+':'+toJS(obj[key], level+1));
            else:
                raise ValueError("Non-compatible Value Encountered for Object Key:", type(key));
        return '{' + ','.join(pairs) + '}';

    raise ValueError("Non-convertible Value Encountered:", type(obj));

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

    for key in py_dict:
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




