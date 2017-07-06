"""
A module for functions which convert
FastProject datatypes to/from JSON
"""

from .DataTypes import ExpressionData
from .Signatures import Signature
from .SigScoreMethods import SignatureScores
import json
import numpy as np


def expressionData_to_JSON(expressionData):
    """Converts a Datatypes.ExpressionData into its
    JSON representation

    Only uses the data values, row labels, and column labels.

    An expressionData in JSON form is specified as:

        {'col_labels': ['sample1', 'sample2', 'sample3', ...],
         'row_labels': ['feature1', 'feature2', 'feature3', ...],
         'data': [[1.2302, 3.43, 2.34, 5.432, ...],
                    [3.23, 5.93, 1.03, 4.34, ...],
                    ...
                   ]
        }


    Parameters
    ----------
    expressionData : DataTypes.ExpressionData
        The object to convert

    Returns
    -------
    str
        The JSON representation of the input object

    """
    if(not isinstance(expressionData, ExpressionData)):
        raise ValueError("Input is not of type FastProject.DataTypes.ExpressionData");

    obj = {
        "row_labels": expressionData.row_labels,
        "col_labels": expressionData.col_labels,
        "data": expressionData.base.tolist()
    }

    return json.dumps(obj);


def signatures_to_JSON(signatures):
    """Returns the JSON representation of the
    input signatures

    Parameters
    ----------
    signatures : list of Signatures.Signature
        The object to convert

    Returns
    -------
    str
        The JSON representation of the input object

    """
    result = []
    for signature in signatures:
        if(not isinstance(signature, Signature)):
            raise ValueError("Input is not of type FastProject.Signatures.Signature");

        obj = {'sig_dict': signature.sig_dict,
               'signed': signature.signed,
               'source': signature.source,
               'name': signature.name};

        result.append(obj);

    return json.dumps(result);


def precomputed_signatures_to_JSON(precomputed_signatures):
    """Returns the JSON representation of the input precomputed
    signature dictionary

    Parameters
    ----------
    precomputed_signature : dict 
        Key: Name of Signature
        Value: SigScoreMethods.SignatureScores
        The object to convert

    Returns
    -------
    str
        The JSON representation of the input object

    """
    result = {};

    for key in precomputed_signatures:
        precomputed_signature = precomputed_signatures[key]

        if(not isinstance(precomputed_signature, SignatureScores)):
            raise ValueError("Input is not of type FastProject.SigScoreMethods.SignatureScores");

        obj = {
            "scores": precomputed_signature.scores,
            "name": precomputed_signature.name,
            "sample_labels": precomputed_signature.sample_labels,
            "isFactor": precomputed_signature.isFactor,
            "isPrecomputed": precomputed_signature.isPrecomputed,
            "numGenes": precomputed_signature.numGenes,
        }

        # Convert to a python list if its a numpy array
        if(isinstance(obj['scores'], np.ndarray)):
            obj['scores'] = obj['scores'].tolist();

        result[key] = obj;

    return json.dumps(result);


def JSON_to_ExpressionData(json_str):
    """Converts `json_str` into a FastProject ExpressionData object

    Parameters
    ----------
    json_str : str
        String representation of an ExpressionData Object to be converted

    Returns
    -------
    DataTypes.ExpressionData
        Resulting ExpressionData Object

    """
    json_obj = json.loads(json_str);

    data = json_obj['data'];
    data = np.array(data);
    row_labels = json_obj['row_labels'];
    col_labels = json_obj['col_labels'];

    return ExpressionData(data, row_labels, col_labels);


def JSON_to_Signatures(json_str):
    """Converts `json_str` into a list of FastProject Signature objects

    Parameters
    ----------
    json_str : str
        String representation of Signature objects to be converted

    Returns
    -------
    list of Signatures.Signature
        Resulting Signature Objects

    """
    json_obj = json.loads(json_str);

    result = [];

    for sig in json_obj:
        sig_dict = sig['sig_dict'];
        signed = sig['signed'];
        source = sig['source'];
        name = sig['name'];
        sig_obj = Signature(sig_dict, signed, source, name);

        result.append(sig_obj);

    return result;


def JSON_to_SignatureScores(json_str):
    """Converts `json_str` into a FastProject SignatureScore object

    Parameters
    ----------
    json_str : str
        String representation of a dictionary of
        SignatureScore objects to be converted

    Returns
    -------
    dict
        key: Name of Signature
        value: Signatures.SignatureScore

    """
    json_obj = json.loads(json_str);
    result = {};

    for key in json_obj:

        sig_score_dict = json_obj[key]

        scores = sig_score_dict['scores'];
        name = sig_score_dict['name'];
        sample_labels = sig_score_dict['sample_labels'];
        isFactor = sig_score_dict['isFactor'];
        isPrecomputed = sig_score_dict['isPrecomputed'];
        numGenes = sig_score_dict['numGenes'];

        sig_score_obj = SignatureScores(scores, name, sample_labels, isFactor, isPrecomputed, numGenes);

        result[key] = sig_score_obj;

    return result;
