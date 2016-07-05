"""
A module for functions which convert
FastProject datatypes to/from JSON
"""

from .DataTypes import ExpressionData
from .Signatures import Signature, SignatureScores
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


def signature_to_JSON(signature):
    """Returns the JSON representation of the
    input signature

    Parameters
    ----------
    signature : Signatures.Signature
        The object to convert

    Returns
    -------
    str
        The JSON representation of the input object

    """
    if(not isinstance(signature, Signature)):
        raise ValueError("Input is not of type FastProject.Signatures.Signature");

    obj = {'sig_dict': signature.sig_dict,
           'signed': signature.signed,
           'source': signature.source,
           'name': signature.name};

    return json.dumps(obj);


def precomputed_signature_to_JSON(precomputed_signature):
    """Returns the JSON representation of the input precomputed
    signature

    Parameters
    ----------
    precomputed_signature : Signatures.SignatureScores
        The object to convert

    Returns
    -------
    str
        The JSON representation of the input object

    """
    if(not isinstance(precomputed_signature, SignatureScores)):
        raise ValueError("Input is not of type FastProject.Signatures.SignatureScores");

    obj = {
        "scores": precomputed_signature.scores,
        "name": precomputed_signature.name,
        "sample_labels": precomputed_signature.sample_labels,
        "isFactor": precomputed_signature.isFactor,
        "isPrecomputed": precomputed_signature.isPrecomputed,
        "numGenes": precomputed_signature.numGenes,
    }

    return json.dumps(obj);


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


def JSON_to_Signature(json_str):
    """Converts `json_str` into a FastProject Signature object

    Parameters
    ----------
    json_str : str
        String representation of an Signature Object to be converted

    Returns
    -------
    Signatures.Signature
        Resulting Signature Object

    """
    json_obj = json.loads(json_str);

    sig_dict = json_obj['sig_dict'];
    signed = json_obj['signed'];
    source = json_obj['source'];
    name = json_obj['name'];

    return Signature(sig_dict, signed, source, name);


def JSON_to_SignatureScore(json_str):
    """Converts `json_str` into a FastProject SignatureScore object

    Parameters
    ----------
    json_str : str
        String representation of an SignatureScore Object to be converted

    Returns
    -------
    Signatures.SignatureScore
        Resulting SignatureScore Object

    """
    json_obj = json.loads(json_str);

    scores = json_obj['scores'];
    name = json_obj['name'];
    sample_labels = json_obj['sample_labels'];
    isFactor = json_obj['isFactor'];
    isPrecomputed = json_obj['isPrecomputed'];
    numGenes = json_obj['numGenes'];

    return SignatureScores(scores, name, sample_labels, isFactor, isPrecomputed, numGenes);
