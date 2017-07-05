# -*- coding: utf-8 -*-
""" Entry point for console fastproject script

- Parses arguments
- Loads data from files
- Runs analysis
- Writes outputs to file

Then launches the main pipeline
"""

from __future__ import absolute_import, print_function, division;
import argparse;
import os
import logging;
import numpy as np;
from . import FileIO;
from . import Signatures;
from . import Pipelines;
from . import HtmlViewer;
from .DataTypes import ExpressionData;
from .Global import FP_Output
import FastProject


def parseFPArgs():
    """Defines the command-line arguments and parses the FastProject call

    Returns
    -------
    argparse.Namespace

    """
    parser = argparse.ArgumentParser(prog="FastProject", description='Analyze a set of expression data with FastProject.');

    parser.add_argument("data_file", help="Input expression data matrix");

    parser.add_argument("-k", "--housekeeping", metavar="FILE", default="",
                      help="Read list of housekeeping genes from FILE.  Uses default list if not specified");

    parser.add_argument("-s", "--signatures", metavar="FILE", nargs='*',
                      help="Loads signatures from FILE.");

    parser.add_argument("-p", "--precomputed", metavar="FILE", nargs='*',
                      help="Loads precomputed signature scores from FILE.");

    parser.add_argument("-o", "--output", metavar="DIRECTORY",
            help="Name of output directory.  Otherwise, output directory is auto-generated");

    parser.add_argument("--nofilter",   action="store_true",
            help="Project using all genes.");

    parser.add_argument("--nomodel",   action="store_true",
            help="No estimation of expression probability or false negative probability");

    parser.add_argument("--pca_filter", action="store_true",
            help="Filters PC principal components that correlate with a calculated QC metric for each sample");

    parser.add_argument("--qc",         action="store_true",
            help="Performs a quality check on samples, filtering samples that do not pass");

    parser.add_argument("--all_sigs",         action="store_true",
            help="Do not remove insignificant signatures from output.");

    parser.add_argument("--debug",         action="store_true",
            help="Run FastProject in Debug mode");

    parser.add_argument("--lean",         action="store_true",
            help="Run FastProject in Lean mode - subset of analyses to reduce runtime on large data sets");

    parser.add_argument("--subsample_size", type=int, metavar="N", default=1000,
            help="Planned Feature: Number of samples to use when sub_sampling");

    parser.add_argument("--min_signature_genes", type=int, metavar="N", default=5,
            help="Signatures that match less than N genes in the data are discarded");

    parser.add_argument("--projections", metavar="FILE", nargs='*',
                        help="Loads projection coordinates from FILE");

    parser.add_argument("--weights", metavar="FILE",
                        help="Loads weights from FILE. Use these weights instead of FastProject's FNR calculation");

    parser.add_argument("--threshold", metavar="N", type=int,
                        help="Removes transcripts detected in less than N samples. " +
                        "Default is 20%% of total sample count.")

    parser.add_argument("--sig_norm_method",
        choices=["none", "znorm_columns", "znorm_rows",
            "znorm_rows_then_columns", "rank_norm_columns"],
        default="znorm_rows",
        help="Pick a normalization method to be applied to data before evaluating signature scores");

    parser.add_argument("--sig_score_method",
        choices=["naive", "weighted_avg", "imputed", "only_nonzero"],
        default="weighted_avg",
        help="Pick a method to evaluate signature scores");

    args = parser.parse_args();

    args = vars(args);  # Convert to a Dictionary

    return args;


def loadFilesFromDisk(args):
    """Loads files from disk into data structures

    Typically called before calling Pipelines.Analysis()

    Parameters
    ----------
    args : dict
        Contains file locations on disk

    Returns
    -------
    expressionMatrix : ExpressionData
    signatures : list of Signatures.Signature
    precomputed_signatures : dict
        Keys are precomputed signature names (str)
        Values are signature levels/scores (SigScoreMethods.SignatureScores)
    housekeeping_genes : list of str
    input_projections : dict
        Keys are of type str, representing projection names
        Values are of type 2xN pandas.DataFrame with column names matching
        sample names in `expressionMatrix`
    input_weights : pandas.DataFrame or None
        Same size as expressionMatrix
        Values are floats from 0.0 to 1.0

    """
    # Read expression data from file
    filename = args["data_file"];

    if(not os.path.isfile(filename)):
        raise ValueError("\n", filename, "not found.\nExiting...");

    (edata, genes, cells) = FileIO.read_matrix(filename);
    expressionMatrix = ExpressionData(edata, genes, cells);

    FP_Output("Imported ", edata.shape[0], " genes across ", edata.shape[1], " samples");

    # Load Signature files
    signatures = [];
    if(args["signatures"]):
        for sig_file in args["signatures"]:
            if(not os.path.isfile(sig_file)):
                raise ValueError("Option Error: signature file " + sig_file + " not found.\nExiting...");

            signatures += Signatures.read_signatures(sig_file);

    # Load Precomputed Sig file
    precomputed_signatures = {};
    if(args["precomputed"]):
        for precomputed_sig_file in args["precomputed"]:
            if(not os.path.isfile(precomputed_sig_file)):
                raise ValueError("Option Error: precomputed signature file " + precomputed_sig_file + " not found.\nExiting...");
            precomputed_signatures.update(Signatures.load_precomputed(precomputed_sig_file, cells));

    if(not args["signatures"] and not args["precomputed"]):  # Need one or the other here
        raise ValueError(
            "Option Error: Must specify either a signature file or a pre-computed signature file.\nExiting...");

    if(len(signatures) + len(precomputed_signatures) == 0):  # Need one or the other here
        raise ValueError(
            "Option Error: Must specify either a signature file or a pre-computed signature file.\nExiting...");

    # Load housekeeping genes
    housekeeping_genes = FileIO.load_housekeeping_genes(args["housekeeping"])

    # Load projection coordinates (if provided)
    input_projections = {};
    if(args["projections"]):
        input_projections = FileIO.load_input_projections(args["projections"], cells);

    # Load input weights (if provided)
    input_weights = None;
    if(args["weights"]):
        input_weights = FileIO.load_input_weights(args["weights"], genes, cells);

    return (expressionMatrix, signatures, precomputed_signatures,
            housekeeping_genes, input_projections, input_weights);


def createOutputDirectories(args):
    """
    Creates the output directory structure
    """

    # Create directory for all outputs
    if(args["output"]):
        dir_name = args["output"];
    else:
        default_dir_name = 'FastProject_Output';
        if(os.path.isdir(default_dir_name)):
            i = 1;
            while(True):
                dir_name = default_dir_name + str(i);
                if(not os.path.isdir(dir_name)):
                    break;
                else:
                    i = i + 1;
        else:
            dir_name = default_dir_name;

    FileIO.make_dirs(dir_name);

    logger = logging.getLogger("FastProject")
    logger.setLevel(logging.INFO);
    fh = logging.FileHandler(os.path.join(dir_name, 'fastproject.log'))
    fh.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
    logger.addHandler(fh)

    logger.info("Running FastProject version " + FastProject.__version__);
    logger.info("Using numpy version " + np.__version__);

    for key in args:
        logger.info(key + ": " + str(args[key]));

    return dir_name;


def entry():
    """Entry point for the fastproject command-line script
    """

    args = parseFPArgs();

    (expressionMatrix, signatures, precomputed_signatures,
     housekeeping_genes, input_projections,
     input_weights) = loadFilesFromDisk(args);

    try:
        dir_name = createOutputDirectories(args);  # Needs to be created first so logging can write here

        models, qc_info = Pipelines.Analysis(expressionMatrix, signatures, precomputed_signatures,
            housekeeping_genes, input_projections, input_weights, args);

        FileIO.saveResultstoDisk(models, signatures, qc_info, dir_name);
    except:
        import traceback;
        import sys;
        traceback.print_exc();

        tb_type, value, tb = sys.exc_info();

        if(args["debug"]):
            import pdb;
            pdb.post_mortem(tb);


