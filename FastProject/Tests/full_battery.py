import os;
from collections import namedtuple;
from ..Pipelines import Analysis;
from ..CLI import loadFilesFromDisk, createOutputDirectories, saveResultstoDisk;
this_directory = os.path.dirname(os.path.abspath(__file__));


def get_default_args():
    args = namedtuple("defalt_args",  ["data_file",  "housekeeping",
                                       "signatures",  "precomputed",  "output",
                                       "nofilter",  "nomodel",  "pca_filter",
                                       "qc",  "subsample_size",
                                       "min_signature_genes",  "projections",
                                       "weights",  "threshold"
                                       "all_sigs", "debug", "lean",
                                       "sig_norm_method", "sig_score_method"]);

    args.data_file = "";
    args.housekeeping = "";
    args.signatures = [];
    args.precomputed = [];
    args.output = "";
    args.nofilter = False;
    args.nomodel = False;
    args.pca_filter = False;
    args.qc = False;
    args.subsample_size = None;
    args.min_signature_genes = 5;
    args.projections = [];
    args.weights = None;
    args.threshold = None;
    args.all_sigs = False;
    args.debug = False;
    args.lean = False;
    args.sig_norm_method = "znorm_rows";
    args.sig_score_method = "weighted_avg";

    return args;


def run_simple():

    args = get_default_args();

    # Input Data Files
    args.data_file = os.path.join(this_directory, "TestFiles", "smallData.txt");
    args.signatures = [os.path.join(this_directory, "TestFiles", "sigsSmall.txt")];
    args.precomputed = [os.path.join(this_directory, "TestFiles", "precomputed_sigs.txt")];

    (expressionMatrix, signatures, precomputed_signatures,
     housekeeping_genes, input_projections,
     input_weights) = loadFilesFromDisk(args);

    dir_name = createOutputDirectories(args);  # Needs to be created first so logging can write here

    models, qc_info = Analysis(expressionMatrix, signatures, precomputed_signatures,
        housekeeping_genes, input_projections, input_weights, args);

    saveResultstoDisk(models, signatures, qc_info, dir_name);

    # Cleanup
    import shutil
    import time
    for x in range(10):  # Solution to Dropbox locking files.
        try:
            shutil.rmtree(dir_name);
        except:
            pass;
        else:
            break;
        time.sleep(0.1);
