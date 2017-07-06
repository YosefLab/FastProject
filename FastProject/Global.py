# -*- coding: utf-8 -*-
"""Contains information that is useful to share
across all modules.

Right now, mainly command-line arguments, logging
and resource locations.
"""
from __future__ import absolute_import, print_function, division;
import logging;
import sys
import os
import random;

# Defaults for arguments are defined here
# These are all overwritten if called from the command line
# However, by putting defaults here, it makes it easier to
#    call FastProject from within other scripts
# These should be kept in sync with arguments in __main__.py
# args = namedtuple("defalt_args",  ["data_file",  "housekeeping",
#                                    "signatures",  "precomputed",  "output",
#                                    "nofilter",  "nomodel",  "pca_filter",
#                                    "qc",  "subsample_size",
#                                    "min_signature_genes",  "projections",
#                                    "weights",  "threshold"]);
# 
# args.data_file = "";
# args.housekeeping = "";
# args.signatures = [];
# args.precomputed = [];
# args.output = "";
# args.nofilter = False;
# args.nomodel = False;
# args.pca_filter = False;
# args.qc = False;
# args.subsample_size = None;
# args.min_signature_genes = 5;
# args.projections = [];
# args.weights = "";
# args.threshold = None;


logger = logging.getLogger("FastProject");
def FP_Output(*args):
    """
    Used to have finer control over outputs.
    """
    print(*args);
    logmessage = ' '.join([str(a) for a in args]);
    if(logmessage.startswith("\n")):
        logmessage = logmessage[1:];
    logger.info(logmessage);

# This section for finding resource files
# Behavior is different depending on whether or not we are running frozen

if getattr(sys, 'frozen', False):
    this_directory = sys._MEIPASS;
else:
    this_directory = os.path.dirname(os.path.abspath(__file__));


def get_viewer_resource_dir():
    return os.path.join(this_directory, "Viewer Resources");


def get_housekeeping_dir():
    return os.path.join(this_directory, "Housekeeping Genes");

# Chosen by roll of a 2147483648-sided die
# Guaranteed to be random
RANDOM_SEED = 1335607828;
random.seed(RANDOM_SEED);
