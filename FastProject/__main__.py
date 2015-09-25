# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 15:16:22 2015

@author: David
"""

from optparse import OptionParser;
import Global;

def parseFPArgs():
    parser = OptionParser(usage = 'usage: FastProject [options] data_file');

    parser.add_option("-k", "--housekeeping", metavar="FILE", default="",
                      help="Read list of housekeeping genes from FILE.  Uses default list if not specified");

    parser.add_option("-s", "--signatures", metavar="FILE",
                      help="Loads signatures from FILE.");

    parser.add_option("-p", "--precomputed", metavar="FILE",
                      help="Loads precomputed signature scores from FILE.");

    parser.add_option("-o", "--output", metavar="DIRECTORY", help="Name of output directory.  Otherwise, output directory is auto-generated");
    parser.add_option("--nofilter",   action="store_true", default=False, help="Project using all genes.");
    parser.add_option("--nomodel",   action="store_true", default=False, help="No estimation of expression probability or false negative probability");
    parser.add_option("--pca_filter", action="store_true", default=False, help="Filters PC principal components that correlate with a calculated QC metric for each sample");
    parser.add_option("--qc",         action="store_true", default=False, help="Performs a quality check on samples, filtering samples that do not pass");
    parser.add_option("--subsample_size", type="int", metavar="N", default=1000, help="Number of samples to use when sub_sampling. Default is 1000");

    (options, args) = parser.parse_args();

    return options, args;

def entry():
    options, args = parseFPArgs();
    Global.options = options;
    Global.args = args;

    from FastProject import Pipelines;
    Pipelines.FullOutput();

if(__name__ == "__main__"):
    entry();
