# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 15:16:22 2015

@author: David
"""

from optparse import OptionParser;
from FastProject import Pipelines;

def parseFPArgs():
    parser = OptionParser(usage = 'usage: FastProject [options] data_file');
    parser.add_option("-k", "--housekeeping", metavar="FILE",
                      help="Read list of housekeeping genes from FILE.  Uses default list if not specified");
    parser.add_option("-s", "--signatures", metavar="FILE",
                      help="Loads signatures from FILE.  Otherwise, signature analysis is skipped in unless in interactive mode.");
    parser.add_option("--precomputed", metavar="FILE",
                      help="Loads precomputed signature scores from FILE.");
    parser.add_option("-f","--filters",default="0",help="""Specifies filters to be used on genes\n\n1. Remove Housekeeping
    2. Threshold (Gene expressed in at least 20% of samples)
    3. Bimodal (Using Hartigans Dip Test p<0.05

    e.g. -f 1, -f 3, -f 13, -f 123""");
    parser.add_option("-p","--probability", action="store_true", default=False, help="Projects using probability of expression rather than log expression level");
    parser.add_option("-c","--pca", action="store_true", default=False, help="Reduces to a smaller number of principal components before projection");
    parser.add_option("--pca_filter", action="store_true", default=False, help="Filters PC principal components that correlate with a calculated QC metric for each sample");
    parser.add_option("-o", "--output", action="store_true", default=False, help="Outputs data after filtering and any transforms");
    parser.add_option("-q", "--qc", action="store_true", default=False, help="Performs a quality check on samples, filtering samples that do not pass");
    parser.add_option("-i", "--interactive", action="store_true", default=False, help="Prompts options via command line instead");
    parser.add_option("-a", "--all", action="store_true", default=False, help="Generate viewer with all projection combinations");

    (options, args) = parser.parse_args();

    return options, args;

def entry():
    options, args = parseFPArgs();

    if(options.all):
        Pipelines.FullOutput(options, args);
    else:
        Pipelines.SingleOutput(options, args);

if(__name__ == "__main__"):
    entry();
