from __future__ import division, print_function;
import logging;

options = dict();
args = [];


def FP_Output(*args):
    """
    Used to have finer control over outputs.
    """
    print(*args);
    logmessage = ' '.join([str(a) for a in args]);
    logging.info(logmessage);
