from __future__ import division, print_function;
import logging;

args = dict();


def FP_Output(*args):
    """
    Used to have finer control over outputs.
    """
    print(*args);
    logmessage = ' '.join([str(a) for a in args]);
    logging.info(logmessage);
