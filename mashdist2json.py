#!/usr/bin/env python3


"""
Purpose
-------

This module is intended to generate a json output for mash dist results that
can be imported in pATLAS.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``mash_output`` : String with the name of the mash screen output file.
    - e.g.: ``'fastaFileA_mashdist.txt'``


Code documentation
------------------

"""

__version__ = "1.0.0"
__build__ = "08022018"
__template__ = "mashsdist2json-nf"

import sys
import os
import json
import traceback

try:
    sys.path.append(os.environ["ASSEMBLERFLOW_UTILS"])
except KeyError:
    pass

from utils.assemblerflow_base import get_logger, _log_error

logger = get_logger(__file__)

def build_versions():
    logger.debug("Checking module versions")

    ver = [{
        "program": __template__,
        "version": __version__,
        "build": __build__
    }]
    logger.debug("Versions list set to: {}".format(ver))

    with open(".versions", "w") as fh:
        fh.write(json.dumps(ver, separators=(",", ":")))

if __file__.endswith(".command.sh"):
    MASH_TXT = '$depthFile'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("MASH_TXT: {}".format(MASH_TXT))

if __file__.endswith(".command.sh"):
    MASH_TXT = '$mashtxt'

def main(mash_output):
    '''
    Main function that allows to dump a mash dist txt file to a json file

    Parameters
    ----------
    mash_output: str
        A string with the input file.

    '''
    out_file = open(" ".join(mash_output.split(".")[:-1]) + ".json", "w")
    input_f = open(mash_output, 'r')

    master_dict = {}
    for line in input_f:
        tab_split = line.split("\t")
        ref_accession = "_".join(tab_split[0].strip().split("_")[0:3])
        mash_dist = tab_split[2].strip()
        ## there is no need to store all values since we are only interested in
        # representing the significant ones
        ## and those that correlate well with ANI (mashdist<=0.1)
        master_dict[ref_accession] = 1 - float(mash_dist)
    ## writes output json
    out_file.write(json.dumps(master_dict))
    out_file.close()

if __name__ == "__main__":
    try:
        build_versions()
        # a variable from nextflow process
        main(MASH_TXT)
    except Exception:
        logger.error("Module exited unexpectedly with error:\\n{}".format(
            traceback.format_exc()))
        _log_error()