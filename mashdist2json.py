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

__version__ = "1.0.1"
__build__ = "20022018"
__template__ = "mashsdist2json-nf"

import os
import json

from utils.assemblerflow_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    MASH_TXT = '$mashtxt'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("MASH_TXT: {}".format(MASH_TXT))

@MainWrapper
def main(mash_output):
    '''
    Main function that allows to dump a mash dist txt file to a json file

    Parameters
    ----------
    mash_output: str
        A string with the input file.

    '''
    out_file = open(".".join(mash_output.split(".")[:-1]) + ".json", "w")
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

    main(MASH_TXT)
