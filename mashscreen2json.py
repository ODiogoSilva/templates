#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to generate a json output for mash screen results that
can be imported in pATLAS.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``mash_output`` : String with the name of the mash screen output file.
    - e.g.: ``'sortedMashScreenResults_SampleA.txt'``


Code documentation
------------------

"""

__version__ = "1.0.0"
__build__ = "08022018"
__template__ = "mashscreen2json-nf"

from statistics import median
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
    MASH_TXT = '$mashtxt'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("MASH_TXT: {}".format(MASH_TXT))

def main(mash_output):
    '''
    converts top results from mash screen txt output to json format

    Parameters
    ----------
    mash_output: str
        this is a string that stores the path to this file, i.e, the name of
        the file

    '''
    logger.info("Reading file : {}".format(mash_output))
    read_mash_output = open(mash_output)

    dic = {}
    median_list = []

    logger.info("Generating dictionary and list to pre-process the final json")
    for line in read_mash_output:
        tab_split = line.split("\t")
        identity = tab_split[0]
        #shared_hashes = tab_split[1]
        median_multiplicity = tab_split[2]
        #p_value = tab_split[3]
        query_id = tab_split[4]
        #query-comment should not exist here and it is irrelevant

        # here identity is what in fact interests to report to json but
        # median_multiplicity also is important since it gives an rough
        # estimation of the coverage depth for each plasmid.
        # Plasmids should have higher coverage depth due to their increased
        # copy number in relation to the chromosome.
        dic[query_id] = [identity, median_multiplicity]
        median_list.append(float(median_multiplicity))

    output_json = open(" ".join(mash_output.split(".")[:-1]) + ".json", "w")

    # median cutoff is twice the median of all median_multiplicity values
    # reported by mash screen. In the case of plasmids, since the database
    # has 9k entries and reads shouldn't have that many sequences it seems ok...
    if len(median_list) > 0:
        # this statement assures that median_list has indeed any entries
        median_cutoff = median(median_list)
        logger.info("Generating final json to dump to a file")
        filtered_dic = {}
        for k, v in dic.items():
            # estimated copy number
            copy_number = int(float(v[1]) / median_cutoff)
            # assure that plasmid as at least twice the median coverage depth
            if float(v[1]) > median_cutoff:
                filtered_dic["_".join(k.split("_")[0:3])] = [v[0],
                                                             str(copy_number)]
        logger.info(
            "Exported dictionary has {} entries".format(len(filtered_dic)))
        output_json.write(json.dumps(filtered_dic))
    else:
        # if no entries were found raise an error
        logger.error("No matches were found using mash screen for the queried reads")

    output_json.close()

if __name__ == "__main__":
    try:
        build_versions()
        # a variable from nextflow process
        main(MASH_TXT)
    except Exception:
        logger.error("Module exited unexpectedly with error:\\n{}".format(
            traceback.format_exc()))
        _log_error()