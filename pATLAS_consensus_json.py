#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to generate a json output from the consensus results from
all the approaches available through options (mapping, assembly, mash screen)

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``mapping_json`` : String with the name of the json file with mapping results.
    - e.g.: ``'mapping_SampleA.json'``
- ``dist_json`` : String with the name of the json file with mash dist results.
    - e.g.: ``'mash_dist_SampleA.json'``
- ``screen_json`` : String with the name of the json file with mash screen results.
    - e.g.: ``'mash_screen_sampleA.json'``


Code documentation
------------------

"""

__version__ = "0.1.0"
__build__ = "24022018"
__template__ = "pATLAS_consensus_json-nf"

import os
import json

from utils.assemblerflow_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    MAPPING_JSON = '$mappingOutputFile'
    DIST_JSON = '$mashDistOutputFile'
    SCREEN_JSON = '$mashScreenOutputFile'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("MAPPING_JSON: {}".format(MAPPING_JSON))
    logger.debug("DIST_JSON: {}".format(DIST_JSON))
    logger.debug("SCREEN_JSON: {}".format(SCREEN_JSON))

@MainWrapper
def main(mapping_json, dist_json, screen_json):
    """

    Parameters
    ----------
    mapping_json: str
        The path to the json file with mapping results
    dist_json: str
        The path to the json file with the mash distance results
    screen_json: str
        The path to the json file with the mash screen results

    """

    print(mapping_json, dist_json, screen_json)

    # TODO continue implementation of consensus json


if __name__ == "__main__":
    main(MAPPING_JSON, DIST_JSON, SCREEN_JSON)