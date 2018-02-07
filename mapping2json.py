#!/usr/bin/env python3

import json

if __file__.endswith(".command.sh"):
    DEPTH_TXT = '$depthFile'
    JSON_LENGTH = '$lengthJson'
    CUTOFF = '$params.cov_cutoff'

def depthfilereader(depth_file, plasmid_length, cutoff):
    '''
    Function that parse samtools depth file and creates 3 dictionaries that
    will be useful to make the outputs of this script, both the tabular file
    and the json file that may be imported by pATLAS

    Parameters
    ----------
    depth_file: str
        the path to depth file for each sample
    plasmid_length: dict
        a dictionary that stores length of all plasmids in fasta given as input
    cutoff: str
        the cutoff used to trim the unwanted matches for the minimum coverage
        results from mapping. This is then converted into a float within this
        function in order to compare with the value returned from the perc_value_per_ref.

    Returns
    -------
    percentage_basescovered: dict
            stores the percentage of the total sequence of a
            reference/accession (plasmid) in a dictionary
    '''
    depth_dic_coverage = {}
    for line in depth_file:
        tab_split = line.split()    # split by any white space
        reference = "_".join(tab_split[0].strip().split("_")[0:3])  # store
        # only the gi for the reference
        position = tab_split[1]
        numreadsalign = float(tab_split[2].rstrip("\n"))
        if reference not in depth_dic_coverage:
            depth_dic_coverage[reference] = {}
        depth_dic_coverage[reference][position] = numreadsalign

    percentage_basescovered = {}
    for ref in depth_dic_coverage:
        # calculates the percentage value per each reference
        perc_value_per_ref = float(len(depth_dic_coverage[ref])) / \
                                       float(plasmid_length[ref])
        # checks if percentage value is higher or equal to the cutoff defined
        if perc_value_per_ref >= float(cutoff):
            percentage_basescovered[ref] = perc_value_per_ref

    return percentage_basescovered

def main(depth_file, json_dict, cutoff):
    '''
    Function that handles the inputs required to parse depth files from bowtie
    and dumps a dict to a json file that can be imported into pATLAS.

    Parameters
    ----------
    depth_file: str
         the path to depth file for each sample
    json_dict: str
        the file that contains the dictionary with keys and values for accessions
        and their respective lengths
    cutoff: str
        the cutoff used to trim the unwanted matches for the minimum coverage
        results from mapping


    '''
    # loads dict from file, this file is provided in docker image

    plasmid_length = json.load(open(json_dict))

    # read depth file
    depth_file_reader = open(depth_file)

    # first reads the depth file and generates dictionaries to handle the input
    # to a simpler format
    percentage_basescovered = depthfilereader(depth_file_reader, plasmid_length, cutoff)

    # then dump do file
    output_json = open(depth_file + ".json", "w")
    output_json.write(json.dumps(percentage_basescovered))
    output_json.close()


if __name__ == "__main__":
    # a variable from nextflow process
    main(DEPTH_TXT, JSON_LENGTH, CUTOFF)
