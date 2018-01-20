#!/usr/bin/env python3

"""
Purpose
-------

This module is intended execute trimmomatic on paired-end FastQ files.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``fastq_id`` : Pair of FastQ file paths.
    - e.g.: ``'SampleA'``
- ``fastq_pair`` : Pair of FastQ file paths.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'``
- ``trim_range`` : Crop range detected using FastQC.
    - e.g.: ``'15 151'``
- ``opts`` : List of options for trimmomatic
    - e.g.: ``'["5:20", "3", "3", "55"]'``
    - e.g.: ``'[trim_sliding_window, trim_leading, trim_trailing, trim_min_length]'``
- ``phred`` : List of guessed phred values for each sample
    - e.g.: ``'[SampleA: 33, SampleB: 33]'``

Generated output
----------------

The generated output are output files that contain an object, usually a string.
(Values within ``${}`` are substituted by the corresponding variable.)

- ``${fastq_id}_*P*``: Pair of paired FastQ files generated by Trimmomatic
    - e.g.: ``'SampleA_1_P.fastq.gz SampleA_2_P.fastq.gz'``
- ``trimmomatic_status``: Stores the status of the trimmomatic run. If it was\
    successfully executed, it stores 'pass'. Otherwise, it stores the \
    ``STDERR`` message.
    - e.g.: ``'pass'``

Code documentation
------------------

"""

# TODO: More control over read trimming
# TODO: Add option to remove adapters
# TODO: What to do when there is encoding failure

__version__ = "1.0.0"
__build__ = "16012018"
__template__ = "trimmomatic-nf"

import os
import json
import logging
import subprocess

from subprocess import PIPE
from collections import OrderedDict


# create logger
logger = logging.getLogger('simple_example')
logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)


def build_versions():

    def get_trimmomatic_version():

        try:

            cli = ["java", "-jar", TRIM_PATH, "-version"]
            p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
            stdout, _ = p.communicate()

            version = stdout.strip().decode("utf8")

        except Exception as e:
            logger.debug(e)
            version = "undefined"

        return {
            "program": "Trimmomatic",
            "version": version,
        }

    logger.debug("Checking module versions")

    ver = [{
        "program": __template__,
        "version": __version__,
        "build": __build__
    }, get_trimmomatic_version()]
    logger.debug("Versions list set to: {}".format(ver))

    with open(".versions", "w") as fh:
        fh.write(json.dumps(ver, separators=(",", ":")))


if __file__.endswith(".command.sh"):
    FASTQ_ID = '$fastq_id'
    FASTQ_PAIR = '$fastq_pair'.split()
    TRIM_RANGE = '$trim_range'.split()
    TRIM_OPTS = [x.strip() for x in '$opts'.strip("[]").split(",")]
    PHRED = '$phred'

    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("FASTQ_ID: {}".format(FASTQ_ID))
    logger.debug("FASTQ_PAIR: {}".format(FASTQ_PAIR))
    logger.debug("TRIM_RANGE: {}".format(TRIM_RANGE))
    logger.debug("TRIM_OPTS: {}".format(TRIM_OPTS))
    logger.debug("PHRED: {}".format(PHRED))

TRIM_PATH = "/NGStools/Trimmomatic-0.36/trimmomatic.jar"


def _log_error():
    """Nextflow specific function that logs an error upon unexpected failing
    """

    import traceback

    with open(".status", "w") as status_fh:
        logger.error("Module exited unexpectedly with error:\\n{}".format(
            traceback.format_exc()))
        status_fh.write("error")


def parse_log(log_file):
    """Retrieves some statistics from a single Trimmomatic log file.

    This function parses Trimmomatic's log file and stores some trimming
    statistics in an :py:class:`OrderedDict` object. This object contains
    the following keys:

        - ``clean_len``: Total length after trimming.
        - ``total_trim``: Total trimmed base pairs.
        - ``total_trim_perc``: Total trimmed base pairs in percentage.
        - ``5trim``: Total base pairs trimmed at 5' end.
        - ``3trim``: Total base pairs trimmed at 3' end.

    Parameters
    ----------
    log_file : str
        Path to trimmomatic log file.

    Returns
    -------
    x : :py:class:`OrderedDict`
        Object storing the trimming statistics.

    """

    template = OrderedDict([
        # Total length after trimming
        ("clean_len", 0),
        # Total trimmed base pairs
        ("total_trim", 0),
        # Total trimmed base pairs in percentage
        ("total_trim_perc", 0),
        # Total trimmed at 5' end
        ("5trim", 0),
        # Total trimmed at 3' end
        ("3trim", 0),
        # Bad reads (completely trimmed)
        ("bad_reads", 0)
    ])

    with open(log_file) as fh:

        for line in fh:
            # This will split the log fields into:
            # 0. read length after trimming
            # 1. amount trimmed from the start
            # 2. last surviving base
            # 3. amount trimmed from the end
            fields = [int(x) for x in line.strip().split()[-4:]]

            if not fields[0]:
                template["bad_reads"] += 1

            template["5trim"] += fields[1]
            template["3trim"] += fields[3]
            template["total_trim"] += fields[1] + fields[3]
            template["clean_len"] += fields[0]

        total_len = template["clean_len"] + template["total_trim"]

        if total_len:
            template["total_trim_perc"] = round(
                (template["total_trim"] / total_len) * 100, 2)
        else:
            template["total_trim_perc"] = 0

    return template


def write_report(storage_dic, output_file):
    """ Writes a report from multiple samples.

    Parameters
    ----------
    storage_dic : dict or :py:class:`OrderedDict`
        Storage containing the trimming statistics. See :py:func:`parse_log`
        for its generation.
    output_file : str
        Path where the output file will be generated.
    """

    with open(output_file, "w") as fh, open(".report.json", "w") as json_rep:

        # Write header
        fh.write("Sample,Total length,Total trimmed,%,5end Trim,3end Trim,"
                 "bad_reads\\n")

        # Write contents
        for sample, vals in storage_dic.items():
            fh.write("{},{}\\n".format(
                sample, ",".join([str(x) for x in vals.values()])))

            json_dic = {
                "tableRow": [
                    {"header": "trimmed",
                     "value": vals["total_trim_perc"],
                     "table": "assembly",
                     "columnBar": True},
                    ],
                "plotData": {
                    "sparkline": vals["clean_len"]
                },
                "badReads": vals["bad_reads"]
            }
            json_rep.write(json.dumps(json_dic, separators=(",", ":")))


def trimmomatic_log(log_file):

    log_storage = OrderedDict()

    log_id = log_file.rstrip("_trimlog.txt")

    log_storage[log_id] = parse_log(log_file)

    os.remove(log_file)

    write_report(log_storage, "trimmomatic_report.csv")


def clean_up():
    """Cleans the working directory of unwanted temporary files"""

    # Find unpaired fastq files
    unpaired_fastq = [f for f in os.listdir(".")
                      if f.endswith("_U.fastq.gz")]

    # Remove unpaired fastq files, if any
    for fpath in unpaired_fastq:
        os.remove(fpath)


def main(fastq_id, fastq_pair, trim_range, trim_opts, phred):
    """ Main executor of the trimmomatic template.

    Parameters
    ----------
    fastq_id : str
        Sample Identification string.
    fastq_pair : list
        Two element list containing the paired FastQ files.
    trim_range : list
        Two element list containing the trimming range.
    trim_opts : list
        Four element list containing several trimmomatic options:
        [*SLIDINGWINDOW*; *LEADING*; *TRAILING*; *MINLEN*]
    phred : int
        Guessed phred score for the sample. The phred score is a generated
        output from :py:class:`templates.integrity_coverage`.

    """

    logger.info("Starting trimmomatic")

    # Create base CLI
    cli = [
        "java",
        "-Xmx{}".format("$task.memory"[:-1].lower().replace(" ", "")),
        "-jar",
        TRIM_PATH.strip(),
        "PE",
        "-threads",
        "$task.cpus"
    ]

    # If the phred encoding was detected, provide it
    try:
        # Check if the provided PHRED can be converted to int
        phred = int(phred)
        phred_flag = "-phred{}".format(str(phred))
        cli += [phred_flag]
    # Could not detect phred encoding. Do not add explicit encoding to
    # trimmomatic and let it guess
    except ValueError:
        pass

    # Add input samples to CLI
    cli += fastq_pair

    # Add output file names
    output_names = []
    for i in range(len(fastq_pair)):
        output_names.append("{}_{}_P.fastq.gz".format(FASTQ_ID, str(i + 1)))
        output_names.append("{}_{}_U.fastq.gz".format(FASTQ_ID, str(i + 1)))
    cli += output_names

    if trim_range != ["None"]:
        cli += [
            "CROP:{}".format(trim_range[1]),
            "HEADCROP:{}".format(trim_range[0]),
        ]

    # Add trimmomatic options
    cli += [
        "SLIDINGWINDOW:{}".format(trim_opts[0]),
        "LEADING:{}".format(trim_opts[1]),
        "TRAILING:{}".format(trim_opts[2]),
        "MINLEN:{}".format(trim_opts[3]),
        "TOPHRED33",
        "-trimlog",
        "{}_trimlog.txt".format(fastq_id)
    ]

    logger.debug("Running trimmomatic subprocess with command: {}".format(cli))

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # Attempt to decode STDERR output from bytes. If unsuccessful, coerce to
    # string
    try:
        stderr = stderr.decode("utf8")
    except (UnicodeDecodeError, AttributeError):
        stderr = str(stderr)

    logger.info("Finished trimmomatic subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished trimmomatic subprocesswith STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished trimmomatic with return code: {}".format(
        p.returncode))

    trimmomatic_log("{}_trimlog.txt".format(fastq_id))

    clean_up()

    # Check if trimmomatic ran successfully. If not, write the error message
    # to the status channel and exit.
    with open(".status", "w") as status_fh:
        if p.returncode != 0:
            status_fh.write("fail")
            return
        else:
            status_fh.write("pass")


if __name__ == '__main__':

    try:
        build_versions()
        main(FASTQ_ID, FASTQ_PAIR, TRIM_RANGE, TRIM_OPTS, PHRED)
    except Exception:
        _log_error()
