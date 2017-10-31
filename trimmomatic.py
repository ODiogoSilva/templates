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

"""

# TODO: More control over read trimming
# TODO: Add option to remove adapters
# TODO: What to do when there is encoding failure

import subprocess

from subprocess import PIPE


if __file__.endswith(".command.sh"):
    FASTQ_ID = '$fastq_id'
    FASTQ_PAIR = '$fastq_pair'.split()
    TRIM_RANGE = '$trim_range'.split()
    TRIM_OPTS = [x.strip() for x in '$opts'.strip("[]").split(",")]
    PHRED = '$phred'


TRIM_PATH = "/NGStools/Trimmomatic-0.36/trimmomatic.jar"


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

    # Add trimmomatic options
    cli += [
        "CROP:{}".format(trim_range[1]),
        "HEADCROP:{}".format(trim_range[0]),
        "SLIDINGWINDOW:{}".format(trim_opts[0]),
        "LEADING:{}".format(trim_opts[1]),
        "TRAILING:{}".format(trim_opts[2]),
        "MINLEN:{}".format(trim_opts[3]),
        "TOPHRED33",
        "-trimlog",
        "{}_trimlog.txt".format(fastq_id)
    ]

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # Check if trimmomatic ran successfully. If not, write the error message
    # to the status channel and exit.
    with open("trimmomatic_status", "w") as fh:
        if p.returncode != 0:
            fh.write(str(stderr))
            return
        else:
            fh.write("pass")


if __name__ == '__main__':
    main(FASTQ_ID, FASTQ_PAIR, TRIM_RANGE, TRIM_OPTS, PHRED)
