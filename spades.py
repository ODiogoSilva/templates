#!/usr/bin/env python3

"""
Purpose
-------

This module is intended execute Spades on paired-end FastQ files.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``fastq_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fastq_pair`` : Pair of FastQ file paths.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'``
- ``kmers`` : Setting for Spades kmers. Can be either ``'auto'``, \
    ``'default'`` or a user provided list.
    - e.g.: ``'auto'`` or ``'default'`` or ``'55 77 99 113 127'``
- ``opts`` : List of options for spades execution.
    1. The minimum number of reads to consider an edge in the de Bruijn \
    graph during the assembly.
        - e.g.: ``'5'``
    2. Minimum contigs k-mer coverage.
        - e.g.: ``['2' '2']``

Generated output
----------------

- ``contigs.fasta`` : Main output of spades with the assembly
    - e.g.: ``contigs.fasta``
- ``spades_status`` :  Stores the status of the spades run. If it was \
    successfully executed, it stores ``'pass'``. Otherwise, it stores the\
    ``STDERR`` message.
    - e.g.: ``'pass'``

Code documentation
------------------

"""

__version__ = "1.0.0"
__build__ = "16012018"
__template__ = "spades-nf"

import os
import json
import logging
import subprocess

from subprocess import PIPE


# create logger
logger = logging.getLogger(os.path.basename(__file__))
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

    def get_spades_version():

        try:

            cli = ["spades.py", "--version"]
            p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
            stdout, _ = p.communicate()

            version = stdout.strip().split()[-1][1:].decode("utf8")

        except Exception as e:
            logger.debug(e)
            version = "undefined"

        return {
            "program": "SPAdes",
            "version": version,
        }

    logger.debug("Checking module versions")

    ver = [{
        "program": __template__,
        "version": __version__,
        "build": __build__
    }, get_spades_version()]
    logger.debug("Versions list set to: {}".format(ver))

    with open(".versions", "w") as fh:
        fh.write(json.dumps(ver, separators=(",", ":")))


if __file__.endswith(".command.sh"):
    FASTQ_ID = '$fastq_id'
    FASTQ_PAIR = '$fastq_pair'.split()
    MAX_LEN = int('$max_len'.strip())
    KMERS = '$kmers'.strip()
    OPTS = [x.strip() for x in '$opts'.strip("[]").split(",")]
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("FASTQ_ID: {}".format(FASTQ_ID))
    logger.debug("FASTQ_PAIR: {}".format(FASTQ_PAIR))
    logger.debug("MAX_LEN: {}".format(MAX_LEN))
    logger.debug("KMERS: {}".format(KMERS))
    logger.debug("OPTS: {}".format(OPTS))


def _log_error():
    """Nextflow specific function that logs an error upon unexpected failing
    """

    import traceback

    with open(".status", "w") as status_fh:
        logger.error("Module exited unexpectedly with error:\\n{}".format(
            traceback.format_exc()))
        status_fh.write("error")


def set_kmers(kmer_opt, max_read_len):
    """Returns a kmer list based on the provided kmer option and max read len.

    Parameters
    ----------
    kmer_opt : str
        The k-mer option. Can be either ``'auto'``, ``'default'`` or a
        sequence of space separated integers, ``'23, 45, 67'``.
    max_read_len : int
        The maximum read length of the current sample.

    Returns
    -------
    kmers : list
        List of k-mer values that will be provided to Spades.

    """

    logger.debug("Kmer option set to: {}".format(kmer_opt))

    # Check if kmer option is set to auto
    if kmer_opt == "auto":

        if max_read_len >= 175:
            kmers = [55, 77, 99, 113, 127]
        else:
            kmers = [21, 33, 55, 67, 77]

        logger.debug("Kmer range automatically selected based on max read"
                     "length of {}: {}".format(max_read_len, kmers))

    # Check if manual kmers were specified
    elif len(kmer_opt.split()) > 1:

        kmers = kmer_opt.split()
        logger.debug("Kmer range manually set to: {}".format(kmers))

    else:

        kmers = []
        logger.debug("Kmer range set to empty (will be automatically "
                     "determined by SPAdes")

    return kmers


def main(fastq_id, fastq_pair, max_len, kmer, opts):
    """Main executor of the spades template.

    Parameters
    ----------
    fastq_id : str
        Sample Identification string.
    fastq_pair : list
        Two element list containing the paired FastQ files.
    max_len : int
        Maximum read length. This value is determined in
        :py:class:`templates.integrity_coverage`
    kmer : str
        Can be either ``'auto'``, ``'default'`` or a
        sequence of space separated integers, ``'23, 45, 67'``.
    opts : List of options for spades execution. See above.

    """

    logging.info("Starting spades")

    min_coverage, min_kmer_coverage = opts

    logging.info("Setting SPAdes kmers")
    kmers = set_kmers(kmer, max_len)
    logging.info("SPAdes kmers set to: {}".format(kmers))

    cli = [
        "spades.py",
        "--careful",
        "--only-assembler",
        "--threads",
        "$task.cpus",
        "--cov-cutoff",
        min_coverage,
        "-o",
        "."
    ]

    # Add kmers, if any were specified
    if kmers:
        cli += ["-k {}".format(",".join([str(x) for x in kmers]))]

    # Add FastQ files
    cli += [
        "-1",
        fastq_pair[0],
        "-2",
        fastq_pair[1]
    ]

    logger.debug("Running SPAdes subprocess with command: {}".format(cli))

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # Attempt to decode STDERR output from bytes. If unsuccessful, coerce to
    # string
    try:
        stderr = stderr.decode("utf8")
        stdout = stdout.decode("utf8")
    except (UnicodeDecodeError, AttributeError):
        stderr = str(stderr)
        stdout = str(stdout)

    logger.info("Finished SPAdes subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished SPAdes subprocesswith STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished SPAdes with return code: {}".format(
        p.returncode))

    with open(".status", "w") as fh:
        if p.returncode != 0:
            fh.write("error")
            return
        else:
            fh.write("pass")

    # Change the default contigs.fasta assembly name to a more informative one
    assembly_file = "{}_spades.assembly.fasta".format(fastq_id)
    os.rename("contigs.fasta", assembly_file)
    logger.info("Setting main assembly file to: {}".format(assembly_file))


if __name__ == '__main__':

    try:
        build_versions()
        main(FASTQ_ID, FASTQ_PAIR, MAX_LEN, KMERS, OPTS)
    except:
        _log_error()
