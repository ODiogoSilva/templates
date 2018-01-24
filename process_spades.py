#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to process the output of Spades from a single sample.
The main input is an assembly file produced by spades, which will then be
filtered according to user-specified parameters.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``fastq_id``: Sample Identification string.
    - e.g.: ``'SampleA'``
- ``assembly``: Fasta file with the assembly from SPAdes.
    - e.g.: ``'contigs.fasta'``
- ``opts``: List of options for processing spades assembly.
    1. Minimum contig length.
        - e.g.: ``'150'``
    2. Minimum k-mer coverage.
        - e.g.: ``'2'``
    3. Maximum number of contigs per 1.5Mb.
        - e.g.: ``'100'``

Generated output
----------------

(Values within ``${}`` are substituted by the corresponding variable.)

- ``'${fastq_id}.assembly.fasta'`` : Fasta file with the filtered assembly.
    - e.g.: ``'Sample1.assembly.fasta'``
- ``${fastq_id}.report.fasta`` : CSV file with the results of the filters for\
    each contig.
    - e.g.: ``'Sample1.report.csv'``

Code documentation
------------------

"""

__version__ = "1.0.0"
__build__ = "16012018"
__template__ = "process_spades-nf"

import os
import json
import logging
import operator


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
    FASTQ_ID = '$fastq_id'
    ASSEMBLY_FILE = '$assembly'
    GSIZE = float('$gsize')
    OPTS = [x.strip() for x in '$opts'.strip("[]").split(",")]
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("FASTQ_ID: {}".format(FASTQ_ID))
    logger.debug("GSIZE: {}".format(GSIZE))
    logger.debug("OPTS: {}".format(OPTS))


def _log_error():
    """Nextflow specific function that logs an error upon unexpected failing
    """

    import traceback

    with open(".status", "w") as status_fh:
        logger.error("Module exited unexpectedly with error:\\n{}".format(
            traceback.format_exc()))
        status_fh.write("error")


class Assembly:
    """Class that parses and filters a Spades Fasta assembly file

    This class parses a SPAdes assembly fasta file, collects a number
    of summary statistics and metadata from the contigs, filters
    contigs based on user-defined metrics and writes filtered assemblies
    and reports.

    Parameters
    ----------
    assembly_file : str
        Path to SPAdes output assembly file.
    min_contig_len : int
        Minimum contig length when applying the initial assembly filter.
    min_kmer_cov : int
        Minimum k-mer coverage when applying the initial assembly.
        filter.
    sample_id : str
        Name of the sample for the current assembly.
    """

    def __init__(self, assembly_file, min_contig_len, min_kmer_cov,
                 sample_id):

        self.contigs = {}
        """
        dict: Dictionary storing data for each contig.
        """

        self.filtered_ids = []
        """
        list: List of filtered contig_ids.
        """

        self.min_gc = 0.05
        """
        float: Sets the minimum GC content on a contig.
        """

        self.sample = sample_id
        """
        str: The name of the sample for the assembly.
        """

        self.report = {}
        """
        dict: Will contain the filtering results for each contig.
        """

        self.filters = [
            ["length", ">=", min_contig_len],
            ["kmer_cov", ">=", min_kmer_cov]
        ]
        """
        list: Setting initial filters to check when parsing the assembly file.
        This can be later changed using the 'filter_contigs' method.
        """

        # Parse assembly and populate self.contigs
        self._parse_assembly(assembly_file)

        # Perform first contig filtering using min_contig_len, min_kmer_cov,
        # and gc content
        self.filter_contigs(*self.filters)

    def _parse_assembly(self, assembly_file):
        """Parse a Spades assembly fasta file.

        This is a Fasta parsing method that populates the
        :py:attr:`~Assembly.contigs` attribute with data for each contig in the
        assembly.

        The insertion of data on the self.contigs is done by the
        :py:meth:`Assembly._populate_contigs` method, which also calculates
        GC content and proportions.

        Parameters
        ----------
        assembly_file : str
            Path to the assembly fasta file.

        """

        # Temporary storage of sequence data
        seq_temp = []
        # Id counter for contig that will serve as key in self.contigs
        contig_id = 0
        # Initialize kmer coverage and header
        cov, header = None, None

        with open(assembly_file) as fh:

            logger.debug("Starting iteration of assembly file: {}".format(
                assembly_file))
            for line in fh:
                # Skip empty lines
                if not line.strip():
                    continue
                else:
                    # Remove whitespace surrounding line for further processing
                    line = line.strip()

                if line.startswith(">"):
                    # If a sequence has already been populated, save the
                    # previous contig information
                    if seq_temp:
                        # Use join() to convert string list into the full
                        # contig string. This is generally much more efficient
                        # than successively concatenating strings.
                        seq = "".join(seq_temp)

                        logger.debug("Populating contig with contig_id '{}', "
                                     "header '{}' and cov '{}'".format(
                                        contig_id, header, cov))
                        self._populate_contigs(contig_id, header, cov, seq)

                        # Reset temporary sequence storage
                        seq_temp = []
                        contig_id += 1

                    header = line[1:]
                    cov = float(line.split("_")[-1])

                else:
                    seq_temp.append(line)

            # Populate last contig entry
            logger.debug("Populating contig with contig_id '{}', "
                         "header '{}' and cov '{}'".format(
                            contig_id, header, cov))
            seq = "".join(seq_temp)
            self._populate_contigs(contig_id, header, cov, seq)

    def _populate_contigs(self, contig_id, header, cov, sequence):
        """ Inserts data from a single contig into\
         :py:attr:`~Assembly.contigs`.

        By providing a contig id, the original header, the coverage that
        is parsed from the header and the sequence, this method will
        populate the :py:attr:`~Assembly.contigs` attribute.

        Parameters
        ----------
        contig_id : int
            Arbitrary unique contig identifier.
        header : str
            Original header of the current contig.
        cov : float
            The contig coverage, parsed from the fasta header
        sequence : str
            The complete sequence of the contig.

        """

        # Get AT/GC/N counts and proportions.
        # Note that self._get_gc_content returns a dictionary with the
        # information on the GC/AT/N counts and proportions. This makes it
        # much easier to add to the contigs attribute using the ** notation.
        gc_kwargs = self._get_gc_content(sequence, len(sequence))
        logger.debug("Populate GC content with: {}".format(gc_kwargs))

        self.contigs[contig_id] = {
            "header": header,
            "sequence": sequence,
            "length": len(sequence),
            "kmer_cov": cov,
            **gc_kwargs
        }

    @staticmethod
    def _get_gc_content(sequence, length):
        """Get GC content and proportions.

        Parameters
        ----------
        sequence : str
            The complete sequence of the contig.
        length : int
            The length of the sequence contig.

        Returns
        -------
        x : dict
            Dictionary with the at/gc/n counts and proportions

        """

        # Get AT/GC/N counts
        at = sum(map(sequence.count, ["A", "T"]))
        gc = sum(map(sequence.count, ["G", "C"]))
        n = length - (at + gc)

        # Get AT/GC/N proportions
        at_prop = at / length
        gc_prop = gc / length
        n_prop = n / length

        return {"at": at, "gc": gc, "n": n,
                "at_prop": at_prop, "gc_prop": gc_prop, "n_prop": n_prop}

    @staticmethod
    def _test_truth(x, op, y):
        """ Test the truth of a comparisong between x and y using an \
        ``operator``.

        If you want to compare '100 > 200', this method can be called as::

            self._test_truth(100, ">", 200).

        Parameters
        ----------
        x : int
            Arbitrary value to compare in the left
        op : str
            Comparison operator
        y : int
            Arbitrary value to compare in the rigth

        Returns
        -------
        x : bool
            The 'truthness' of the test
        """

        ops = {
            ">": operator.gt,
            "<": operator.lt,
            ">=": operator.ge,
            "<=": operator.le,
        }

        return ops[op](x, y)

    def filter_contigs(self, *comparisons):
        """Filters the contigs of the assembly according to user provided\
        comparisons.

        The comparisons must be a list of three elements with the
        :py:attr:`~Assembly.contigs` key, operator and test value. For
        example, to filter contigs with a minimum length of 250, a comparison
        would be::

            self.filter_contigs(["length", ">=", 250])

        The filtered contig ids will be stored in the
        :py:attr:`~Assembly.filtered_ids` list.

        The result of the test for all contigs will be stored in the
        :py:attr:`~Assembly.report` dictionary.

        Parameters
        ----------
        comparisons : list
            List with contig key, operator and value to test.

        """

        # Reset list of filtered ids
        self.filtered_ids = []
        self.report = {}

        gc_filters = [
            ["gc_prop", ">=", self.min_gc],
            ["gc_prop", "<=", 1 - self.min_gc]
        ]

        self.filters = list(comparisons) + gc_filters

        logger.debug("Filtering contigs using filters: {}".format(
            self.filters))

        for contig_id, contig in self.contigs.items():
            for key, op, value in list(comparisons) + gc_filters:
                if not self._test_truth(contig[key], op, value):
                    self.filtered_ids.append(contig_id)
                    self.report[contig_id] = "{}/{}/{}".format(key,
                                                               contig[key],
                                                               value)
                    break
                else:
                    self.report[contig_id] = "pass"

    def get_assembly_length(self):
        """Returns the length of the assembly, without the filtered contigs.

        Returns
        -------
        x : int
            Total length of the assembly.

        """

        return sum(
            [vals["length"] for contig_id, vals in self.contigs.items()
             if contig_id not in self.filtered_ids])

    def write_assembly(self, output_file, filtered=True):
        """Writes the assembly to a new file.

        The ``filtered`` option controls whether the new assembly will be
        filtered or not.

        Parameters
        ----------
        output_file : str
            Name of the output assembly file.
        filtered : bool
            If ``True``, does not include filtered ids.
        """

        logger.debug("Writing the filtered assembly into: {}".format(
            output_file))
        with open(output_file, "w") as fh:

            for contig_id, contig in self.contigs.items():
                if contig_id not in self.filtered_ids and filtered:
                    fh.write(">{}_{}\\n{}\\n".format(self.sample,
                                                     contig["header"],
                                                     contig["sequence"]))

    def write_report(self, output_file):
        """Writes a report with the test results for the current assembly

        Parameters
        ----------
        output_file : str
            Name of the output assembly file.

        """

        logger.debug("Writing the assembly report into: {}".format(
            output_file))
        with open(output_file, "w") as fh:

            for contig_id, vals in self.report.items():
                fh.write("{}, {}\\n".format(contig_id, vals))


def main(fastq_id, assembly_file, gsize, opts):
    """Main executor of the process_spades template.

    Parameters
    ----------
    fastq_id : str
        Sample Identification string.
    assembly_file : str
        Path to the assembly file generated by Spades.
    gsize : int
        Estimate of genome size.
    opts : list
        List of options for processing spades assembly.

    """

    logger.info("Starting SPAdes processing")
    warnings = []
    fails = None

    min_contig_len, min_kmer_cov, max_contigs = [int(x) for x in opts]
    logger.debug("Setting minimum conting length to: {}".format(
        min_contig_len))
    logger.debug("Setting minimum kmer coverage: {}".format(min_kmer_cov))

    # Parse the spades assembly file and perform the first filtering.
    logger.info("Starting assembly parsing")
    spades_assembly = Assembly(assembly_file, min_contig_len, min_kmer_cov,
                               fastq_id)

    with open(".warnings", "w") as warn_fh:
        t_80 = gsize * 1000000 * 0.8
        t_150 = gsize * 1000000 * 1.5
        # Check if assembly size of the first assembly is lower than 80% of the
        # estimated genome size. If True, perform the filtering without the
        # k-mer coverage filter
        assembly_len = spades_assembly.get_assembly_length()
        logger.debug("Checking assembly length: {}".format(assembly_len))
        if assembly_len < t_80:

            logger.warning("Assembly size ({}) smaller than the minimum "
                           "threshold of 80% of expected genome size. "
                           "Applying contig filters without the k-mer "
                           "coverage filter".format(assembly_len))
            spades_assembly.filter_contigs(*[
                ["length", ">=", min_contig_len]
            ])

            assembly_len = spades_assembly.get_assembly_length()
            logger.debug("Checking updated assembly length: "
                         "{}".format(assembly_len))
            if assembly_len < t_80:

                warn_msg = "Assembly size smaller than the minimum" \
                           " threshold of 80% of expected genome size.".format(
                                assembly_len)
                logger.warning(warn_msg)
                warn_fh.write(warn_msg)
                fails = "Small_genome_size_({})".format(assembly_len)

        if assembly_len > t_150:

            warn_msg = "Assembly size ({}) smaller than the maximum" \
                       " threshold of 150% of expected genome size.".format(
                            assembly_len)
            logger.warning(warn_msg)
            warn_fh.write(warn_msg)
            fails = "Large_genome_size_({})".format(assembly_len)

        logger.debug("Checking number of contigs: {}".format(
            len(spades_assembly.contigs)))
        contig_threshold = (max_contigs * gsize) / 1.5
        if len(spades_assembly.contigs) > contig_threshold:

            warn_msg = "The number of contigs ({}) exceeds the threshold of " \
                       "100 contigs per 1.5Mb ({})".format(
                            spades_assembly.contigs, contig_threshold)

            logger.warning(warn_msg)
            warn_fh.write(warn_msg)
            warnings.append("excessive_contigs:moderate")

    # Write filtered assembly
    output_assembly = "{}.assembly.fasta".format(fastq_id)
    spades_assembly.write_assembly(output_assembly)
    # Write report
    output_report = "{}.report.csv".format(fastq_id)
    spades_assembly.write_report(output_report)
    # Write json report
    with open(".report.json", "w") as json_report:
        json_dic = {
            "warnings": {
                "process": "Spades",
                "value": warnings
            },
            "fail": {
                "process": "Spades",
                "value": fails
            }
        }
        json_report.write(json.dumps(json_dic, separators=(",", ":")))

    with open(".status", "w") as status_fh:
        status_fh.write("pass")


if __name__ == '__main__':

    try:
        build_versions()
        main(FASTQ_ID, ASSEMBLY_FILE, GSIZE, OPTS)
    except:
        _log_error()
