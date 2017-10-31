#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to provide a summary report for a given assembly
in Fasta format.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``fastq_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``assembly`` : Path to assembly file in Fasta format.
    - e.g.: ``'assembly.fasta'``

Generated output
----------------

- ``${fastq_id}_assembly_report.csv`` : CSV with summary information of the \
    assembly.
    - e.g.: ``'SampleA_assembly_report.csv'``

"""

from collections import OrderedDict


if __file__.endswith(".command.sh"):
    FASTQ_ID = '$fastq_id'
    ASSEMBLY_FILE = '$assembly'


class Assembly:
    """Class that parses and filters an assembly file in Fasta format.

    This class parses an assembly file, collects a number
    of summary statistics and metadata from the contigs and reports.

    Parameters
    ----------
    assembly_file : str
        Path to assembly file.
    sample_id : str
        Name of the sample for the current assembly.
    """

    def __init__(self, assembly_file, sample_id):

        self.summary_info = OrderedDict([
            ("ncontigs", 0),
            ("avg_contig_size", []),
            ("n50", 0),
            ("total_len", 0),
            ("avg_gc", []),
            ("missing_data", 0)
        ])
        """
        OrderedDict: Initialize summary information dictionary. Contains keys:
        
            - ``ncontigs``: Number of contigs
            - ``avg_contig_size``: Average size of contigs
            - ``n50``: N50 metric
            - ``total_len``: Total assembly length
            - ``avg_gc``: Average GC proportion
            - ``missing_data``: Count of missing data characters
        """

        self.contigs = OrderedDict()
        """
        OrderedDict: Object that maps the contig headers to the corresponding
        sequence
        """

        self.sample = sample_id
        """
        str: Sample id
        """

        self._parse_assembly(assembly_file)

    def _parse_assembly(self, assembly_file):
        """Parse an assembly file in fasta format.

        This is a Fasta parsing method that populates the
        :py:attr:`Assembly.contigs` attribute with data for each contig in the
         assembly.

        Parameters
        ----------
        assembly_file : str
            Path to the assembly fasta file.

        """

        with open(assembly_file) as fh:

            header = None

            for line in fh:

                # Skip empty lines
                if not line.strip():
                    continue

                if line.startswith(">"):
                    # Add contig header to contig dictionary
                    header = line[1:].strip()
                    self.contigs[header] = []

                else:
                    # Add sequence string for the current contig
                    self.contigs[header].append(line.strip())

            # After populating the contigs dictionary, convert the values
            # list into a string sequence
            self.contigs = OrderedDict(
                (header, "".join(seq)) for header, seq in self.contigs.items())

    def get_summary_stats(self, output_csv):
        """Generates a CSV report with summary statistics about the assembly

        The calculated statistics are:

            - Number of contigs
            - Average contig size
            - N50
            - Total assembly length
            - Average GC content
            - Amount of missing data

        Parameters
        ----------
        output_csv: str
            Name of the output CSV file.
        """

        contig_size_list = []

        self.summary_info["ncontigs"] = len(self.contigs)

        for contig_id, sequence in self.contigs.items():

            # Get contig sequence size
            contig_len = len(sequence)

            # Add size for average contig size
            contig_size_list.append(contig_len)

            # Add to total assembly length
            self.summary_info["total_len"] += contig_len

            # Add to average gc
            self.summary_info["avg_gc"].append(
                sum(map(sequence.count, ["G", "C"])) / contig_len
            )

            # Add to missing data
            self.summary_info["missing_data"] += sequence.count("N")

        # Get average contig size
        self.summary_info["avg_contig_size"] = \
            sum(contig_size_list) / len(contig_size_list)

        # Get average gc content
        self.summary_info["avg_gc"] = \
            sum(self.summary_info["avg_gc"]) / len(self.summary_info["avg_gc"])

        # Get N50
        cum_size = 0
        for l in sorted(contig_size_list, reverse=True):
            cum_size += l
            if cum_size >= self.summary_info["total_len"] / 2:
                self.summary_info["n50"] = l
                break

        # Write summary info to CSV
        with open(output_csv, "w") as fh:
            summary_line = "{}, {}\\n".format(
                self.sample, ",".join(
                    [str(x) for x in self.summary_info.values()]))
            fh.write(summary_line)


def main(fastq_id, assembly_file):
    """Main executor of the assembly_report template.

    Parameters
    ----------
    fastq_id : str
        Sample Identification string.
    assembly_file : str
        Path to assembly file in Fasta format

    """

    assembly_obj = Assembly(assembly_file, fastq_id)

    assembly_obj.get_summary_stats("{}_assembly_report.csv".format(fastq_id))


if __name__ == '__main__':
    main(FASTQ_ID, ASSEMBLY_FILE)
