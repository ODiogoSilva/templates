#!/usr/bin/env python3

"""
process_abricate template for nextflow

Purpose
-------

This module is intended parse the results of the Abricate for a given sample

Expected input
--------------
fastq_id: Pair of FastQ file paths
    .: 'SampleA'
abr_file: Path to abricate output file
    .: 'abr_resfinder.tsv

Generated output
----------------

"""


import os
import logging
import operator


if __file__.endswith(".command.sh"):
    FASTQ_ID = '$fastq_id'
    ABRICATE_FILE = '$abricate_file'


class Abricate:

    def __init__(self, fls):
        """Main parser for Abricate file

        Parameters
        ----------
        fls : list
           List of paths to Abricate output files
        """

        self.storage = {}

        self._key = 0
        """
        Arbitrary key for unique entries in the storage attribute
        """

        self.parse_files(fls)

    def parse_files(self, fls):
        """Public method for parsing abricate output files

        This method is called at class init for the provided output files.
        Additional abricate output files can be added using this method after
        the class instantiation.

        Parameters
        ----------
        fls : list
            List of paths to Abricate files

        """

        for f in fls:
            # Make sure paths exists
            if os.path.exists(f):
                self._parser(f)
            else:
                logging.warning("File {} does not exist".format(f))

    def _parser(self, fl):
        """Parser for a single abricate output file

        Parameters
        ----------
        fl : str
            Path to abricate output file

        Notes
        -----
        This method will populate the `storage` attribute with all compliant
        lines in the abricate output file. Entries are inserted using an
        arbitrary key that is set by the `_key` attribute.

        """

        with open(fl) as fh:

            for line in fh:
                # Skip header and comment lines
                if line.startswith("#") or line.strip() == "":
                    continue

                fields = line.strip().split("\t")

                try:
                    coverage = float(fields[8])
                except ValueError:
                    coverage = None
                try:
                    identity = float(fields[9])
                except ValueError:
                    identity = None

                try:
                    accession = fields[11]
                except IndexError:
                    print(line)
                    accession = None

                self.storage[self._key] = {
                    "infile": fields[0],
                    "reference": fields[1],
                    "seq_range": (int(fields[2]), int(fields[3])),
                    "gene": fields[4],
                    "accession": accession,
                    "database": fields[10],
                    "coverage": coverage,
                    "identity": identity
                }

                self._key += 1

    @staticmethod
    def _test_truth(x, op, y):
        """ Test the truth of a comparison between x and y using an operator.

        If you want to compare '100 > 200', this method can be called as
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
        _ : bool
            The truthness of the test
        """

        ops = {
            ">": operator.gt,
            "<": operator.lt,
            ">=": operator.ge,
            "<=": operator.le,
            "==": operator.eq,
            "!=": operator.ne
        }

        return ops[op](x, y)

    def iter_filter(self, filters, databases=None, fields=None):
        """General purpose filter iterator

        This general filter iterator allows the filtering of entries based
        on one or more custom filters. These filters must contain
        an entry of the `storage` attribute, a comparison operator, and the
        test value. For example, to filter out entries with coverage below 80:

            my_filter = ["coverage", ">=", 80]

        Filters should always be provide as a list of lists and using the '*'
        notation when calling the function:

            iter_filter(*[["coverage", ">=", 80]])
            # or
            my_filters = [["coverage", ">=", 80],
                          ["identity", ">=", 50]]
            iter_filter(*my_filters)

        As a convenience, a list of the desired databases can be directly
        specified using the `database` argument, which will only report
        entries for the specified databases:

            iter_filter(*my_filters, databases=["plasmidfinder"])

        By default, this method will yield the complete entry record. However,
        the returned filters can be specified using the `fields` option.

            iter_filter(*my_filters, fields=["reference", "coverage"])

        Parameters
        ----------
        filters : list
            List of lists with the custom filter. Each list should have three
            elements. (1) the key from the entry to be compared; (2) the
            comparison operator; (3) the test value. Example:
                [["identity", ">", 80]]
        databases : list
            List of databases that should be reported.
        fields : list
            List of fields from each individual entry that are yielded.

        yields
        ------
        dic : dict
            Dictionary object containing a `storage` entry that passed filters.

        """

        for key, dic in self.storage.items():

            # Turns False when a filter test fails
            flag = True

            # Filter for databases
            if databases:
                # Skip entry if not in specified database
                if dic["database"] not in databases:
                    continue

            # Apply filters
            for f in filters:
                # Get value of current filter
                val = dic[f[0]]
                if not self._test_truth(val, f[1], f[2]):
                    flag = False

            if flag:
                if fields:
                    yield dict((x, y) for x, y in dic.items() if x in fields)
                else:
                    yield dic

    def get_filter(self, *args, **kwargs):
        """ Wrapper of the iter_filter method that returns a list with results

        It should be called exactly as in the `iter_filter`

        Returns
        -------
        _ : list
            List of dictionary entries that passed the filters in the
            `iter_filter` method.

        See Also
        --------
        iter_filter
        """

        return list(self.iter_filter(*args, **kwargs))


if __name__ == '__main__':

    def main(fastq_id, abr_file):

        abr = Abricate([abr_file])

    main(FASTQ_ID, ABRICATE_FILE)
