#!/usr/bin/env python3

"""
Purpose
-------

This module is intended parse the results of the Abricate for one or more
samples.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``abricate_files`` : Path to abricate output file.
    - e.g.: ``'abr_resfinder.tsv'``

Generated output
----------------

None


Code documentation
------------------

"""


import re
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


if __file__.endswith(".command.sh"):
    ABRICATE_FILES = '$abricate_file'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("ABRICATE_FILE: {}".format(ABRICATE_FILES))


class Abricate:
    """Main parser for Abricate output files.

    This class parses one or more output files from Abricate, usually from
    different databases. In addition to the parsing methods, it also provides
    a flexible method to filter and re-format the content of the abricate
    files.

    Parameters
    ----------
    fls : list
       List of paths to Abricate output files.
    """

    def __init__(self, fls):

        self.storage = {}
        """
        dic: Main storage of Abricate's file content. Each entry corresponds
        to a single line and contains the keys::
        
            - ``infile``: Input file of Abricate.
            - ``reference``: Reference of the query sequence.
            - ``seq_range``: Range of the query sequence in the database
             sequence.
            - ``gene``: AMR gene name.
            - ``accession``: The genomic source of the sequence.
            - ``database``: The database the sequence came from.
            - ``coverage``: Proportion of gene covered.
            - ``identity``: Proportion of exact nucleotide matches.
        """

        self._key = 0
        """
        int: Arbitrary key for unique entries in the storage attribute
        """

        self.parse_files(fls)

    def parse_files(self, fls):
        """Public method for parsing abricate output files.

        This method is called at at class instantiation for the provided
        output files. Additional abricate output files can be added using
        this method after the class instantiation.

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
        """Parser for a single abricate output file.

        This parser will scan a single Abricate output file and populate
        the :py:attr:`Abricate.storage` attribute.

        Parameters
        ----------
        fl : str
            Path to abricate output file

        Notes
        -----
        This method will populate the :py:attr:`Abricate.storage` attribute
        with all compliant lines in the abricate output file. Entries are
        inserted using an arbitrary key that is set by the
        :py:attr:`Abricate._key` attribute.

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
            Arbitrary value to compare in the left.
        op : str
            Comparison operator.
        y : int
            Arbitrary value to compare in the right.

        Returns
        -------
        x : bool
            The 'truthness' of the test.
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

    def iter_filter(self, filters, databases=None, fields=None,
                    filter_behavior="and"):
        """General purpose filter iterator.

        This general filter iterator allows the filtering of entries based
        on one or more custom filters. These filters must contain
        an entry of the `storage` attribute, a comparison operator, and the
        test value. For example, to filter out entries with coverage below 80::

            my_filter = ["coverage", ">=", 80]

        Filters should always be provide as a list of lists::

            iter_filter([["coverage", ">=", 80]])
            # or
            my_filters = [["coverage", ">=", 80],
                          ["identity", ">=", 50]]

            iter_filter(my_filters)

        As a convenience, a list of the desired databases can be directly
        specified using the `database` argument, which will only report
        entries for the specified databases::

            iter_filter(my_filters, databases=["plasmidfinder"])

        By default, this method will yield the complete entry record. However,
        the returned filters can be specified using the `fields` option::

            iter_filter(my_filters, fields=["reference", "coverage"])

        Parameters
        ----------
        filters : list
            List of lists with the custom filter. Each list should have three
            elements. (1) the key from the entry to be compared; (2) the
            comparison operator; (3) the test value. Example:
                ``[["identity", ">", 80]]``.
        databases : list
            List of databases that should be reported.
        fields : list
            List of fields from each individual entry that are yielded.
        filter_behavior : str
            options: ``'and'`` ``'or'``
            Sets the behaviour of the filters, if multiple filters have been
            provided. By default it is set to ``'and'``, which means that an
            entry has to pass all filters. It can be set to ``'or'``, in which
            case one one of the filters has to pass.

        yields
        ------
        dic : dict
            Dictionary object containing a :py:attr:`Abricate.storage` entry
            that passed the filters.

        """

        if filter_behavior not in ["and", "or"]:
            raise ValueError("Filter behavior must be either 'and' or 'or'")

        for key, dic in self.storage.items():

            # This attribute will determine whether an entry will be yielded
            # or not
            _pass = False

            # Stores the flags with the test results for each filter
            # The results will be either True or False
            flag = []

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
                    flag.append(False)
                else:
                    flag.append(True)

            # Test whether the entry will pass based on the test results
            # and the filter behaviour
            if filter_behavior == "and":
                if all(flag):
                    _pass = True
            elif filter_behavior == "or":
                if any(flag):
                    _pass = True

            if _pass:
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


class AbricateReport(Abricate):
    """Report generator for single Abricate output files

    This class is intended to parse an Abricate output file from a single
    sample and database and generates a JSON report for the report webpage.

    Parameters
    ----------
    fls : list
       List of paths to Abricate output files.
    database : (optional) str
        Name of the database for the current report. If not provided, it will
        be inferred based on the first entry of the Abricate file.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_plot_data(self):
        """ Generates the JSON report to plot the gene boxes

        Following the convention of the reports platform, this method returns
        a list of JSON/dict objects with the information about each entry in
        the abricate file. The information contained in this JSON is::

            {contig_id: <str>,
             seqRange: [<int>, <int>],
             gene: <str>,
             accession: <str>,
             coverage: <float>,
             identity: <float>
             }

        Note that the `seqRange` entry contains the position in the
        corresponding contig, not the absolute position in the whole assembly.

        Returns
        -------
        json_dic : list
            List of JSON/dict objects with the report data.
        """

        json_dic = {"plotData": {}}

        for key, entry in self.storage.items():
            # Get contig ID using the same regex as in `assembly_report.py`
            # template
            contig_id = re.search(
                ".*_NODE_([0-9]*)_.*", entry["reference"]).group(1)
            # Get database
            database = entry["database"]
            if database not in json_dic["plotData"]:
                json_dic["plotData"][database] = []

            json_dic["plotData"][database].append(
                {"contig": contig_id,
                 "seqRange": entry["seq_range"],
                 "gene": entry["gene"],
                 "accession": entry["accession"],
                 "coverage": entry["coverage"],
                 "identity": entry["identity"]}
            )

        return json_dic

    def write_report_data(self):
        """Writes the JSON report to a json file
        """

        json_dic = self.get_plot_data()

        with open(".report.json", "w") as json_report:
            json_report.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == '__main__':

    def main(abr_file):

        abr = AbricateReport(fls=[abr_file])
        abr.write_report_data()

    main(ABRICATE_FILES)
