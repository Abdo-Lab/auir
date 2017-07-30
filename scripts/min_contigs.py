__author__ = "Christopher Dean"
__copyright__ = ""
__credits__ = ["Christopher Dean"]
__version__ = ""
__maintainer__ = "Christopher Dean"
__email__ = "cdean11@colostate.edu"
__status__ = "I'm doing fine."

import sys
import argparse

def parse_cmdline_params(cmdline_params):
    
    info = "Remove minimum contigs from FASTA formatted reference file"
    
    parser = argparse.ArgumentParser(description=info)
    
    parser.add_argument('-r', '--reference_sequence', type=str, required=True,
    			help='Please provide a FASTA formatted reference file')

    parser.add_argument('-m', '--min_contig', type=int, required=True,
                        help='Please provide an integer value for the min contig length')

    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Please provide an output file name in FASTA format')

    return parser.parse_args(cmdline_params)

def read_fasta(filename, min_length):
    """
    Removes minimum contigs from FASTA file
    :param (str) filename: FASTA file
    :param (int) min_length: Minimum contig length
    :return (dict) records: A dictionary of FASTA records
    """
    with open(filename, 'r') as fp:
        records = {}
        for line in fp:
            key = line
            value = fp.next()
            if len(value) > min_length:
                records[key] = value
    fp.close()
    return records

def write_fasta(filename, records):
    """
    Writes a dictionary of records to FASTA file
    :param (str) filename: Output file
    :param (dict) records: A dictionary of FASTA records
    :return (void): Void method
    """
    handler = open(filename, 'w')
    for k, v in records.items():
        handler.write(k)
        handler.write(v)
    handler.close()

if __name__ == "__main__":
    opts = parse_cmdline_params(sys.argv[1:])
    records = read_fasta(opts.reference_sequence, opts.min_contig)
    write_fasta(opts.output, records)
