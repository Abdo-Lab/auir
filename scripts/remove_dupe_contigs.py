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
    
    info = "Remove duplicate annotations from FASTA formatted reference file"
    
    parser = argparse.ArgumentParser(description=info)
    
    parser.add_argument('-r', '--reference_sequence', type=str, required=True,
    			help='Please provide a FASTA formatted reference file')

    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Please provide a FASTA formatted reference file')

    return parser.parse_args(cmdline_params)

def read_fasta(filename):
    """
    Removes duplicate annotations from FASTA file.
    :param (str) filename: FASTA file
    :return (dict) records: A dictionary of FASTA records
    """
    with open(filename, 'r') as fp:
        records = {}
        for line in fp:
            key = line
            value = fp.next()
            if key in records:
                if len(value) > len(records[key]):
                    records[key] = value.rstrip()
            else:
                records[key] = value
    fp.close()
    return records

def write_fasta(filename, records):
    """
    Writes a dictionary of FASTA records to file.
    :param (str) filename: Output file.
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
    records = read_fasta(opts.reference_sequence)
    write_fasta(opts.output, records)
