__author__ = "Christopher Dean"
__copyright__ = ""
__credits__ = ["Christopher Dean"]
__version__ = ""
__maintainer__ = "Christopher Dean"
__email__ = "cdean11@colostate.edu"
__status__ = ""

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

import os
import sys
import argparse

ALIGNMENTS = 1
DATABASE = "nt"
DESCRIPTIONS = 1
FORMAT = "fasta"
HITLIST_SIZE = 1
PROGRAM = "blastn"

def parse_cmdline_params(cmdline_params):
    	info = "Runs a BLAST search to identify genomic features of interest from a FASTA file"

    	parser = argparse.ArgumentParser(description=info)

    	parser.add_argument('-i', '--input_files', nargs='+', type=str, required=True,
    			    help='Please provide a FASTA formatted reference file')

	return parser.parse_args(cmdline_params)

def sample_name(fasta_file):
	"""
    	Returns sample id of file
    	:param (str) filename: FASTA file
    	:return (str) sample_id: sample name
    	"""
	base = os.path.basename(fasta_file)
	sample_id = base.rsplit('.', 2)[0]
	return sample_id

def count_lines(fasta_file):
	"""
        Counts number of lines in a fasta file
        :param (str) filename: FASTA file
        :return (int): number of lines in file
        """
	return sum(1 for line in open(fasta_file))

def reformat_alignment_title(alignment_title):
	"""
	Replaces white space with _ character in BLAST annotation
	:param (str) alignment_title: annotation
	:return (str): reformatted annotation
	"""
	return alignment_title.replace(" ", "_")

def rename_header_lines(alignment_titles, fasta_file):
	"""
        Prefixes each FASTA header with the sample name
        :param (list) alignment_titles: list of annotations
	:param (str) FASTA file
        :return (void): void method
        """
	lines = open(fasta_file).read().splitlines()
	line_count = count_lines(fasta_file)
	sample_id = sample_name(fasta_file)
	index = 0
	for i in range(0, line_count, 2):
		lines[i] = '>' + sample_id + '_' + alignment_titles[index]
		index += 1
	open(sample_id + '.contigs.annotated.fa', 'w').write('\n'.join(lines))

def blast_search(fasta_files):
	"""
	Identifies genomic features of interest using BLAST and annotates each contig accordingly
        :param (list) fasta_files: list of fasta files
        :return (void): void method
        """
	for fasta in fasta_files:
		fasta_record = open(fasta).read()
		blast_handle = NCBIWWW.qblast(PROGRAM, DATABASE, fasta_record, descriptions = DESCRIPTIONS, alignments = ALIGNMENTS, hitlist_size = HITLIST_SIZE)
	
		alignment_titles = []

		blast_records = NCBIXML.parse(blast_handle)
        	for blast_record in blast_records:
                	for alignment in blast_record.alignments:
                       		alignment_title = reformat_alignment_title(alignment.title.encode("utf-8"))
                       		alignment_titles.append(alignment_title)
        	rename_header_lines(alignment_titles, fasta)

if __name__ == "__main__":
	opts = parse_cmdline_params(sys.argv[1:])
	blast_search(opts.input_files)
