import os
import sys
import argparse

segment_types = {'segment_1': 'PB2', 
		 'segment_2': 'PB1', 
		 'segment_3': 'PA',
	         'segment_4': 'HA', 
		 'segment_5': 'NP', 
		 'segment_6': 'NA',
		 'segment_7': 'M1_M2', 
		 'segment_8': 'NS1_NEP', 
		 'PB2': 'PB2', 
		 'PB1': 'PB1',
		 'PA': 'PA', 
		 'HA': 'HA', 
		 'NP': 'NP', 
		 'NA': 'NA', 
		 'M1': 'M1_M2',
		 'M2': 'M1_M2',
		 'NS1': 'NS1_NEP',
		 'NEP': 'NS1_NEP'
		}

def parse_cmdline_params(cmdline_params):
    	info = "Splits contigs from a FASTA file based on their BLAST annotations"

    	parser = argparse.ArgumentParser(description=info)

    	parser.add_argument('-i', '--input_files', nargs='+', type=str, required=True,
    			    help='Please provide a FASTA formatted assembly file')

	return parser.parse_args(cmdline_params)

class Record(object):
    	def __init__(self, header, seq):
		self.header = header
		self.seq = seq

def sample_name(fasta_file):
	base = os.path.basename(fasta_file)
	sample_id = base.rsplit('.', 3)[0]
	return sample_id

def count_lines(fasta_file):
	return sum(1 for line in open(fasta_file))

def reformat_alignment_title(alignment_title):
	return alignment_title.replace(" ", "_")

def split_contigs(fasta_files):
	for fasta in fasta_files:
		records = open(fasta, 'r').read().splitlines()
		fasta_dict = {}
		line_count = count_lines(fasta)
		sample_id = sample_name(fasta)
		for i in range(0, line_count, 2):
			record = Record(records[i], records[i+1])
			for segment in segment_types:
				if segment in record.header:
					annotation = segment_types[segment]
					if annotation in fasta_dict and len(record.seq) > len(fasta_dict[annotation].seq):
						fasta_dict[annotation] = record
					if annotation not in fasta_dict:
						fasta_dict[annotation] = record
					break

		for k, v in fasta_dict.items():
			new_header = '>' + sample_id + '_' + v.header[1:len(v.header)]
			open(sample_id + '.' + k + '.fa', 'w').write(new_header + '\n' + v.seq)

if __name__ == "__main__":
	opts = parse_cmdline_params(sys.argv[1:])
	split_contigs(opts.input_files)
