from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser(description='''This script can be used to generate one vs one multimer
		combinations, using a primary input, with one amino acid sequence and a secondary input, with multiple 
		amino acid sequences.  Be warned the outfile is amended, which is useful for generating .fasta containing
		all desired combinations for input to tsl_af_runner, but unhelpful if you want separate outfiles, so remember
		to use different outfile names if required.''')
parser.add_argument("-i1", "--input1", help="File containing primary amino acid sequence.")
parser.add_argument("-i2", "--input2", help="File containing secondary amino acid sequences.")
parser.add_argument("-o", "--outfile", help="Outfile for combined sequences.")
args = parser.parse_args()


primary_aa = args.input1
secondary_aa = args.input2

def combine_aa(p_aa, s_aa):
	
	for i, s_record in enumerate(SeqIO.parse(s_aa, "fasta")):
			for p_record in SeqIO.parse(p_aa, "fasta"):
				new_id_s = SeqRecord(s_record.seq, f"{p_record.id}_{s_record.id}_{s_record.id}", "","")
				new_id_p = SeqRecord(p_record.seq, f"{new_id_s.id}", "","")
				
				with open(args.outfile, "a") as output:
					# This appends, if you use the same outfile name for other protein combinations
					# it will be added to that file
					SeqIO.write(new_id_s, output, "fasta")
					SeqIO.write(new_id_p, output, "fasta")
				
				
combine_aa(primary_aa, secondary_aa)
	