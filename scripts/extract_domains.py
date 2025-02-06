from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help=".fasta containing sequences with domain of interest.")
parser.add_argument("-d", "--domaincoords", help="File containing the information for the domain coordinates.")
parser.add_argument("-o", "--outfile", help="Outfile to save domain sequences to.")
args = parser.parse_args()


def get_coords(coords):
    """Generate dictionary from the domain coordinate file generated from 'parse_hmmscan_results.py'"""
    
    df_coords = pd.read_csv(coords, sep='\t')
    df_coords.to_dict()
            
    return df_coords

def extract_domain(sequence, start, end):
    """Extract domain sequence given start and end coordinates (1-based)"""
    return sequence[start-1:end]

def extract_domains_from_fasta(input_file, output_file, domain_coords):
    """
    Extract domains from sequences using provided coordinates
    domain_coords should be a dict with sequence IDs as keys and (start, end) tuples as values
    """
    domain_sequences = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        if record.id in domain_coords:
            start, end = domain_coords[record.id]
            domain_seq = extract_domain(record.seq, start, end)
            
            new_record = SeqRecord(
                Seq(str(domain_seq)),
                id=record.id,
                description=f"domain_{start}_{end}"
            )
            domain_sequences.append(new_record)
    
    SeqIO.write(domain_sequences, output_file, "fasta")

extract_domains_from_fasta(args.fasta, args.outfile, get_coords(args.domaincoords))