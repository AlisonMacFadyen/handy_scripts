# Script to process hmmscan output

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--domtblout", help="hmmscan parseable table of per-domain hits")
parser.add_argument("-p", "--pfamdomain", help="Specific Pfam domain you want to extract coordinates for")
parser.add_argument("-o", "--outfile", help="Outfile to save results to.")
args = parser.parse_args()

def parse_hmmscan_output(domtblout_file):
    domain_coords = {}
    with open(domtblout_file) as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines
                continue
            fields = line.split()
            protein_name = fields[3]
            domain_name = fields[1]
            env_start = int(fields[19])  # Envelope start
            env_end = int(fields[20])    # Envelope end
            
            # Store coordinates (you might want to add domain name to the key if proteins have multiple domains)
            if protein_name not in domain_coords:
                domain_coords[protein_name] = []
            domain_coords[protein_name].append((env_start, env_end, domain_name))
            
    df_all = pd.DataFrame(data = domain_coords)
    df_all.to_csv(args.outfile, sep = '\t', index = False)

def parse_hmmscan_output_domain(domtblout_file, domain_of_interest):
    domain_coords = {}
    with open(domtblout_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            protein_name = fields[3]
            domain_name = fields[1]
            
            if domain_name == domain_of_interest:
                env_start = int(fields[19])
                env_end = int(fields[20])
                if protein_name not in domain_coords:
                    domain_coords[protein_name] = (env_start, env_end)
    
    df_domains = pd.DataFrame(data = domain_coords)
    df_domains.to_csv(args.outfile, sep = '\t', index = False)


if args.pfamdomain:
    parse_hmmscan_output_domain(args.domtblout, args.pfamdomain)

else:
    parse_hmmscan_output(args.domtblout)

