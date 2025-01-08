import argparse
import pandas as pd
import os

# Parse the arguments
parser = argparse.ArgumentParser(description='Generate ENA manifest files')
parser.add_argument("-i", "--infile", help=".tsv file with the information for manifest file", required=True)
args = parser.parse_args()

def create_manifest(csv_file):
    # Read in the csv file, converting NaN values to an empty string
    df = pd.read_csv(csv_file, keep_default_na=False, sep = '\t')
    
    for row in df.itertuples():
        # Convert to string and check if it's non-empty before splitting
        fastq_forward = str(row.fastq_forward).strip()
        
        # Check if there are multiple fastq files for a sample
        if fastq_forward and ',' in fastq_forward:
            fastq_forward_list = [f.strip() for f in fastq_forward.split(',')]
            fastq_reverse = str(row.fastq_reverse).strip()
            fastq_reverse_list = [r.strip() for r in fastq_reverse.split(',')]
            
            for i in range(len(fastq_forward_list)):
                # Extract the required information
                study = row.study
                sample_accession = row.sample_accession
                sample_name = f"{row.sample_name}_{i + 1}"
                instrument = row.instrument
                library_source = row.library_source
                library_selection = row.library_selection
                library_strategy = row.library_strategy
                fastq_forward_file = fastq_forward_list[i]
                fastq_reverse_file = fastq_reverse_list[i]

                manifest = f"""STUDY	{study}
SAMPLE	{sample_accession}
NAME	{sample_name}
INSTRUMENT	{instrument}
LIBRARY_SOURCE	{library_source}
LIBRARY_SELECTION	{library_selection}
LIBRARY_STRATEGY {library_strategy}
FASTQ	{fastq_forward_file}
FASTQ   {fastq_reverse_file}
"""
                # Make a directory for the sample accession to save manifest file
                if not os.path.exists(sample_accession):
                    os.mkdir(sample_accession)

                # Write the manifest file
                with open(f"{sample_accession}/{sample_accession}_{sample_name}_manifest.txt", "w") as manifest_save:
                    manifest_save.write(manifest)
        else:
            # Extract the required information
            study = row.study
            sample_accession = row.sample_accession
            sample_name = row.sample_name
            instrument = row.instrument
            library_source = row.library_source
            library_selection = row.library_selection
            library_strategy = row.library_strategy
            fastq_forward = str(row.fastq_forward).strip()
            fastq_reverse = str(row.fastq_reverse).strip()

            manifest = f"""STUDY	{study}
SAMPLE	{sample_accession}
NAME	{sample_name}
INSTRUMENT	{instrument}
LIBRARY_SOURCE	{library_source}
LIBRARY_SELECTION	{library_selection}
LIBRARY_STRATEGY {library_strategy}
FASTQ	{fastq_forward}
FASTQ   {fastq_reverse}
"""
    
            # Make a directory for the sample accession to save manifest file
            if not os.path.exists(sample_accession):
                os.mkdir(sample_accession)
            # Write the manifest file
            with open(f"{sample_accession}/{sample_accession}_manifest.txt", "w") as manifest_save:
                manifest_save.write(manifest)

# Call the function
create_manifest(args.infile)