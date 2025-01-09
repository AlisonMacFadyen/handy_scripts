import argparse
import pandas as pd
import os

# Parse the arguments
parser = argparse.ArgumentParser(description='Add fastq paths to manifest file')
parser.add_argument("-i", "--infile", help=".tsv file with the information for manifest file", required=True)
parser.add_argument("-p", "--path", help="File with fastq paths", required=True)
parser.add_argument("-o", "--outfile", help="Output file for manifest_info with paths", required=True)
args = parser.parse_args()

def add_fastq_paths(tsv_file, path_file, out_file):
    # Read in the manifest file, converting NaN values to an empty string
    df_manifest = pd.read_csv(tsv_file, keep_default_na=False, sep='\t')
    
    # Read in the paths file
    df_paths = pd.read_csv(path_file, keep_default_na=False)
    
    #Edit the extension to match that of the reads you are processing
    forward_ext = "1.fq.gz"
    reverse_ext = "2.fq.gz"
    
    # Create lists to store forward and reverse paths for each sample
    forward_paths = {}
    reverse_paths = {}
    
    # Iterate through paths to collect forward and reverse reads for each sample
    for _, path_row in df_paths.iterrows():
        mac_path = path_row['Mac_Path']
        
        # Check for each sample in the manifest
        for sample_name in df_manifest['sample_name']:
            # Note you need to edit the sample_name if it contains more than just the sample name e.g. sample_name total genomic DNA
            # For the example abvove you could use the below:

            # sample_name = sample_name.split(" ")[0]

            # Ensure sample name matches and is part of the full path
            if sample_name in mac_path:
                if forward_ext in mac_path:
                    # Add to forward paths, creating list if not exists
                    if sample_name not in forward_paths:
                        forward_paths[sample_name] = []
                    forward_paths[sample_name].append(mac_path)
                
                if reverse_ext in mac_path:
                    # Add to reverse paths, creating list if not exists
                    if sample_name not in reverse_paths:
                        reverse_paths[sample_name] = []
                    reverse_paths[sample_name].append(mac_path)
    
    # Function to convert list to comma-separated string, or empty string if no paths
    def paths_to_str(sample_name, paths_dict):
        # Note if splitting the sample name, the same needs to be applied here.
        return ','.join(paths_dict.get(sample_name, [])) if sample_name in paths_dict else ''
    
    # Add paths to the manifest dataframe
    df_manifest['fastq_forward'] = df_manifest['sample_name'].apply(
        lambda x: paths_to_str(x, forward_paths)
    )
    df_manifest['fastq_reverse'] = df_manifest['sample_name'].apply(
        lambda x: paths_to_str(x, reverse_paths)
    )
    
    # Write the updated dataframe to a new file
    df_manifest.to_csv(out_file, sep='\t', index=False)


add_fastq_paths(args.infile, args.path, args.outfile)

