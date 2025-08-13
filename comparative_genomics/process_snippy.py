#!/usr/bin/env python3
"""
SNP Lineage Processor
Compares SNPs between evolved strains and their ancestors to identify unique mutations
in lineage evolution or within-host evolution studies
Uses output from Snippy
"""

import argparse
import csv
import pandas as pd
from pathlib import Path
from collections import defaultdict


class SNPLineageProcessor:
    """
    Class for processing and comparing SNP data between ancestor and evolved strains
    """
    
    def __init__(self, evolved_file, ancestor_file, output_file):
        self.evolved_file = Path(evolved_file)
        self.ancestor_file = Path(ancestor_file)
        self.output_file = output_file
        self.ancestor_snps = {}
        self.evolved_snps = {}
        self.unique_snps = []
        self.shared_snps = []
        self.modified_snps = []
        self.lost_snps = []
        self.reverted_snps = []
        
    def load_snp_data(self, file_path):
        """
        Load SNP data from CSV file into a dictionary
        Assumes the format: [sample, position, ...other_columns...] based on Snippy output .csv
        Returns dict with position as key and full row data as value
        """
        snp_data = {}
        
        try:
            with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
                # Try to detect if first row is header
                first_line = csvfile.readline()
                csvfile.seek(0)
                
                reader = csv.reader(csvfile)
                
                # Skip header if it exists (check if second column contains 'POS' or similar)
                first_row = next(reader)
                if any(header.upper() in ['POS', 'POSITION', 'CHROM', 'CHROMOSOME'] for header in first_row):
                    # This looks like a header, continue with next row
                    pass
                else:
                    # This is data, reset and process from beginning
                    csvfile.seek(0)
                    reader = csv.reader(csvfile)
                
                for row in reader:
                    if len(row) >= 2:  # Ensure we have at least sample and position
                        position = row[1]  # Position is in second column
                        snp_data[position] = row
                        
        except FileNotFoundError:
            print(f"Error: File {file_path} not found")
            return {}
        except Exception as e:
            print(f"Error reading {file_path}: {str(e)}")
            return {}
            
        return snp_data
    
    def compare_snps(self):
        """
        Compare SNPs between evolved and ancestor strains
        Identifies: unique SNPs, shared SNPs, modified SNPs, lost SNPs and reverted SNPs
        
        Categories:
        - unique_to_evolved: New SNPs in evolved strain
        - shared: Same SNP in both strains  
        - modified: Same position, different allele
        - lost: SNP in ancestor but completely absent from evolved
        - reverted: SNP in ancestor but evolved has reference allele (reversion) NOTE ONLY WORKS IF MAPPING TO COMPARING ANCESTOR
        """
        print("Loading ancestor SNP data...")
        self.ancestor_snps = self.load_snp_data(self.ancestor_file)
        print(f"Loaded {len(self.ancestor_snps)} ancestor SNPs")
        
        print("Loading evolved SNP data...")
        self.evolved_snps = self.load_snp_data(self.evolved_file)
        print(f"Loaded {len(self.evolved_snps)} evolved SNPs")
        
        print("Comparing SNPs...")
        
        for position, evolved_row in self.evolved_snps.items():
            if position not in self.ancestor_snps:
                # SNP unique to evolved strain
                self.unique_snps.append({
                    'position': position,
                    'type': 'unique_to_evolved',
                    'evolved_data': evolved_row,
                    'ancestor_data': None
                })
            else:
                # SNP present in both - check if it's the same
                ancestor_row = self.ancestor_snps[position]
                
                # Compare relevant columns (assuming column 4 is ALT allele)
                # Column 3 is typically REF, column 4 is typically ALT
                if len(evolved_row) > 4 and len(ancestor_row) > 4:
                    evolved_alt = evolved_row[4] if len(evolved_row) > 4 else ""
                    ancestor_alt = ancestor_row[4] if len(ancestor_row) > 4 else ""
                    evolved_ref = evolved_row[3] if len(evolved_row) > 3 else ""
                    ancestor_ref = ancestor_row[3] if len(ancestor_row) > 3 else ""
                    
                    if evolved_alt == ancestor_alt and evolved_ref == ancestor_ref:
                        # Same mutation
                        self.shared_snps.append({
                            'position': position,
                            'type': 'shared',
                            'evolved_data': evolved_row,
                            'ancestor_data': ancestor_row
                        })
                    elif evolved_alt == ancestor_ref:
                        # Evolved strain has reverted to reference allele
                        self.reverted_snps.append({
                            'position': position,
                            'type': 'reverted',
                            'evolved_data': evolved_row,
                            'ancestor_data': ancestor_row
                        })
                    else:
                        # Different mutation at same position
                        self.modified_snps.append({
                            'position': position,
                            'type': 'modified',
                            'evolved_data': evolved_row,
                            'ancestor_data': ancestor_row
                        })
        
        # Check for SNPs that are only in ancestor (lost/absent from evolved strain)
        for position, ancestor_row in self.ancestor_snps.items():
            if position not in self.evolved_snps:
                self.lost_snps.append({
                    'position': position,
                    'type': 'lost',
                    'evolved_data': None,
                    'ancestor_data': ancestor_row
                })
    
    def write_results(self, output_format='detailed'):
        """
        Write comparison results to output file
        output_format: 'detailed', 'unique_only', or 'summary'
        """
        
        if output_format == 'unique_only':
            # Write only SNPs unique to evolved strain
            with open(f"{self.output_file}.csv", 'w', newline='') as outfile:
                writer = csv.writer(outfile)
                
                # Write header
                if self.unique_snps:
                    sample_row = next((snp['evolved_data'] for snp in self.unique_snps 
                                     if snp['evolved_data'] is not None), None)
                    if sample_row:
                        writer.writerow(sample_row[0:1] + ['Position'] + sample_row[2:])
                
                # Write unique SNPs from evolved strain
                for snp in self.unique_snps:
                    if snp['type'] == 'unique_to_evolved' and snp['evolved_data']:
                        writer.writerow(snp['evolved_data'])
        
        elif output_format == 'detailed':
            # Write detailed comparison with all categories
            with open(f"{self.output_file}_detailed.csv", 'w', newline='') as outfile:
                writer = csv.writer(outfile)
                
                # Write header
                writer.writerow(['Category', 'Chromosome', 'Position', 'Type', 'Ref', 'Alt', 
                               'Evidence', 'FType', 'Strand', 'NT_POS', 'AA_POS', 'Effect', 'Locus_Tag', 'Gene', 'Product'])
                
                # Write all SNPs with categories
                all_snps = (self.unique_snps + self.shared_snps + self.modified_snps + 
                           self.lost_snps + self.reverted_snps)
                
                for snp in all_snps:
                    # For lost SNPs, use ancestor data since evolved data is None
                    row_data = snp['evolved_data'] if snp['evolved_data'] else snp['ancestor_data']
                    if row_data:
                        row = [snp['type']] + row_data
                        writer.writerow(row)
                        
                        # For modified and reverted SNPs, also write the ancestor version for comparison
                        if snp['type'] in ['modified', 'reverted'] and snp['ancestor_data']:
                            ancestor_row = [f"{snp['type']}_ancestor"] + snp['ancestor_data']
                            writer.writerow(ancestor_row)
        
        elif output_format == 'summary':
            # Write summary statistics
            with open(f"{self.output_file}_summary.txt", 'w') as outfile:
                outfile.write("SNP Comparison Summary\n")
                outfile.write("=" * 50 + "\n\n")
                outfile.write(f"Ancestor SNPs: {len(self.ancestor_snps)}\n")
                outfile.write(f"Evolved SNPs: {len(self.evolved_snps)}\n\n")
                
                outfile.write("SNP Categories:\n")
                outfile.write(f"  Unique to evolved strain: {len(self.unique_snps)}\n")
                outfile.write(f"  Lost in evolved strain: {len(self.lost_snps)}\n") 
                outfile.write(f"  Reverted to reference: {len(self.reverted_snps)}\n")
                outfile.write(f"  Shared (identical): {len(self.shared_snps)}\n")
                outfile.write(f"  Modified (same position, different allele): {len(self.modified_snps)}\n\n")
                
                # List positions for each category
                if len(self.unique_snps) > 0:
                    outfile.write("Positions unique to evolved strain:\n")
                    for snp in self.unique_snps:
                        outfile.write(f"  {snp['position']}\n")
                    outfile.write("\n")
                
                if len(self.lost_snps) > 0:
                    outfile.write("Positions lost in evolved strain:\n") 
                    for snp in self.lost_snps:
                        outfile.write(f"  {snp['position']}\n")
                    outfile.write("\n")
                    
                if len(self.reverted_snps) > 0:
                    outfile.write("Positions reverted to reference:\n")
                    for snp in self.reverted_snps:
                        outfile.write(f"  {snp['position']}\n")
                    outfile.write("\n")
                    
                if len(self.modified_snps) > 0:
                    outfile.write("Positions with modified alleles:\n")
                    for snp in self.modified_snps:
                        ancestor_alt = snp['ancestor_data'][4] if len(snp['ancestor_data']) > 4 else "?"
                        evolved_alt = snp['evolved_data'][4] if len(snp['evolved_data']) > 4 else "?"
                        outfile.write(f"  {snp['position']}: {ancestor_alt} -> {evolved_alt}\n")
    
    def get_statistics(self):
        """
        Return statistics about the SNP comparison
        """
        return {
            'ancestor_snps': len(self.ancestor_snps),
            'evolved_snps': len(self.evolved_snps),
            'unique_to_evolved': len(self.unique_snps),
            'lost_snps': len(self.lost_snps),
            'reverted_snps': len(self.reverted_snps),
            'shared_snps': len(self.shared_snps),
            'modified_snps': len(self.modified_snps)
        }
    
    def create_comparison_matrix(self):
        """
        Create a matrix showing SNP presence/absence for visualization
        """
        all_positions = set(self.ancestor_snps.keys()) | set(self.evolved_snps.keys())
        
        matrix_data = []
        for position in sorted(all_positions):
            row = {
                'Position': position,
                'Ancestor': 1 if position in self.ancestor_snps else 0,
                'Evolved': 1 if position in self.evolved_snps else 0
            }
            
            # Add classification
            if position in self.ancestor_snps and position in self.evolved_snps:
                row['Status'] = 'Shared'
            elif position in self.evolved_snps:
                row['Status'] = 'Evolved_only'
            else:
                row['Status'] = 'Ancestor_only'
            
            matrix_data.append(row)
        
        return pd.DataFrame(matrix_data)


def main():
    parser = argparse.ArgumentParser(description='Compare SNPs between evolved strain and ancestor')
    parser.add_argument("-e", "--evolved", required=True, 
                       help="SNPs of evolved strain CSV file")
    parser.add_argument("-a", "--ancestor", required=True,
                       help="SNPs of ancestor strain CSV file")
    parser.add_argument("-o", "--output", required=True,
                       help="Output file prefix")
    parser.add_argument("-f", "--format", choices=['unique_only', 'detailed', 'summary', 'all'],
                       default='unique_only', 
                       help="Output format (default: unique_only)")
    parser.add_argument("-s", "--stats", action='store_true',
                       help="Print comparison statistics")
    parser.add_argument("-m", "--matrix", action='store_true',
                       help="Create comparison matrix")
    
    args = parser.parse_args()
    
    # Create processor instance
    processor = SNPLineageProcessor(args.evolved, args.ancestor, args.output)
    
    # Perform comparison
    print("Starting SNP comparison...")
    processor.compare_snps()
    
    # Write output based on format choice
    if args.format == 'all':
        processor.write_results('unique_only')
        processor.write_results('detailed') 
        processor.write_results('summary')
    else:
        processor.write_results(args.format)
    
    # Create comparison matrix if requested
    if args.matrix:
        matrix_df = processor.create_comparison_matrix()
        matrix_df.to_csv(f"{args.output}_matrix.csv", index=False)
        print(f"Comparison matrix written to: {args.output}_matrix.csv")
    
    # Print statistics
    if args.stats or args.format == 'summary':
        stats = processor.get_statistics()
        print("\n" + "="*50)
        print("SNP COMPARISON STATISTICS")
        print("="*50)
        print(f"Ancestor SNPs: {stats['ancestor_snps']}")
        print(f"Evolved SNPs: {stats['evolved_snps']}")
        print(f"Unique to evolved: {stats['unique_to_evolved']}")
        print(f"Lost in evolved: {stats['lost_snps']}")
        print(f"Shared (identical): {stats['shared_snps']}")
        print(f"Modified (same position): {stats['modified_snps']}")
    
    print(f"\nProcessing complete! Output written to: {args.output}*")


if __name__ == "__main__":
    main()