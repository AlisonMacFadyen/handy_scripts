#!/usr/bin/env python3
"""
Ariba Results Processor
Processes Phandango CSV files from Ariba results (CARD, VFDB, PlasmidFinder)
"""

import csv
import re
import argparse
from pathlib import Path
from collections import defaultdict
import pandas as pd


class AribaProcessor:
    """
    Class for processing Ariba results from Phandango CSV files
    """
    
    def __init__(self, input_file, output_file=None):
        self.input_file = Path(input_file)
        self.output_file = output_file or self.input_file.stem + "_processed.tsv"
        self.data = []
        self.gene_columns = []
        self.all_genes = set()
        
    def extract_isolate_name(self, full_path):
        """
        Extract isolate name from full path
        e.g., './isolate_name/report.tsv' -> 'Isolate'
        """
        # Remove leading ./ if present
        # Consider update for absolute paths
        clean_path = full_path.lstrip('./')
        
        # Extract isolate name before '_' and database suffix
        # Handle patterns like: isolate_card_out/report.tsv -> isolate
        # Consider update for different naming conventions
        isolate_match = re.match(r'([^/_]+)(?:_(?:card|vfdb|plasmid))?', clean_path)
        
        if isolate_match:
            return isolate_match.group(1)
        else:
            # Fallback: use the directory name or filename stem
            parts = Path(clean_path).parts
            if len(parts) > 1:
                return parts[0].split('_')[0]
            else:
                return Path(clean_path).stem.split('_')[0]
    
    def identify_gene_columns(self, headers):
        """
        Automatically identify gene columns (we do not want the colour columns)
        """
        gene_columns = []
        
        for i, header in enumerate(headers):
            # Skip the name column (first column)
            if i == 0:
                continue
            # Skip color columns (contain ':colour' or end with ':color')
            if ':colour' in header or ':color' in header:
                continue
            
            gene_columns.append((i, header))
        
        return gene_columns
    
    def process_csv(self):
        """
        Process the Phandango CSV file
        """
        with open(self.input_file, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            headers = next(reader)
            
            # Identify gene columns automatically
            self.gene_columns = self.identify_gene_columns(headers)
            
            print(f"Found {len(self.gene_columns)} gene columns")
            print("Processing isolates...")
            
            for row in reader:
                if len(row) < len(headers):
                    continue  # Skip incomplete rows
                
                # Extract isolate name
                isolate_name = self.extract_isolate_name(row[0])
                
                # Find genes present (value = 'yes')
                genes_present = []
                for col_idx, gene_name in self.gene_columns:
                    if col_idx < len(row) and row[col_idx].lower() == 'yes':
                        genes_present.append(gene_name)
                        self.all_genes.add(gene_name)
                
                self.data.append({
                    'isolate': isolate_name,
                    'genes': genes_present
                })
        
        print(f"Processed {len(self.data)} isolates")
        print(f"Found {len(self.all_genes)} unique genes")
    
    def create_binary_matrix(self):
        """
        Create a binary matrix with isolates as rows and genes as columns
        """
        # Sort genes for consistent output
        sorted_genes = sorted(self.all_genes)
        
        # Create matrix
        matrix_data = []
        
        for entry in self.data:
            row = {'Isolate': entry['isolate']}
            for gene in sorted_genes:
                row[gene] = 1 if gene in entry['genes'] else 0
            matrix_data.append(row)
        
        return pd.DataFrame(matrix_data)
    
    def create_gene_list_format(self):
        """
        Create output in the format: Isolate\tgene_list
        """
        output_lines = []
        for entry in self.data:
            gene_list = ', '.join(sorted(entry['genes'])) if entry['genes'] else 'None'
            output_lines.append(f"{entry['isolate']}\t{gene_list}")
        
        return output_lines
    
    def write_output(self, format_type='matrix'):
        """
        Write processed data to output file
        format_type: 'matrix' or 'list'
        """
        if format_type == 'matrix':
            df = self.create_binary_matrix()
            df.to_csv(self.output_file, sep='\t', index=False)
            print(f"Binary matrix written to: {self.output_file}")
            
            # Print summary statistics
            print("\nGene frequency summary:")
            gene_counts = df.drop('Isolate', axis=1).sum().sort_values(ascending=False)
            print(gene_counts.head(10))
            
        elif format_type == 'list':
            lines = self.create_gene_list_format()
            with open(self.output_file, 'w') as outfile:
                outfile.write("Isolate\tGenes\n")
                for line in lines:
                    outfile.write(line + '\n')
            print(f"Gene list format written to: {self.output_file}")
    
    def get_gene_statistics(self):
        """
        Calculate and return gene statistics
        """
        df = self.create_binary_matrix()
        gene_counts = df.drop('Isolate', axis=1).sum().sort_values(ascending=False)
        
        stats = {
            'total_isolates': len(self.data),
            'total_genes': len(self.all_genes),
            'gene_frequencies': gene_counts.to_dict(),
            'most_common_genes': gene_counts.head(10).to_dict(),
            'isolates_without_common_genes': {}
        }
        
        # Find isolates lacking the most common genes
        for gene in gene_counts.head(5).index:
            isolates_without = df[df[gene] == 0]['Isolate'].tolist()
            stats['isolates_without_common_genes'][gene] = isolates_without
        
        return stats


def main():
    parser = argparse.ArgumentParser(description='Process Ariba Phandango CSV files')
    parser.add_argument('input_file', help='Input Phandango CSV file')
    parser.add_argument('-o', '--output', help='Output file name')
    parser.add_argument('-f', '--format', choices=['matrix', 'list'], 
                       default='matrix', help='Output format (default: matrix)')
    parser.add_argument('-s', '--stats', action='store_true', 
                       help='Print detailed statistics')
    
    args = parser.parse_args()
    
    # Create processor instance
    processor = AribaProcessor(args.input_file, args.output)
    
    # Process the CSV file
    processor.process_csv()
    
    # Write output
    processor.write_output(args.format)
    
    # Print statistics if requested
    # Consider adding functionality to write to file
    if args.stats:
        stats = processor.get_gene_statistics()
        print(f"\n=== STATISTICS ===")
        print(f"Total isolates: {stats['total_isolates']}")
        print(f"Total genes: {stats['total_genes']}")
        print(f"\nTop 10 most common genes:")
        for gene, count in stats['most_common_genes'].items():
            print(f"  {gene}: {count} ({count/stats['total_isolates']*100:.1f}%)")
        
        print(f"\nIsolates lacking common genes:")
        for gene, isolates in stats['isolates_without_common_genes'].items():
            if isolates:
                print(f"  {gene}: {', '.join(isolates)}")


if __name__ == "__main__":
    main()