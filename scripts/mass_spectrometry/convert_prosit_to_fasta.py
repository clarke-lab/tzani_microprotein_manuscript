#!/usr/bin/env python

import csv

def create_fasta(csv_file, fasta_file):
    # Open CSV file for reading
    with open(csv_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        
        # Open Fasta file for writing
        with open(fasta_file, 'w') as fastafile:
            header_id = 1
            for row in reader:
                # Skip empty rows
                if not row:
                    continue
                
                # Extract peptide and header
                peptide = row[0]
                header = row[3]
                
                # Write to fasta file
                fastafile.write(f'>{header}{header_id}\n{peptide}\n')
                
                # Increment header ID
                header_id += 1

# Path to your CSV file
csv_file_path = 'known_proteome/peptides.csv'

# Path to the output fasta file
fasta_file_path = 'known_proteome/known_peptides.fasta'

# Create fasta file
create_fasta(csv_file_path, fasta_file_path)

