#!/usr/bin/env python

import sys
import oktoberfest as ok

def peptide_digestion(fasta, output_dir, specificity, cleavage, min_length,max_length):
    ok.pp.digest(fasta, 
                 output_dir, 
                 "HCD", 
                 specificity, 
                 cleavage, 
                 "target", 
                 "trypsin", 
                 "KR", 
                 min_length, 
                 max_length)

fasta = sys.argv[1]
output_dir = sys.argv[2]
specificity = sys.argv[3]
cleavage = sys.argv[4]
min_length = sys.argv[5]
max_length = sys.argv[6]

peptide_digestion(fasta, output_dir, specificity, cleavage, min_length, max_length)