#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. semi tryptic peptides for microproteins
####    2. make a reference proteome for pepquery
####    3. convert the mass calibrated mzML files to MGF for HCP and lysate
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

mkdir $PWD/proteomics/pepquery && pepquery_dir=$_

mkdir $pepquery_dir/peptide_digestion/ && digest_dir=$_

# remove stop codon asterix from microprotein sequences 
sed 's/*//g' orf_identification/orf_filtered/microproteins.fasta > \
$digest_dir/microproteins.fasta

# in silico digestion
# semi-tryptic, 2 missed cleavages, between 7 & 45 AAs
source activate oktoberfest
python scripts/mass_spectrometry/peptide_digestion.py \
$digest_dir/microproteins.fasta $digest_dir/ \
semi 2 7 45
conda deactivate

# peptide to protein mapping file
grep ',2,' $digest_dir/prosit_input.csv | cut -d, -f1  > \
$digest_dir/microprotein_peptides.txt

mv $digest_dir/prosit_input_with_proteins.csv \
$digest_dir/microprotein_peptide_mapping.csv

#### reference proteome
mkdir $pepquery_dir/known_proteome && ref_proteome=$_

cp reference_proteome/*.fasta $ref_proteome

# all mabs
cp data/protein_sequences/mabs/* $ref_proteome

# all standards
cp data/protein_sequences/standards/* $ref_proteome

find $ref_proteome/*.fasta -type f -exec cat {} \; > $ref_proteome/known_proteins.fasta

#### file conversion
# make MGF format files from mass calibrated files

type=("lysate")

for type in "${type[@]}"; do

    mkdir -p $pepquery_dir/mgf_files/$type && mgf_dir=$_

    find proteomics/metamorpheus/$type -type f -name "*calib*.mzML" > \
    $pepquery_dir/mgf/"$type"_calibrated_file_list.txt

    # source activate microprotein_process_env
    msconvert --mgf \
    -f $pepquery_dir/mgf/"$type"_calibrated_file_list.txt \
    -o $mgf_dir
    # conda deactivate

done
