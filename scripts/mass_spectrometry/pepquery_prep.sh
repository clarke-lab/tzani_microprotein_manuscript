#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. semi tryptic peptides for microproteins
####    2. make a reference proteome for pepquery
####    3. convert the mass calibrated mzML files to MGF for HCP and lysate
####    4. make pepquery index
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
$digest_dir/microprotein_all_peptides.txt

mv $digest_dir/prosit_input_with_proteins.csv \
$digest_dir/microprotein_peptide_mapping.csv

# we need to do two separate pepquery runs for each data type
# one is for all peptides with on ox on M as var mod
# the second is a reduced search of the peptides with the n-term of the protein
# here, we add n-term acetylation and ox on M as var mods

awk '/^>/ {if (seq) {print substr(seq, 1, 7)}; seq=""; next} {seq = seq $0} END {if (seq) print substr(seq, 1, 7)}' \
$digest_dir/microproteins.fasta | awk '{print "^" substr($0, 1, 7)}' | \
grep -f - $digest_dir/microprotein_all_peptides.txt > \
$digest_dir/microprotein_nterm_peptides.txt

# known proteins
mkdir $pepquery_dir/known_proteome && ref_proteome=$_

# uniprot and cRAP
cp reference_proteome/*.fasta $ref_proteome

# all mabs
cp data/protein_sequences/mabs/* $ref_proteome

# all standards
cp data/protein_sequences/standards/* $ref_proteome

find $ref_proteome/*.fasta -type f -exec cat {} \; > $ref_proteome/known_proteins.fasta

# Construct pepquery index 
# 1. make MGF format files from mzML mass calibrated files
# 2. make an index for the drug product and lysate ms data

data_types=("drug_product" "lysate")

pepquery_jar=$CONDA_PREFIX/share/pepquery-2.0.2-0/pepquery-2.0.2.jar

source activate microprotein_process_env

for data_type in "${data_types[@]}"; do

    mkdir -p $pepquery_dir/mgf_files/$data_type && mgf_dir=$_

    # list mass calibrated mzML files
    find proteomics/metamorpheus/$data_type -type f -name "*calib*.mzML" > \
    $pepquery_dir/mgf_files/"$data_type"_calibrated_file_list.txt
    
    # convert mzML to MGF
    msconvert --mgf \
    -f $pepquery_dir/mgf_files/"$data_type"_calibrated_file_list.txt \
    -o $mgf_dir

    # build the pepquery index
    mkdir -p $pepquery_dir/pepquery_index/$data_type && index_dir=$_

    java -Xmx275G -jar $pepquery_jar index \
    -c 70 -f mgf -i $mgf_dir -o $index_dir

done

conda deactivate
