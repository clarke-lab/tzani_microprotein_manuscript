#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. semi tryptic peptides for microproteins
####    2. make a reference proteome for pepquery
####    3. convert the mass calibrated mzML files to MGF for HCP and lysate
####    4. make pepquery index
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

mkdir $PWD/proteomics/pepquery && pepquery_dir=$_

pepquery_dir=$PWD/proteomics/pepquery

mkdir $pepquery_dir/protein_digestion/ && digest_dir=$_

# remove stop codon asterix from microprotein sequences 
sed 's/*//g' orf_identification/orf_filtered/microproteins.fasta > \
$digest_dir/microproteins.fasta

# known proteins
mkdir $pepquery_dir/known_proteome && ref_proteome=$_

# uniprot and cRAP
cp reference_proteome/*.fasta $ref_proteome

# all mabs
cp data/protein_sequences/mabs/* $ref_proteome

# all standards
cp data/protein_sequences/standards/* $ref_proteome

find $ref_proteome/*.fasta -type f -exec cat {} \; > $ref_proteome/known_proteins.fasta

cp $ref_proteome/known_proteins.fasta $digest_dir

# select semi-tryptic peptides from metamorpheus FDR <= 0.1
source activate microprotein_r_env
Rscript scripts/mass_spectrometry/select_st_peptides.R
conda deactivate

source activate microprotein_process_env
# digest both microprotein and known proteins 
chainsaw -c Trypsin --numMissedCleavages 2 --minLength 7 --maxLength  45 --specificity semi $digest_dir/microproteins.fasta
chainsaw -c Trypsin -n 2 -m 7 -M 45 -s fully $digest_dir/known_proteins.fasta
conda deactivate

source activate microprotein_r_env
Rscript scripts/mass_spectrometry/peptide_digestion.R
conda deactivate

# Construct pepquery index 
# 1. make MGF format files from metamorpheus fully tryptic mzML mass calibrated files 
# 2. make an index for the drug product and lysate ms data

data_types=("lysate")
prep_types=("native" "reducing")

source activate microprotein_process_env

# specific the path the jar to allow setting -Xmx 
pepquery_jar=$CONDA_PREFIX/share/pepquery-2.0.2-0/pepquery-2.0.2.jar

for data_type in "${data_types[@]}"; do

    for prep_type in "${prep_types[@]}"; do

    if [ $prep_type == "reducing" ]; then

        mkdir -p $pepquery_dir/mgf_files/$data_type/reducing && mgf_dir=$_
                
        # list mass calibrated mzML files
        find proteomics/metamorpheus/$data_type/*/*/reducing -type f -name "*calib*.mzML" > \
        $pepquery_dir/mgf_files/"$data_type"_"$prep_type"_calibrated_file_list.txt
        
        # convert mzML to MGF
        msconvert --mgf \
        -f $pepquery_dir/mgf_files/"$data_type"_"$prep_type"_calibrated_file_list.txt \
        -o $mgf_dir

        # build the pepquery index
        mkdir -p $pepquery_dir/index/$data_type/reducing && index_dir=$_

        # make the index
        java -Xmx275G -jar $pepquery_jar index \
        -c 70 -f mgf -i $mgf_dir -o $index_dir

    elif [ "$prep_type" == "native" && "$data_type" == "drug_product" ]; then

        mkdir -p $pepquery_dir/mgf_files/$data_type/native && mgf_dir=$_
                
        # list mass calibrated mzML files
        find proteomics/metamorpheus/$data_type/*/*/$prep_type -type f -name "*calib*.mzML" > \
        $pepquery_dir/mgf_files/"$data_type"_"$prep_type"_calibrated_file_list.txt
        
        # # convert mzML to MGF
        msconvert --mgf \
        -f $pepquery_dir/mgf_files/"$data_type"_"$prep_type"_calibrated_file_list.txt \
        -o $mgf_dir

        # build the pepquery index
        mkdir -p $pepquery_dir/index/$data_type/native && index_dir=$_

        # # make the index
        java -Xmx275G -jar $pepquery_jar index \
        -c 70 -f mgf -i $mgf_dir -o $index_dir

    fi
    done
done

conda deactivate
