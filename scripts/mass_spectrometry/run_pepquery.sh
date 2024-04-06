#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. semi tryptic peptides for microproteins
####    2. make a reference proteome for pepquery
####    3. convert the mass calibrated mzML files to MGF for HCP and lysate
####    4. make pepquery index
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# specific the path the jar to allow setting -Xmx 
pepquery_jar=$CONDA_PREFIX/share/pepquery-2.0.2-0/pepquery-2.0.2.jar

pepquery_dir=proteomics/pepquery
mkdir $pepquery_dir/pepquery_results && $pepquery_result

data_types=("drug_product" "lysate")
search_types=("all" "nterm")

for data_type in "${data_types[@]}"; do
    
    for search_type in "${search_types[@]}"; do

        mkdir $pepquery_result/$data_type/$search_type & pepquery_out=$_

        if [ "$search_type" = "all" ]; then
            modification="2"
        else
            modification="2,5"
        fi

        java -Xmx250G -jar $pepquery_jar \
        -i $pepquery_dir/peptide_digestion/microprotein_"$search_type"_peptides.txt \
        -t peptide -s 1 -c 2 \
        -varMod $modification \
        -tol 10 \
        -itol 0.05 \
        -minLength 7 \
        -maxLength 45 \
        -db $pepquery_dir/known_proteome/known_proteins.fasta \
        -ms $pepquery_dir/pepquery_index/$data_type \
        -cpu 70 \
        -o $pepquery_out \
        1> $pepquery_dir/pepquery_search_"$data_type"_"$search_type"_log.txt

    done
done
