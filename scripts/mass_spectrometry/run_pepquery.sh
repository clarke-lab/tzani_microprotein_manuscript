#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. use pepquery to search the drug product and lysate data
####    2. search all microprotein peptides and the n-term with acetylation
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

pepquery_dir=proteomics/pepquery
mkdir $pepquery_dir/results && pepquery_result=$_

data_types=("lysate")
search_types=("nterm")

source activate microprotein_process_env

# specific the path the jar to allow setting -Xmx 
pepquery_jar=$CONDA_PREFIX/share/pepquery-2.0.2-0/pepquery-2.0.2.jar

for data_type in "${data_types[@]}"; do
    
    for search_type in "${search_types[@]}"; do

        mkdir -p $pepquery_result/$data_type/$search_type && pepquery_out=$_

        # set the variable modifications based on search type
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
        -ms $pepquery_dir/index/$data_type \
        -cpu 70 -o $pepquery_out \
        1> $pepquery_dir/pepquery_search_"$data_type"_"$search_type"_log.txt

    done
done

conda deactivate