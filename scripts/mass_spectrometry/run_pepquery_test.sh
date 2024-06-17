#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. use pepquery to search the drug product and lysate data
####    2. search all microprotein peptides and the n-term with acetylation
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

pepquery_dir=proteomics/pepquery
mkdir $pepquery_dir/results  && pepquery_result=$_

data_types=("drug_product")
search_types=("all" "nterm")
prep_type=("native" "reducing")

source activate microprotein_process_env

# specific the path the jar to allow setting -Xmx 
pepquery_jar=../bin/pepquery-2.0.5-beta/pepquery-2.0.5-beta.jar

for prep_type in "${prep_types[@]}"; do
    
    for data_type in "${data_types[@]}"; do

        for search_type in "${search_types[@]}"; do

            mkdir -p $pepquery_result/$data_type/$search_type/first_pass/$prep_type && pepquery_out=$_

            # set the variable modifications based on search type
            if [ "$search_type" = "all" ]; then
                var_modification="2"
            else
                var_modification="2,5"
            fi

            if [ "$prep_type" = "native" ]; then
                fix_modification="0"
            else
                fix_modification="1"
            fi

            
            # semi-tryptic v semi-tryptic
            java -Xmx250G -jar $pepquery_jar \
            -i data/protein_digestion/microproteins/microprotein_"$search_type"_peptides.txt \
            -t peptide -s 1 -c 1 \
            -e NoDigestion \
            -fixMod $fix_modification \
            -varMod $var_modification \
            -tol 10 \
            -itol 0.05 \
            -minLength 7 \
            -maxLength 45 \
            -db data/protein_digestion/known_proteome/known_peptides.fasta \
            -ms $pepquery_dir/index/$data_type/$prep_type \
            -cpu 70 -o $pepquery_out \
            1> $pepquery_dir/pepquery_search_"$data_type"_first_pass_"$search_type"_"$prep_type"_log.txt

            # select peptides that are significant from first pass semi-tryptic search
            # Rcode 
            Rscript ./scripts/mass_spectrometry/pepquery_first_pass_filter.R $data_type $search_type $prep_type

            # second pass pepquery for unrestricted search using tryptic digestion of
            # the known protein to mitigate the large search space
            mkdir -p $pepquery_result/$data_type/$search_type/second_pass && pepquery_out=$_ 
            java -Xmx250G -jar $pepquery_jar \
            -i proteomics/pepquery/results_test_1/$data_type/$search_type/first_pass/selected_peptides.txt \
            -t peptide -s 1 -c 1 \
            -e 1 \
            -fixMod $fix_modification \
            -varMod $var_modification \
            -tol 10 \
            -itol 0.05 \
            -minLength 7 \
            -maxLength 45 \
            -db $pepquery_dir/known_proteome/known_proteins.fasta \
            -ms $pepquery_dir/index/$data_type \
            -cpu 70 -o $pepquery_out \
            1> $pepquery_dir/second_pass_search_"$data_type"_"$search_type"_log.txt

        done
    done
done

conda deactivate