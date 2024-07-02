#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. use pepquery to search the drug product and lysate data
####    2. search all microprotein peptides and the n-term with acetylation
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

pepquery_dir=proteomics/pepquery
mkdir $pepquery_dir/results && pepquery_result=$_

pepquery_result=$pepquery_dir/results

data_types=("drug_product" "lysate")
search_types=("all" "nterm")

# specific the path the jar to allow setting -Xmx 
pepquery_jar=../bin/pepquery-2.0.5-beta/pepquery-2.0.5-beta.jar

for data_type in "${data_types[@]}"; do

    for search_type in "${search_types[@]}"; do

        if [ $data_type == "drug_product" ]; then
            prep_types=("reducing" "native") 
        else
            prep_types=("reducing") # only reducing for lysate
        fi

        for prep_type in "${prep_types[@]}"; do

            mkdir -p $pepquery_result/$data_type/$search_type/$prep_type && pepquery_out=$_

            # set the variable modifications based on search type
            if [ "$search_type" = "all" ]; then
                var_modification="2" # ox of met
            else
                var_modification="2,5" # ox of met, n-term aceylation
            fi

            if [ "$prep_type" = "native" ]; then
                fix_modification="0" # no fixed mod on C for native digestion
            else
                fix_modification="1"
            fi
            
            source activate microprotein_process_env
            
            java -Xmx250G -jar $pepquery_jar \
            -i $pepquery_dir/protein_digestion/microprotein_"$search_type"_peptides.txt \
            -t peptide -s 1 -c 2 \
            -e NoDigestion \
            -fixMod $fix_modification \
            -varMod $var_modification \
            -tol 10 \
            -itol 0.05 \
            -minLength 7 \
            -maxLength 45 \
            -db  $pepquery_dir/protein_digestion/known_peptides_"$data_type".fasta \
            -ms $pepquery_dir/index/$data_type/$prep_type \
            -cpu 70 -o $pepquery_out \
            1> $pepquery_dir/"$data_type"_"$search_type"_"$prep_type"_log.txt
            conda deactivate            
        done
    done
done
