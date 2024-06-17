#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. use pepquery to search the drug product and lysate data
####    2. search all microprotein peptides and the n-term with acetylation
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

pepquery_dir=proteomics/pepquery
mkdir $pepquery_dir/results  && pepquery_result=$_

pepquery_result=$pepquery_dir/results

data_types=("drug_product")
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

            mkdir -p $pepquery_result/$data_type/$search_type/first_pass/$prep_type && pepquery_out=$_

            # set the variable modifications based on search type
            if [ "$search_type" = "all" ]; then
                var_modification="2"
            else
                var_modification="2,5" # n-term aceylation
            fi

            if [ "$prep_type" = "native" ]; then
                fix_modification="0" # no fixed mod on C for native digestion
            else
                fix_modification="1"
            fi
            
            source activate microprotein_process_env
            # semi-tryptic v semi-tryptic
            java -Xmx250G -jar $pepquery_jar \
            -i $pepquery_dir/protein_digestion/microprotein_"$search_type"_peptides.txt \
            -t peptide -s 1 -c 1 \
            -e NoDigestion \
            -fixMod $fix_modification \
            -varMod $var_modification \
            -tol 10 \
            -itol 0.05 \
            -minLength 7 \
            -maxLength 35 \
            -db  $pepquery_dir/protein_digestion/known_peptides.fasta \
            -ms $pepquery_dir/index/$data_type/$prep_type \
            -cpu 70 -o $pepquery_out \
            1> $pepquery_dir/"$data_type"_first_pass_"$search_type"_"$prep_type"_log.txt
            conda deactivate

            # select peptides that are significant from first pass semi-tryptic search
            # Rcode 

            source activate microprotein_r_env            
            Rscript ./scripts/mass_spectrometry/pepquery_first_pass_filter.R $data_type $search_type $prep_type
            conda deactivate

            # second pass pepquery for unrestricted search using tryptic digestion of
            # the known protein to mitigate the large search space
            mkdir -p $pepquery_result/$data_type/$search_type/second_pass/$prep_type && pepquery_out=$_

            first_pass_peptides=$pepquery_dir/results/$data_type/$search_type/first_pass/$prep_type/selected_peptides.txt
            
            source activate microprotein_process_env
            if [[ $(wc -l < $first_pass_peptides) -ge 2 ]]; then
                
                java -Xmx250G -jar $pepquery_jar \
                -i $first_pass_peptides \
                -t peptide -s 1 -c 1 \
                -e 1 \
                -fixMod $fix_modification \
                -varMod $var_modification \
                -tol 10 \
                -itol 0.05 \
                -minLength 7 \
                -maxLength 35 \
                -db $pepquery_dir/known_proteome/known_proteins.fasta \
                -ms $pepquery_dir/index/$data_type/$prep_type \
                -cpu 70 -o $pepquery_out \
                1> $pepquery_dir/"$data_type"_second_pass_"$search_type"_"$prep_type"_log.txt
            fi
            conda deactivate
        done
    done
done