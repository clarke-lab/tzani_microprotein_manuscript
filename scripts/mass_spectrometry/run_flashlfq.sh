#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. use pepquery to search the drug product and lysate data
####    2. search all microprotein peptides and the n-term with acetylation
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

source activate microprotein_process_env

flashlfq_dir=proteomics/flashlfq
metamorpheus_dir=proteomics/metamorpheus

types=("drug_product" "lysate")

for type in "${types[@]}"; do 

    

    if [ $type == "drug_product" ]; then
            
            studies=("tzani")
                    
    elif [ "$type" == "lysate" ]; then
            
            studies=("tzani")
            
    fi

    for study in "${studies[@]}"; do
       
        if [ $study == "tzani" ] && [ $type == "drug_product" ]; then

            samples=("adalimumab" "denosumab" "etanercept" "nivolumab" "pertuzumab" "vedolizumab")


        elif [ "$study" == "tzani" ] && [ "$type" == "lysate" ]; then

            samples=("d4d7" "tempshift")

        fi

        for sample in "${samples[@]}"; do

        mkdir -p proteomics/flashlfq/result/$type/$sample && out_dir=$_
    
        grep -v '|' $flashlfq_dir/psm/$type/"$study"_"$sample".psms.tsv > \
        $flashlfq_dir/psm/$type/"$study"_"$sample".psms.filtered.tsv

        FlashLFQ \
        --idt $flashlfq_dir/psm/$type/"$study"_"$sample".psms.filtered.tsv \
        --rep $metamorpheus_dir/$type/$study/$sample/Task2CalibrationTask \
        --ppm 10 \
        --thr 40 \
        --mbr true \
        --mrt 0.7 \
        --out $out_dir
        done
    done

done

conda deactivate