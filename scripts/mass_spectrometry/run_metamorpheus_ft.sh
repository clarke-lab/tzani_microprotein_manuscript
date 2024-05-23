#!/bin/bash
#### Description: Make coverage tracks for the manuscript figures
####             1. Merged Data for Figure 2B
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie 

### 1 merged transcriptome tracks for new annotations 

# # 1. Merged Data for Figure 2B
# ## Demonstrate the coverage at transcript-level
# ## a) CHX full coverage
# ## b) CHX p-site offset
# ## c) CHX a-site offset
# ## d) HARR-ND p-site

source activate metamorpheus

# Download the uniprot and cRAP contaminant fasta
if ! [ -e reference_proteome ]; then
mkdir reference_proteome
cat data/download/reference_proteome_files.txt | parallel -j 4 wget -P reference_proteome {}
gunzip reference_proteome/*.gz
fi

uniprot=reference_proteome/UP000001075_10029.fasta
crap=reference_proteome/crap.fasta

mkdir proteomics/metamorpheus
metamorpheus_dir=proteomics/metamorpheus

types=("drug_product")

for type in "${types[@]}"; do 

    mkdir proteomics/metamorpheus/$type

    if [ $type == "drug_product" ]; then
        
        studies=("pythoud")
                
    elif [ "$type" == "lysate" ]; then
        
        studies=("tzani")
        
    fi
    
    for study in "${studies[@]}"; do
       
        mkdir -p proteomics/raw/$type/$study/

        if [ $study == "tzani" ] && [ $type == "drug_product" ]; then

            samples=("adalimumab" "denosumab" "etanercept" "nivolumab" "pertuzumab" "vedolizumab")

        elif [ "$study" == "pythoud" ] && [ "$type" == "drug_product" ]; then

            samples=("adalimumab" "bevacizumab" "nivolumab" "trastuzumab")

        elif [ "$study" == "tzani" ] && [ "$type" == "lysate" ]; then

            samples=("d4d7" "tempshift")

        fi

        for sample in "${samples[@]}"; do

            mkdir -p $metamorpheus_dir/$type/$study/$sample && out_dir=$_
            
            if [ $study == "tzani" ] && [ $type == "drug_product" ]; then
                
                mkdir -p $metamorpheus_dir/$type/$study/$sample/reducing 

                raw_dir=proteomics/raw_files/$type/$study/$sample
                quant_standard=data/protein_sequences/standards/hi3_standard.fasta    
                mab=data/protein_sequences/mabs/"$sample".fasta

                metamorpheus \
                -t data/metamorpheus/Task1-SearchTaskconfig_reducing.toml \
                data/metamorpheus/Task2-CalibrateTaskconfig_reducing.toml \
                data/metamorpheus/Task3-SearchTaskconfig_reducing.toml \
                -s $raw_dir \
                -d $uniprot $crap $mab $quant_standard \
                -o $metamorpheus_dir/$type/$study/$sample/reducing

            elif [ "$study" == "pythoud" ] && [ "$type" == "drug_product" ]; then

                prep_types=("native" "reducing")

                for prep_type in "${prep_types[@]}"; do

                mkdir -p $metamorpheus_dir/$type/$study/$sample/$prep_type

                raw_dir=proteomics/raw_files/$type/$study/$sample/$prep_type
                quant_standard=data/protein_sequences/standards/massprep_standard.fasta
                irt_standard=data/protein_sequences/standards/irt_standard.fasta
                mab=data/protein_sequences/mabs/"$sample".fasta

                metamorpheus \
                -t data/metamorpheus/Task1-SearchTaskconfig_"$prep_type".toml \
                data/metamorpheus/Task2-CalibrateTaskconfig_"$prep_type".toml \
                data/metamorpheus/Task3-SearchTaskconfig_"$prep_type".toml \
                -s $raw_dir \
                -d $uniprot $crap $mab $quant_standard $irt_standard \
                -o $metamorpheus_dir/$type/$study/$sample/$prep_type
                
                done

            elif [ "$study" == "tzani" ] && [ "$type" == "lysate" ]; then
                
                raw_dir=proteomics/raw_files/$type/$study/$sample/reducing

                metamorpheus \
                -t data/metamorpheus/Task1-SearchTaskconfig_reducing.toml \
                data/metamorpheus/Task2-CalibrateTaskconfig_reducing.toml \
                data/metamorpheus/Task3-SearchTaskconfig_reducing.toml \
                -s $raw_dir \
                -d $uniprot $crap \
                -o $metamorpheus_dir/$type/$study/$sample/reducing
            fi
        done
    done

done

conda deactivate 