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

# Download the uniprot and contaminant databases
if ! [ -e reference_proteome ]; then
mkdir reference_proteome
cat data/download/reference_proteome_files.txt | parallel -j 4 wget -P reference_proteome {}
gunzip reference_proteome/*.gz
fi

uniprot=reference_proteome/UP000001075_10029.fasta
crap=reference_proteome/crap.fasta

mkdir proteomics/metamorpheus
metamorpheus_dir=proteomics/metamorpheus

types=("drug_product" "lysate")

for type in "${types[@]}"; do 

    mkdir proteomics/metamorpheus/$type

    if [ $type == "drug_product" ]; then
        
        studies=("tzani" "pythoud")
                
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
                
                raw_dir=proteomics/raw_files/$type/$study/$sample
                quant_standard=data/protein_sequences/standards/hi3_standard.fasta    
                mab=data/protein_sequences/mabs/"$sample".fasta

                metamorpheus \
                -t data/metamorpheus/Task1-SearchTaskconfig.toml \
                data/metamorpheus/Task2-CalibrateTaskconfig.toml \
                data/metamorpheus/Task3-SearchTaskconfig.toml \
                -s $raw_dir \
                -d $uniprot $crap $mab $quant_standard \
                -o $metamorpheus_dir/$type/$study/$sample

            elif [ "$study" == "pythoud" ] && [ "$type" == "drug_product" ]; then

                raw_dir=proteomics/raw_files/$type/$study/$sample
                quant_standard=data/protein_sequences/standards/massprep_standard.fasta
                irt_standard=data/protein_sequences/standards/irt_standard.fasta
                mab=data/protein_sequences/mabs/"$sample".fasta

                metamorpheus \
                -t data/metamorpheus/Task1-SearchTaskconfig.toml \
                data/metamorpheus/Task2-CalibrateTaskconfig.toml \
                data/metamorpheus/Task3-SearchTaskconfig.toml \
                -s $raw_dir \
                -d $uniprot $crap $mab $quant_standard $irt_standard \
                -o $metamorpheus_dir/$type/$study/$sample

            elif [ "$study" == "tzani" ] && [ "$type" == "lysate" ]; then
                
                raw_dir=proteomics/raw_files/$type/$study/$sample

                metamorpheus \
                -t data/metamorpheus/Task1-SearchTaskconfig.toml \
                data/metamorpheus/Task2-CalibrateTaskconfig.toml \
                data/metamorpheus/Task3-SearchTaskconfig.toml \
                -s $raw_dir \
                -d $uniprot $crap \
                -o $metamorpheus_dir/$type/$study/$sample
            fi
        done
    done

done

conda deactivate 