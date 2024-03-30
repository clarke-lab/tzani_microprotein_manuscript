

source activate microprotein_process_env

mkdir proteomics/metamorpheus

metamorpheus_dir=proteomics/metamorpheus
raw_dir=proteomics/raw_files

types=("drug_product" "lysate")

for type in "${types[@]}"; do 

    mkdir proteomics/metamorpheus/$type

    if [ "$type" == "drug_product" ]; then
        studies=("tzani" "pythoud")
    elif [ "$type" == "lysate" ]; then
        studies=("tzani")
    fi

    for study in "${studies[@]}"; do

        mkdir -p proteomics/raw/$type/$study/

        for sample in "${samples[@]}"; do

            mkdir -p $metamorpheus_dir/$type/$study/$sample && out_dir=$_
            
            if [ "$study" == "tzani" ] && [ "$type" == "drug_product" ]; then
                
                raw_dir=proteomics/raw_files/$study/$sample
                quant_standard=hi3_standard.fasta    
                mab=data/protein_sequences/mabs/"$sample".fasta

                metamorpheus \
                -t data/metamorpheus/Task1-SearchTaskconfig.toml \
                data/metamorpheus/Task2-CalibrateTaskconfig.toml \
                data/metamorpheus/Task3-SearchTaskconfig.toml \
                -s $raw_dir/$type/$study/$sample \
                -d $uniprot $crap $mab $quant_standard \
                -o $metamorpheus_dir/$type/$study/$sample

            elif [ "$study" == "pythoud" ] && [ "$type" == "drug_product" ]; then

                raw_dir=proteomics/raw_files/$type/$study/$sample

                quant_standard=data/protein_sequences/standards/massprep.fasta
                irt_standard=data/protein_sequences/stanards/irt_standard.fasta
                mab=data/protein_sequences/mabs/"$sample".fasta

                metamorpheus \
                -t data/metamorpheus/Task1-SearchTaskconfig.toml \
                data/metamorpheus/Task2-CalibrateTaskconfig.toml \
                data/metamorpheus/Task3-SearchTaskconfig.toml \
                -s $raw_dir/$type/$study/$sample \
                -d $uniprot $crap $mab $quant_standard $irt_standard \
                -o $metamorpheus_dir/$type/$study/$sample

            elif [ "$study" == "tzani" ] && [ "$type" == "lysate" ]; then
                
                raw_dir=proteomics/raw_files/$type/$study/$sample

                metamorpheus \
                -t data/metamorpheus/Task1-SearchTaskconfig.toml \
                data/metamorpheus/Task2-CalibrateTaskconfig.toml \
                data/metamorpheus/Task3-SearchTaskconfig.toml \
                -s $raw_dir/$type/$study/$sample \
                -d $uniprot $crap \
                -o $metamorpheus_dir/$type/$study/$sample
            fi
        done
    done

done

conda deactivate 