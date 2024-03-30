#!/bin/bash

metamorpheus_dir=proteomics/metamorpheus_us

uniprot=reference_proteome/UP000001075_10029.fasta
crap=reference_proteome/crap.fasta

studies=("nibrt")

for study in "${studies[@]}"; do


    if [ "$study" == "nibrt" ]; then
        samples=("d4d7")

    elif [ "$study" == "beaumal" ]; then
        samples=("nivolumab" "trastuzumab")

    elif [ "$study" == "pythoud" ]; then
        samples=("bevacizumab")

    elif [ "$study" == "hessman" ]; then
        samples=("protein_a")

    fi

    for sample in "${samples[@]}"; do

        mkdir -p $metamorpheus_dir/$study/$sample && out_dir=$_
        
        raw_dir=proteomics/upstream_data/$study/$sample

        metamorpheus \
            -t data/metamorpheus/Task1-SearchTaskconfig.toml \
            data/metamorpheus/Task2-CalibrateTaskconfig.toml \
            data/metamorpheus/Task3-SearchTaskconfig.toml \
            -s $raw_dir \
            -d $uniprot $crap \
            -o $out_dir

    done
done