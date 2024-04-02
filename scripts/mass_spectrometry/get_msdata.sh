#!/bin/bash

types=("drug_product" "lysate")

for type in "${types[@]}"; do

    if [ "$type" == "drug_product" ]; then
        studies=("tzani" "pythoud")
    elif [ "$type" == "lysate" ]; then
        studies=("tzani")
    fi

    for study in "${studies[@]}"; do

        if [ "$study" == "tzani" ] && [ "$type" == "drug_product" ]; then
            samples=("adalimumab" "denosumab" "pertuzumab" "vedolizumab")
        elif [ "$study" == "pythoud" ] && [ "$type" == "drug_product" ]; then
            samples=("adalimumab" "bevacizumab" "nivolumab" "trastuzumab")
        elif [ "$study" == "tzani" ] && [ "$type" == "lysate" ]; then
            samples=("d4d7" "tempshift")
        fi
        
        echo $samples
        for sample in "${samples[@]}"; do
            
            echo "Downloading" $study $sample "analysis"

            lcms_data=proteomics/raw_files/$type/$study/$sample

            if [ ! -d "$lcms_data" ]; then
                mkdir -p $lcms_data
            fi 
            
            xargs < $PWD/data/download/msdata_dl/"$study"_"$sample"_download.txt -P 6 -L 1 wget -O

            if [ $study == "pythoud" ]; then
                find $lcms_data -name "*.zip" -exec unzip -d $lcms_data {} \;
                mv $lcms_data/*/*.raw $lcms_data
                rm $lcms_data/*.zip
                rm -r $lcms_data/*Semi-fractionation $lcms_data/*_OLD $lcms_data/*_OND $lcms_data/*_ND

             # # Loop through each line in list.txt
                 while IFS=' ' read -r current_filename new_filename; do
                
                 # Check if the current file exists in the directory
                     if [ -f "$current_filename" ]; then
                         # Rename the file to the new filename
                        mv "$current_filename" "$new_filename"
                        echo "Renamed: $current_filename -> $new_filename"

                     else
                         echo "File not found: $current_filename"
                     fi
                 done < data/download/msdata_dl/pythoud_"$sample"_rename.txt
            fi

        done
    done
done