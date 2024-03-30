#!/bin/bash



studies=("nibrt")

for study in "${studies[@]}"; do

    mkdir -p proteomics/downstream_data/$study/

    if [ "$study" == "tzani" ]; then
        samples=("adalimumab" "denosumab" "pertuzumab" "vedolizumab")
    elif [ "$study" == "pythoud" ]; then
        samples=("adalimumab" "bevacizumab" "nivolumab" "trastuzumab")
    fi

    for sample in "${samples[@]}"; do
        
        echo "Downloading" $study $sample "analysis"

        lcms_data=proteomics/downstream_data/$study/$sample

        if [ ! -d "$lcms_data" ]; then
            mkdir -p $lcms_data
        fi 
        
        xargs < $PWD/data/ds_msdata_dl/"$study"_"$sample"_download.txt -P 6 -L 1 wget -O

        if [ $study == "pythoud" ]; then
            find $lcms_data -name "*.zip" -exec unzip -d $lcms_data {} \;
            mv $lcms_data/*/*.raw $lcms_data
            rm $lcms_data/*.zip
            rm -r $lcms_data/*Semi-fractionation $lcms_data/*_OLD $lcms_data/*_OND $lcms_data/*_ND

        # # Loop through each line in list.txt
            while IFS=' ' read -r current_filename new_filename; do
            # Check if the current file exists in the directory
                if [ -f $lcms_data/"$current_filename" ]; then
                    # Rename the file to the new filename
                    mv $lcms_data/"$current_filename" $lcms_data/"$new_filename"
                    echo "Renamed: $current_filename -> $new_filename"
                else
                    echo "File not found: $current_filename"
                fi
            done < data/ds_msdata_dl/pythoud_"$sample"_rename.txt
        fi

        mkdir $PWD/$lcms_data/mgf
        mkdir $PWD/$lcms_data/mzML

        mono ../bin/ThermoRawFileParser.exe -f 0 -d $lcms_data -o $PWD/$lcms_data/mgf \
        2> $PWD/$lcms_data/mgf_convert_error.txt \
        1> $PWD/$lcms_data/mgf_convert_log.txt 

        mono ../bin/ThermoRawFileParser.exe -f 1 -d $lcms_data -o $PWD/$lcms_data/mzML  \
        2> $PWD/$lcms_data/mzml_convert_error.txt \
        1> $PWD/$lcms_data/mzml_convert_log.txt 

        rm $lcms_data/*.raw

    echo "Completed" $study $sample "analysis"

     done
done