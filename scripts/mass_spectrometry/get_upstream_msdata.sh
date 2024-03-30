
bin_dir=../bin

#study="nibrt"
studies=("roderiguez")
type=dda

echo $study
for study in "${studies[@]}"; do

    mkdir -p proteomics/upstream_data/$study/

    if [ "$study" == "nibrt" ]; then
        samples=("d4d7" "tempshift")
    elif [ "$study" == "roderiguez" ]; then
        samples=("microsome")
    elif [ "$study" == "pythoud" ]; then
        samples=("adalimumab" "bevacizumab" "nivolumab" "trastuzumab")
    elif [ "$study" == "hessman" ]; then
        samples=("protein_a")
    fi

    for sample in "${samples[@]}"; do
        
         echo "Downloading" $study $sample "analysis"

         lcms_data=proteomics/upstream_data/$study/$sample/$type

        if [ ! -d "$lcms_data" ]; then
             mkdir -p $lcms_data
         fi 
        
        if [ ! -d "$lcms_data/*mgf" ]; then
        
             xargs < $PWD/data/us_msdata_dl/"$study"_"$sample"_"$type"_download.txt -P 6 -L 1 wget -O
        
            mkdir $PWD/$lcms_data/mgf
            mkdir $PWD/$lcms_data/mzML

            mono $bin_dir/ThermoRawFileParser.exe -f 0 -d $lcms_data -o $PWD/$lcms_data/mgf \
            2> $PWD/$lcms_data/mgf_convert_error.txt \
            1> $PWD/$lcms_data/mgf_convert_log.txt 

            mono $bin_dir/ThermoRawFileParser.exe -f 1 -d $lcms_data -o $PWD/$lcms_data/mzML  \
            2> $PWD/$lcms_data/mzml_convert_error.txt \
            1> $PWD/$lcms_data/mzml_convert_log.txt 

            rm $lcms_data/*.raw
        fi
    done
    
    echo "Completed" $study $sample "download"

done
