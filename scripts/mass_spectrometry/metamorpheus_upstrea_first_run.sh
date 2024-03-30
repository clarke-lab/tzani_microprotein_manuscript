

metamorpheus_dir=proteomics/metamorpheus_ds


find proteomics/metamorpheus_us/*/*/Task2CalibrationTask/ -type f -name "*calib*.mzML" > \
$pepquery_dir/lysate_ms_file_list.txt

pepquery_mgf=proteomics/ds_pepquery_analysis/mgf_lysate
mkdir $pepquery_dir/lysate_index && pepquery_index=$_
pepquery_jar=/mnt/HDD2/bin/miniconda/envs/microprotein_process_env/share/pepquery-2.0.2-0/pepquery-2.0.2.jar
java -Xmx275G -jar $pepquery_jar \
index \
-c 70 \
-i $pepquery_mgf \
-f mgf \
-o $pepquery_index

mkdir $pepquery_dir/lysate_output && pepquery_out=$_
nohup java -Xmx250G -jar $pepquery_jar \
        -i  $pepquery_dir/microprotein_peptides.txt \
        -t peptide \
        -varMod 2 \
        -s 1 \
        -c 2 \
        -tol 10 \
        -itol 0.05 \
        -db $pepquery_dir/reference_proteome.fasta \
        -ms $pepquery_index \
        -minLength 7 \
        -maxLength 45 \
        -cpu 70 \
        -o $pepquery_out \
        1> $pepquery_dir/pepquery_lysate_log.txt &

mkdir $pepquery_dir/lysate_start_output && pepquery_out=$_
        nohup java -Xmx250G -jar $pepquery_jar \
                -i  $pepquery_dir/microprotein_start_peptides.txt \
                -t peptide \
                -varMod 2,5 \
                -s 1 \
                -c 2 \
                -tol 10 \
                -itol 0.05 \
                -db $pepquery_dir/reference_proteome.fasta \
                -ms $pepquery_index \
                -minLength 7 \
                -maxLength 45 \
                -cpu 70 \
                -o $pepquery_out \
                1> $pepquery_dir/pepquery_lysate_start_log.txt &
