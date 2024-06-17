
# remove stop codon asterix from microprotein sequences 
sed 's/*//g' orf_identification/orf_filtered/microproteins.fasta > \
$digest_dir/microproteins.fasta

# in silico digestion
# semi-tryptic, 2 missed cleavages, between 7 & 45 AAs
source activate oktoberfest
python scripts/mass_spectrometry/peptide_digestion.py \
$digest_dir/microproteins.fasta $digest_dir/ \
semi 2 7 45
conda deactivate

# peptide to protein mapping file
grep ',2,' $digest_dir/prosit_input.csv | cut -d, -f1  > \
$digest_dir/microprotein_all_peptides.txt

mv $digest_dir/prosit_input_with_proteins.csv \
$digest_dir/microprotein_peptide_mapping.csv

# we need to do two separate pepquery runs for each data type
# one is for all peptides with on ox on M as var mod
# the second is a reduced search of the peptides with the n-term of the protein
# here, we add n-term acetylation and ox on M as var mods

awk '/^>/ {if (seq) {print substr(seq, 1, 7)}; seq=""; next} {seq = seq $0} END {if (seq) print substr(seq, 1, 7)}' \
$digest_dir/microproteins.fasta | awk '{print "^" substr($0, 1, 7)}' | \
grep -f - $digest_dir/microprotein_all_peptides.txt > \
$digest_dir/microprotein_nterm_peptides.txt

# known proteins
mkdir $pepquery_dir/known_proteome && ref_proteome=$_

# uniprot and cRAP
cp reference_proteome/*.fasta $ref_proteome

# all mabs
cp data/protein_sequences/mabs/* $ref_proteome

# all standards
cp data/protein_sequences/standards/* $ref_proteome

find $ref_proteome/*.fasta -type f -exec cat {} \; > $ref_proteome/known_proteins.fasta
