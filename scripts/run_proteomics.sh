#!/bin/bash

#### Description: Execute mass spectrometry analysis pipeline
#### and quantitification 
####    1. Run a fully tryptic search of drug product and lysate data 
####    2. Run a semi tryptic search of drug product and lysate data
####    3. Preparate for pepquery analysis
####    4. Run PepQuery2 
####    5. Process the PepQuery2 results
####    6. Label free quantitation 
####    7. Quantitate microproteins in drug prodct
####    8. Differential protein abundance in lysate
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# 1. identify canonical proteins in antibody drug products 
# Full tryptic digest 

./scripts/mass_spectrometry/run_metamorpheus_ft.sh

# 2. Semi tryptic search to generate 10% FDR known peptide set 

./scripts/mass_spectrometry/run_metamorpheus_st.sh

# 3. Analyse the drug product and lysate data with PepQuer
./scripts/mass_spectrometry/prep_pepquery.sh

# 4. Analyse the drug product and lysate data with PepQuery 
./scripts/mass_spectrometry/run_pepquery.sh

# 5. Process the PepQuery output, merge PSMs with canonical 
./scripts/mass_spectrometry/run_pepquery.sh

# 6. Run FlashLFQ
./scripts/mass_spectrometry/run_flashlfq.sh

# 7. Analyse the drug product and lysate data with PepQuery
source activate microprotein_r_env
./scripts/mass_spectrometry/quantify_tzani_mps.R
conda deactivate

# 8. Analyse the drug product and lysate data with PepQuery
source activate microprotein_r_env
./scripts/mass_spectrometry/run_proda.R
conda deactivate




