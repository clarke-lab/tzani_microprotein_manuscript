#!/bin/bash

#### Description: preparation for the pepuery analysis of microproteins 
####    1. use pepquery to search the drug product and lysate data
####    2. search all microprotein peptides and the n-term with acetylation
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

# Stage 1 identify canonical proteins in antibody drug products 
# Full tryptic digest 

./scripts/mass_spectrometry/run_metamorpheus_ft.sh

# Stage 2 identify canonical proteins in antibody drug products 
# Semi-tryptic digest 

./scripts/mass_spectrometry/run_metamorpheus_st.sh

# Analyse the drug product and lysate data with PepQuery 
./scripts/mass_spectrometry/run_pepquery.sh





