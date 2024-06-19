#!/bin/bash

#### Description: Prepare figures/tables for each results section of the paper
#### Written by: NIBRT Clarke Lab. - colin.clarke@nibrt.ie

source activate microprotein_r_env


# Section 2.1
Rscript ./scripts/manuscript_prep/section2.1.R

# Section 2.2
Rscript ./scripts/prep_manuscript/section2.2.R

# Section 2.3
Rscript ./scripts/prep_manuscript/section2.3.R

# Section 2.4
Rscript ./scripts/prep_manuscript/section2.4.R

# Section 2.5
Rscript ./scripts/prep_manuscript/section2.5.R

# Section 2.6
Rscript ./scripts/prep_manuscript/section2.6.R

conda deactivate