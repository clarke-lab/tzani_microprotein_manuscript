# import the pepquery microprotein search results
import_pepquery <- function(pepquery_path, sample_type, region, pass, sample_prep,
                            pep_to_prot_map) {
  
  if (region == "nterm"){
  
    psms <- read_delim(paste0(pepquery_path, "results/", sample_type, "/", 
                              region, "/", pass, "/", sample_prep, "/psm_rank.txt"), 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE, show_col=F) %>%
      filter(str_detect(modification, "Acetylation")) %>% #only acetylated nterm
      filter(confident == "Yes") %>% 
      left_join(pep_to_prot_map, by="peptide")
    
  } else {
    
    psms <- read_delim(paste0(pepquery_path, "results/", sample_type, "/", 
                              region, "/", pass, "/", sample_prep, "/psm_rank.txt"),
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE, show_col=F) %>%
      filter(confident == "Yes") %>%
      left_join(pep_to_prot_map, by="peptide") 
    }
  
  if (sample_type == "drug_product"){

    # assign study, drug product, and replicate
   psms <- psms  %>% # only Ox on M variable mod
     mutate(study = case_when(
       str_detect(spectrum_title,"tzani") ~ "tzani",
       str_detect(spectrum_title,"pythoud") ~ "pythoud"
     )) %>%
     mutate(product = case_when(
       str_detect(spectrum_title,"adalimumab") ~ "adalimumab",
       str_detect(spectrum_title,"denosumab") ~ "denosumab",
       str_detect(spectrum_title,"pertuzumab") ~ "pertuzumab",
       str_detect(spectrum_title,"vedolizumab") ~ "vedolizumab",
       str_detect(spectrum_title,"trastuzumab") ~ "trastuzumab",
       str_detect(spectrum_title,"nivolumab") ~ "nivolumab",
       str_detect(spectrum_title,"etanercept") ~ "etanercept",
       str_detect(spectrum_title,"bevacizumab") ~ "bevacizumab",
       str_detect(spectrum_title,"trastuzumab") ~ "trastuzumab")) %>%
     mutate(rep=case_when(
       study == "tzani" & str_detect(spectrum_title,"b1_t1") ~ "r1",
       study == "tzani" & str_detect(spectrum_title,"b1_t2") ~ "r2",
       study == "tzani" & str_detect(spectrum_title,"b1_t3") ~ "r3",
       study == "tzani" & str_detect(spectrum_title,"b2_t1") ~ "r4",
       study == "tzani" & str_detect(spectrum_title,"b2_t2") ~ "r5",
       study == "tzani" & str_detect(spectrum_title,"b2_t3") ~ "r6",
       study == "pythoud" & str_detect(spectrum_title,"r1") ~ "r1",
       study == "pythoud" & str_detect(spectrum_title,"r2") ~ "r2",
       study == "pythoud" & str_detect(spectrum_title,"r3") ~ "r3",
       study == "pythoud" & str_detect(spectrum_title,"r4") ~ "r4",
       study == "pythoud" & str_detect(spectrum_title,"r5") ~ "r5"
     )) %>%
     mutate(rep=case_when(
       study == "tzani" & str_detect(spectrum_title,"b1_t1") ~ "r1",
       study == "tzani" & str_detect(spectrum_title,"b1_t2") ~ "r2",
       study == "tzani" & str_detect(spectrum_title,"b1_t3") ~ "r3",
       study == "tzani" & str_detect(spectrum_title,"b2_t1") ~ "r4",
       study == "tzani" & str_detect(spectrum_title,"b2_t2") ~ "r5",
       study == "tzani" & str_detect(spectrum_title,"b2_t3") ~ "r6",
       study == "pythoud" & str_detect(spectrum_title,"r1") ~ "r1",
       study == "pythoud" & str_detect(spectrum_title,"r2") ~ "r2",
       study == "pythoud" & str_detect(spectrum_title,"r3") ~ "r3",
       study == "pythoud" & str_detect(spectrum_title,"r4") ~ "r4",
       study == "pythoud" & str_detect(spectrum_title,"r5") ~ "r5"
     )) %>%
     mutate(sample_prep = case_when(
       study == "tzani" ~ "nibrt_method",
       study == "pythoud" & str_detect(spectrum_title,"ond_r") ~ "ond",
       study == "pythoud" & str_detect(spectrum_title,"old_r") ~ "old",
       study == "pythoud" & str_detect(spectrum_title,"nd_r") ~ "nd",
       study == "pythoud" & str_detect(spectrum_title,"sfp1") ~ "sfp1",
       study == "pythoud" & str_detect(spectrum_title,"sfp2") ~ "sfp2",
       study == "pythoud" & str_detect(spectrum_title,"sfp3") ~ "sfp3"
     )) 
     

  } else {
    # assign study, drug product, and replicate
    psms <- psms %>%
      mutate(experiment = case_when(
        str_detect(spectrum_title,"ts_") ~ "tempshift",
        str_detect(spectrum_title,"^d") ~ "d4d7"
      )) %>%
      mutate(condition=case_when(
        str_detect(spectrum_title,"nts_24hr") ~ "nts_24hr",
        str_detect(spectrum_title,"ts_24hr") ~ "ts_24hr",
        str_detect(spectrum_title,"ts_48hr") ~ "ts_48hr",
        str_detect(spectrum_title,"d4") ~ "day4",
        str_detect(spectrum_title,"d7") ~ "day7")) %>%
      mutate(rep=case_when(
        experiment == "tempshift" & str_detect(spectrum_title,"b1_t1") ~ "r1",
        experiment == "tempshift" & str_detect(spectrum_title,"b1_t2") ~ "r2",
        experiment == "tempshift" & str_detect(spectrum_title,"b2_t1") ~ "r3",
        experiment == "tempshift" & str_detect(spectrum_title,"b2_t2") ~ "r4",
        experiment == "tempshift" & str_detect(spectrum_title,"b3_t1") ~ "r5",
        experiment == "tempshift" & str_detect(spectrum_title,"b3_t2") ~ "r6",
        experiment == "d4d7" & str_detect(spectrum_title,"b1") ~ "r1",
        experiment == "d4d7" & str_detect(spectrum_title,"b2") ~ "r2",
        experiment == "d4d7" & str_detect(spectrum_title,"b3") ~ "r3",
        experiment == "d4d7" & str_detect(spectrum_title,"b4") ~ "r4"
      )) 
    
  }
  return(psms)
}

# import the metamorpheus canonical protein search results
import_canonical <- function(type, study, sample, prep) {
  
  proteins <- read_delim(paste0("proteomics/metamorpheus/", type, "/",study,"/",sample, "/", prep,
                                "/Task3SearchTask/AllQuantifiedProteinGroups.tsv"), 
                         delim = "\t", escape_double = FALSE, show_col=F,
                         col_select = c(`Protein Accession`, `Gene`, 
                                        `Protein Full Name`, `Organism`, 
                                        `Unique Peptides`, `Shared Peptides`,
                                        `Sequence Coverage Fraction`, 
                                        `Number of PSMs`,
                                        `Protein QValue`,
                                        `Number of Peptides`),
                         trim_ws = TRUE) %>% 
    filter(`Number of Peptides` > 1 & `Protein QValue` < 0.01) %>%
    filter(!str_detect(replace_na(`Protein Full Name`, ''), "Keratin")) 
    
    if (type=="drug_product"){ 
       proteins <- proteins %>%
       mutate(study=study, product=sample) %>% 
       dplyr::select(study, product,`Protein Accession`, `Gene`, 
                  `Protein Full Name`, `Organism`, `Number of PSMs`,
                  `Number of Peptides`,`Sequence Coverage Fraction`, 
                  `Protein QValue`) 
    } else {
        proteins <- proteins %>% 
        mutate(study=study, experiment=sample) %>% 
        dplyr::select(study, experiment,`Protein Accession`, `Gene`, 
                  `Protein Full Name`, `Organism`, `Number of PSMs`,
                  `Number of Peptides`,`Sequence Coverage Fraction`, 
                  `Protein QValue`) 
    }
   
  allpsms <- read_delim(paste0("proteomics/metamorpheus/", type, "/", study,"/", sample, "/", prep,
                               "/Task3SearchTask/AllPSMs.psmtsv"), 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE, show_col_types = FALSE) %>%
    filter(`Protein Accession` %in% proteins$`Protein Accession`) %>%
    filter(`PEP_QValue` < 0.01)

    if (type=="drug_product"){  
    allpsms <- allpsms %>%
    mutate(study=study, product=sample) %>%
    dplyr::select(study,product,`File Name`, `Base Sequence`, `Full Sequence`,
                  `Peptide Monoisotopic Mass`, `Scan Retention Time`, 
                  `Precursor Charge`, `Protein Accession`) %>%
    mutate(`Peptide Monoisotopic Mass` = as.numeric(`Peptide Monoisotopic Mass`))
    } else {
    allpsms <- allpsms %>%
    mutate(study=study, experiment=sample) %>%
    dplyr::select(study,experiment,`File Name`, `Base Sequence`, `Full Sequence`,
                  `Peptide Monoisotopic Mass`, `Scan Retention Time`, 
                  `Precursor Charge`, `Protein Accession`) %>%
    mutate(`Peptide Monoisotopic Mass` = as.numeric(`Peptide Monoisotopic Mass`))    
    }
  
  output <- list(proteins=proteins, psms=allpsms)
  
  return(output)
}

# add the modification strings compatiable with flashlfq for pepquery results
add_peptide_mods <- function(mod_string, peptide) {
  if (mod_string != "-") {
    
    # Split the pepquery mod string
    mods <- unlist(strsplit(mod_string, ";"))
    
    # sort the substrings so largest first, to keep the start ref for stri_sub
    # below
    mods <- mods[order(as.numeric(gsub("@", "", str_extract(mods, "@[0-9]+"))), 
                       decreasing = T)]
    
    foreach(mod = mods) %do% {
      
      if (str_detect(mod, "Acetylation")) {
        stri_sub(peptide, 0, 0) <- "[Common Biological:Acetylation on X]"
      } else {
        # regex
        pattern <- "* of ([A-Za-z])@(\\d+)"
        match <- regexec(pattern, mod)
        modification_info <- list()
        modification_type <- regmatches(mod, match)[[1]][2]
        modification_number <- as.numeric(regmatches(mod, match)[[1]][3])
        
        if (modification_type == "C") {
          stri_sub(peptide, modification_number, modification_number) <- "C[Common Fixed:Carbamidomethyl on C]"
        } else if (modification_type == "M") {
          stri_sub(peptide, modification_number, modification_number) <- "M[Common Variable:Oxidation on M]"
        }
      }
      # extact the mod ino
    }
    
    # Return the modified peptide
    return(peptide)
    
  } else {
    
    # Return the unmodified peptide
    return(peptide)
  }
}

# add the retention time required flashlfq for pepquery results
add_retention_time <- function(sps, spectrum, type, filename, mz) {
  
  mgf_parent_dir= "/mnt/HDD2/colin_work/microprotein_analysis/proteomics/pepquery/mgf_files/"
  mgf_file <- paste0(mgf_parent_dir, type, "/reducing/", filename)
  
  sps_sub <- filterDataOrigin(sps,dataOrigin=mgf_file)
  
  # convert to minutes for flashlfq
  retention_time = sps_sub[as.numeric(spectrum),]$rtime/60
  
  # compare the pepquery and mgf files to ensure the masses match
  if(round(mz,4) - round(sps_sub[as.numeric(spectrum),]$precursorMz,4) != 0) {
    print("masses don't match")
    
    print(c(as.numeric(mz), as.numeric(round(sps_sub[as.numeric(spectrum),]$precursorMz,4))))
    
  } else {
    print("masses match")
  }
  # print(paste0("mass_from mgf:",sps_sub[as.numeric(spectrum),]$precursorMz))
  
  return(retention_time)
}