
starting <- orftable_tis_filtered %>%
filter(!`ORF type` != "Annotated") 
dim(starting)[1]

# determine which ORFs are the only ORFs on a transcript
single_orf_per_transcript <- orftable_tis_filtered %>%
filter(!`ORF type` != "Annotated")  %>%
  group_by(`Transcript family`) %>%
  arrange(`Transcript family`) %>%
  mutate(count_orf = dplyr::n()) %>%
  dplyr::filter(count_orf == 1)
dim(single_orf_per_transcript)[1]

# determine instances where more than one novel ORFs are found on a transcript
multiple_orfs_per_transcript <- orftable_tis_filtered %>%
filter(!`ORF type` != "Annotated")  %>%
  group_by(`Transcript family`) %>%
  mutate(count_orf = dplyr::n())  %>%
  dplyr::filter(count_orf > 1)  %>%
  arrange(-`Length (AAs)`) # important as we start with the longest ORF when comparing overlaps below

dim(multiple_orfs_per_transcript)[1]

# for each set of overlapping ORFs divide to pairs to determine overlap
orf_granges <- GRanges(seqnames = multiple_orfs_per_transcript$`Transcript family`, 
        ranges = IRanges(start = multiple_orfs_per_transcript$`Transcript start position`, 
                         end = multiple_orfs_per_transcript$`Transcript stop position`))

names(orf_granges) <- multiple_orfs_per_transcript$`ORF-RATER name` 

non_overlapping_orfs <- c()
orf_pairs <- data.frame()

k=0 # iterator for no overlapping hits
while (length(orf_granges) > 0) {
  # print(length(orf_granges))
  if (length(orf_granges) > 1) {  
    
    # take reference off top of current iteration
    reference_orf <- orf_granges[1] 
    
    # remove ref from comparison set
    target_gr <- orf_granges[-1] 
    orf_granges <- target_gr

    if (length(findOverlaps(reference_orf, target_gr)) > 0) {
  
      # if the ref has overlaps, get the index
      hit_index <- subjectHits(findOverlaps(reference_orf, target_gr))
  
      # loop through hits and classifiy
      for (j in hit_index) {

        # overlap_type <- defineIsoform(reference_orf,target_gr[j]) 
        # print(names(target_gr[target_gr[hit_index[j]]))
        
        ref_orfrater <- orftable_tis_filtered %>%
          filter(`ORF-RATER name` == names(reference_orf)) 
        
         gr_orfrater <- orftable_tis_filtered %>%
          filter(`ORF-RATER name` == names(target_gr[j]))
        
        
        # push onto array
        overlap_def <-  data.frame(
         tid_ref=ref_orfrater$`Transcript family`,
          tid_target=gr_orfrater$`Transcript family`,
          ref=names(reference_orf), 
          target = names(target_gr[j]), 
          ref_start=ref_orfrater$`Transcript start position`,
          target_start=gr_orfrater$`Transcript start position`,
          ref_stop=ref_orfrater$`Transcript stop position`,
          target_stop=gr_orfrater$`Transcript stop position`,
          reference_codon=ref_orfrater$`Start codon`,
          target_codon=gr_orfrater$`Start codon`,
          reference_rhar = ref_orfrater$Rharr_min_Rnd,
          target_rhar =  gr_orfrater$Rharr_min_Rnd)
        
        orf_pairs <- bind_rows(orf_pairs, overlap_def)
      }
      # when finished remove the hits from test as we know the overlap with longest
      orf_granges <- orf_granges[-hit_index]
      
    } else {
        k=k+1
        non_overlapping_orfs[k] <- names(reference_orf[1])
        
         }
  } else{
    # if only left must by definition by non-overlapping
      k=k+1
      non_overlapping_orfs[k] <- names(orf_granges[1])
      orf_granges <- orf_granges[0]
  }
}

dim(single_orf_per_transcript)[1] + length(unique(c(orf_pairs$ref, orf_pairs$target))) + length(non_overlapping_orfs)

# some overlaps have the same start and end, keep one for each
# we remove these overlaps, but come each ref/target but will remove later if they  
# do not meet the criteria of the inside/outside TIS window for other ORFs they overlap with
perfect_match <- orf_pairs %>%
  dplyr::filter((tid_ref == tid_target) & (ref_start==target_start) & (ref_stop == target_stop)) %>%
  distinct(tid_ref,.keep_all=T) %>%
  dplyr::select(ref)
dim(orf_pairs)[1]

# remove perfect matches from the orf_pairs
orf_pairs <- orf_pairs %>%
  dplyr::filter(!((tid_ref == tid_target) & (ref_start==target_start) & (ref_stop == target_stop)))
dim(orf_pairs)[1]

# some overlaps have the same start but different ends
# we remove these overlaps from the potentials now, but remove the if they  
# do not meet the criteria of the inside/outside TIS window filters for other ORFs they overlap with
different_stops <- orf_pairs %>%
  filter((tid_ref == tid_target) & (ref_start==target_start) & (ref_stop != target_stop)) 

# remove different stops
orf_pairs <- orf_pairs %>%
  filter(!((tid_ref == tid_target) & (ref_start==target_start) & (ref_stop != target_stop)))
dim(orf_pairs)[1]

# for the remaining ORF pairs, determine if these ORFs occur within the TIS window (7 nt around the TIS)

with_start_codon  <- orf_pairs %>%
  arrange(ref) %>%
  mutate(tis_diff=abs(target_start-ref_start)) %>%
  mutate(tis_window = case_when(
    (tis_diff < -3) | (tis_diff > 3) ~ "outside",
    (tis_diff >= -3)  | (tis_diff <= 3)  ~ "inside")
  )

table(with_start_codon$tis_window)

inside_window <- with_start_codon %>%
  filter(tis_window == "inside") %>%
  rowwise() %>%
  mutate(ref_atg_tis = ifelse("ATG" %in% `reference_codon`, TRUE, FALSE)) %>%
  mutate(target_atg_tis = ifelse("ATG" %in% `target_codon`, TRUE, FALSE)) %>%
  mutate(max_value = pmax(reference_rhar, target_rhar),
  result = if_else(max_value == reference_rhar, ref, target)) %>%
  mutate(decision = case_when(
    (ref_atg_tis == T & target_atg_tis == F) ~ ref, 
    (target_atg_tis == T & ref_atg_tis == F) ~ target,
    (target_atg_tis == T & ref_atg_tis == T) ~ result,
    (target_atg_tis == F & ref_atg_tis == F) ~ result
  )) 

all_inside <- unique(c(inside_window$ref, inside_window$target))
length(all_inside)

inside_window_removed <- all_inside[!all_inside %in% inside_window$decision]
length(inside_window_removed)

outside_table <- with_start_codon %>%
  filter(tis_window == "outside") 

all_outside_orfs <- unique(c(outside_table$ref, outside_table$target))

outside_table <- with_start_codon %>%
  filter(tis_window == "outside") %>%
  rowwise() %>%
  mutate(ref_atg_tis = ifelse("ATG" %in% `reference_codon`, TRUE, FALSE)) %>%
  mutate(target_atg_tis = ifelse("ATG" %in% `target_codon`, TRUE, FALSE)) %>%
  mutate(decision = case_when(
    (ref_atg_tis == T & target_atg_tis == F) & (reference_rhar*5)  > target_rhar  ~ "remove_target",
    (target_atg_tis == T & ref_atg_tis == F) & (target_rhar*5) > reference_rhar  ~ "remove_reference",
    (target_atg_tis == T & ref_atg_tis == T)   ~ "keep",
    (target_atg_tis == F & ref_atg_tis == F)   ~ "keep"
  ))

table(outside_table$decision)

outside_table_removed <- outside_table %>%
filter(decision != "keep") %>%
mutate(removed_orf = ifelse(decision == "remove_target", target, ref))

outside_tis_window_removed <- unique(outside_table_removed$removed_orf)

overlapping_orf_filter <- c(outside_tis_window_removed, inside_window_removed)

a <- orftable_tis_filtered %>%
filter(!`ORF-RATER name` %in% unique(overlapping_orf_filter)) %>%
filter(`ORF-RATER name` %in% unique(drug_product_psms$protein))


unique(drug_product_psms$protein)[!unique(drug_product_psms$protein) %in% a$`ORF-RATER name` ]