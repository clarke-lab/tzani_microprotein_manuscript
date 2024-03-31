
starting <- orftable_tis_filtered %>%
filter(`ORF type` !=  "Annotated")

# determine which ORFs are the only ORFs on a transcript
sole_novel_orfs <- orftable_tis_filtered %>%
filter(`ORF type` !=  "Annotated") %>%
  group_by(`Transcript family`) %>%
  arrange(`Transcript family`) %>%
  mutate(count_orf = dplyr::n()) %>%
 dplyr::filter(count_orf == 1)

dim(sole_novel_orfs)[1]

# determine instances where more than one novel ORFs are found on a transcript
potential_overlapping_orfs <- orftable_tis_filtered %>%
filter(`ORF type` !=  "Annotated") %>% # 
  group_by(`Transcript family`) %>%
  mutate(count_orf = dplyr::n())  %>%
  dplyr::filter(count_orf > 1)  %>%
  arrange(-`Length (AAs)`) # important as we start with the longest ORF when comparing overlaps below

dim(potential_overlapping_orfs)[1]

# for each set of overlapping ORFs divide to pairs to determine overlap
overlap_orf_granges <- GRanges(seqnames = potential_overlapping_orfs$`Transcript family`, 
        ranges = IRanges(start = potential_overlapping_orfs$`Transcript start position`, 
                         end = potential_overlapping_orfs$`Transcript stop position`))

names(overlap_orf_granges) <- potential_overlapping_orfs$`ORF-RATER name` 


non_overlapping_orfs <- c()
overlap_pairs <- data.frame()

k=0 # iterator for no overlapping hits
while (length(overlap_orf_granges) > 0) {
  # print(length(overlap_orf_granges))
  if (length(overlap_orf_granges) > 1) {  
    
    # take reference off top of current iteration
    reference_orf <- overlap_orf_granges[1] 
    
    # remove ref from comparison set
    target_gr <- overlap_orf_granges[-1] 
    overlap_orf_granges <- target_gr

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
        
        overlap_pairs <- bind_rows(overlap_pairs, overlap_def)
      }
      # when finished remove the hits from test as we know the overlap with longest
      overlap_orf_granges <- overlap_orf_granges[-hit_index]
      
    } else {
        k=k+1
        non_overlapping_orfs[k] <- names(reference_orf[1])
        
         }
  } else{
    # if only left must by definition by non-overlapping
      k=k+1
      non_overlapping_orfs[k] <- names(overlap_orf_granges[1])
      overlap_orf_granges <- overlap_orf_granges[0]
  }
}

dim(sole_novel_orfs)[1] + length(unique(c(overlap_pairs$ref, overlap_pairs$target))) + length(non_overlapping_orfs)

# some overlaps have the same start and end, keep one for each
perfect_match_retained <- overlap_pairs %>%
  dplyr::filter((tid_ref == tid_target) & (ref_start==target_start) & (ref_stop == target_stop)) %>%
  distinct(tid_ref,.keep_all=T) %>%
  dplyr::select(ref)

# some overlaps have the same start but different ends
different_stops <- overlap_pairs %>%
  filter((tid_ref == tid_target) & (ref_start==target_start) & (ref_stop != target_stop)) 

different_stops <- unique(c(different_stops$ref, different_stops$target))

# start codon window
with_start_codon  <- overlap_pairs %>%
  arrange(ref) %>%
  mutate(tis_diff=abs(target_start-ref_start)) %>%
  mutate(tis_window = case_when(
    (tis_diff < -3) | (tis_diff > 3) ~ "outside",
    (tis_diff >= -3)  | (tis_diff <= 3)  ~ "inside")
  )

table(with_start_codon$tis_window)

inside_window <- with_start_codon %>%
  filter(tis_window == "inside") %>%
  mutate(max_value = pmax(reference_rhar, target_rhar),
  result = if_else(max_value == reference_rhar, ref, target)) 

all_inside <- unique(c(inside_window$ref, inside_window$target))
length(all_inside)

length(unique(inside_window$result))

inside_window_removed <- all_inside[!all_inside %in% inside_window$result]

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

remove_target <- outside_table %>% 
  filter(decision=="remove_target")

remove_reference <- outside_table %>% 
  filter(decision=="remove_reference")
