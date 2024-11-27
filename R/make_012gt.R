#' Make '012'-formatted genotype file
#' @description Make the genotype file in designated '012' format for RecView ShinyApp.
#' @usage make_012gt(gt, F0F1_rename, save_filename)
#' @param gt the output file when using vcftools to extract "GT"
#' @param F0F1_rename a vector of the labels of F0 (father's father, father's mother, mother's father, mother's mother) and F1 (father, mother); use "NA" if missing.
#' @param save_filename the name for saving the '012'-formatted output, as .rds .
#' @param readable if TRUE, a CSV file will also be produced.
#' @note The file will be saved as .rds
#'
#' @export
make_012gt <- function(gt, F0F1_rename = "no_rename", save_filename = "gt_012.rds", readable = FALSE) {
  gt <- read_tsv(gt, col_types = cols())
  
  out <- tibble(CHROM = gt$CHROM, POS = gt$POS) %>%
    unite(col = "id", c("CHROM", "POS"), sep = ':', remove = FALSE) %>%
    select(id, CHROM, POS)
  
  gt <- select(gt, -CHROM, -POS)
  col_names <- colnames(gt)
  for (i in 1:length(col_names)) {
    left <- str_sub(get(col_names[i], gt), start = 1, end = 1) %>% as.numeric()
    right <- str_sub(get(col_names[i], gt), start = 3, end = 3) %>% as.numeric()
    
    out <- out %>% mutate(!!col_names[i] := (left + right))
  }
  
  if (F0F1_rename[1] != "no_rename") {
    miss_ind_index <- which(F0F1_rename == "NA")
    out_new <- out %>%
      select(-any_of(F0F1_rename)) %>%
      mutate(A = ifelse(rep("1" %in% miss_ind_index, times = nrow(.)), NA, get(F0F1_rename[1], out)),
             B = ifelse(rep("2" %in% miss_ind_index, times = nrow(.)), NA, get(F0F1_rename[2], out)),
             C = ifelse(rep("3" %in% miss_ind_index, times = nrow(.)), NA, get(F0F1_rename[3], out)),
             D = ifelse(rep("4" %in% miss_ind_index, times = nrow(.)), NA, get(F0F1_rename[4], out)),
             AB = ifelse(rep("5" %in% miss_ind_index, times = nrow(.)), NA, get(F0F1_rename[5], out)),
             CD = ifelse(rep("6" %in% miss_ind_index, times = nrow(.)), NA, get(F0F1_rename[6], out))) %>% 
      nest(everything_else = !any_of(c("A", "B", "C", "D", "AB", "CD"))) %>% 
      mutate(Missing_ind_num = 0, Missing_ind = "")
    
    for(i in 1:nrow(out_new)) {
      miss_ind_index_per_row <- which(is.na(c(out_new$A[i], out_new$B[i], out_new$C[i], out_new$D[i], out_new$AB[i], out_new$CD[i])))
      out_new$Missing_ind_num[i] <- length(miss_ind_index_per_row)
      out_new$Missing_ind[i] <- str_c(c("A", "B", "C", "D", "AB", "CD")[miss_ind_index_per_row], collapse = "_")
    }
    
    out_new <- out_new %>% 
      unnest(cols = "everything_else") %>% 
      select(id, CHROM, POS, Missing_ind_num, Missing_ind, A, B, C, D, AB, CD, everything())
  }
  
  if (F0F1_rename[1] != "no_rename") {
    if (readable) {
      write_csv(out_new, file = str_replace(save_filename, ".rds", ".csv"))
      fst::write_fst(x = out_new, path = save_filename, compress = 100)
    } else {
      fst::write_fst(x = out_new, path = save_filename, compress = 100)
    }
  } else {
    if (readable) {
      write_csv(out, file = str_replace(save_filename, ".rds", ".csv"))
      fst::write_fst(x = out_new, path = save_filename, compress = 100)
    } else {
      fst::write_fst(x = out_new, path = save_filename, compress = 100)
    }
  }
}
