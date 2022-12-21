#' Make '012'-formatted genotype file
#' @description Make the genotype file in designated '012' format for RecView ShinyApp.
#' @usage make_012gt(gt, save_filename)
#' @param gt the output file when using vcftools to extract "GT"
#' @param F0F1_rename a vector of the labels of F0 (father's father, father's mother, mother's father, mother's mother), and F1 (father, mother)
#' @param save_filename the name for saving the '012'-formatted output, as .csv .
#' @note The file will be saved as .csv
#'
#' @export
make_012gt <- function(gt, F0F1_rename = "no_rename", save_filename = "gt_012.csv") {
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

  if ((F0F1_rename != "no_rename") && (length(F0F1_rename) == 6)) {
    out_new <- out %>%
      select(-F0F1_rename) %>%
      mutate(A = get(F0F1_rename[1], out),
             B = get(F0F1_rename[2], out),
             C = get(F0F1_rename[3], out),
             D = get(F0F1_rename[4], out),
             AB = get(F0F1_rename[5], out),
             CD = get(F0F1_rename[6], out)) %>%
      select(id, CHROM, POS, A, B, C, D, AB, CD, everything())

    write_csv(out_new, file = save_filename)
  } else {
    write_csv(out, file = save_filename)
  }
}
