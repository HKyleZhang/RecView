#' Index '012'-formatted genotype file
#' @description Index the genotype file for faster access in RecView ShinyApp.
#' @usage index_gt012(gt012_file, scaffold_file)
#' @param gt012_file the output of make_gt012() or make_012gt_from_vcf()
#' @param scaffold_file a file with the same format as used in RecView ShinyApp.
#' @note The file will be saved with the extension of .idx
#'
#' @export
index_gt012 <- function(gt012_file, scaffold_file) {
  sc <- read_csv(scaffold_file, col_types = cols(CHR = col_character(),
                                                         order = col_double(),
                                                         scaffold = col_character(),
                                                         size = col_double(),
                                                         orientation = col_character()))
  
  header <- read_csv(gt012_file, n_max = 1, col_types = cols()) %>% colnames()
  laf <- LaF::laf_open_csv(gt012_file, column_types = c(rep("character", 2), "numeric", rep("character", 2), rep("numeric", length(header) - 5)), skip = 1)
  num_rows <- LaF::nrow(laf)
  
  chr_list <- list()
  for (i in seq(1, num_rows, 1e4)) {
    chr_list[[length(chr_list)+1]] <- LaF::next_block(laf, nrows = 1e4) %>%
      mutate(row = i + seq(1,nrow(.)) - 1) %>%
      select(V2, row) %>%
      left_join(sc %>% select(scaffold, CHR), by = c("V2" = "scaffold")) %>%
      summarise(POS = str_flatten_comma(row), .by = "CHR")
  }
  out <- reduce(chr_list, bind_rows) %>% 
    summarise(POS = str_flatten_comma(POS), .by = "CHR")
  
  write_tsv(out, file = paste0(gt012_file, ".idx"))
}
