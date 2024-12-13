#' Cumulative Continuity Score (CCS) algorithm
#' @description Cumulative Continuity Score (CCS) algorithm for locating change points.
#' @usage CCS_algorithm(x, value_col = "...", threshold = 50, full_result = FALSE)
#' @param x a data frame or tibble.
#' @param value_col character; the column in which the change points will be analysed. This column should contain binary data, e.g. 0 and 1.
#' @param threshold numeric; the CCS value above which the segments will be included for locating the change points. Default = 50.
#' @param symbol character; the character for assigning negative scores.
#' @param full_result logical; it indicates whether the result having CCS at all rows will be produced, or only the change points otherwise.
#' 
#' @export
CCS_algorithm <- function(x, value_col, threshold = 50, symbol = NULL, full_result = FALSE) {
  if (nrow(x) != 0) {
    data <- x %>%
      mutate("{paste0(value_col, '_+_1')}" := c(0, pull(.,value_col))[1:length(pull(.,value_col))]) %>% 
      mutate(similarity = ifelse(pull(., value_col) == pull(., paste0(value_col, "_+_1")), "Y", "N")) %>% 
      mutate(index = seq(1,nrow(.)))
    
    N_rows <- which(pull(data, similarity) == "N")
    tb_Ns <- data[N_rows,] %>% 
      mutate(CCS = 0, abs_CCS = 0,
             start  = index, end = index,
             gt_threshold = ifelse(abs_CCS >= threshold, "Y", "N"))
    
    tb_Ys_list <- list()
    if (length(N_rows) == 1) {
      tb_Ys_tmp <- data[(N_rows[1] + 1):nrow(data),]
      
      if (unique(pull(tb_Ys_tmp, value_col)) != symbol) {
        tb_Ys_list[[length(tb_Ys_list)+1]] <- tb_Ys_tmp %>% 
          mutate(CCS = seq(1,nrow(.)), abs_CCS = rep(nrow(.), times = nrow(.)),
                 start = min(index) - 1, end = max(index),
                 gt_threshold = ifelse(abs_CCS >= threshold, "Y", "N"))
      } else if (unique(pull(tb_Ys_tmp, value_col)) == symbol) {
        tb_Ys_list[[length(tb_Ys_list)+1]] <- tb_Ys_tmp %>% 
          mutate(CCS = seq(-1,-nrow(.),-1), abs_CCS = rep(nrow(.), times = nrow(.)),
                 start = min(index) - 1, end = max(index),
                 gt_threshold = ifelse(abs_CCS >= threshold, "Y", "N"))
      }
    } else {
      
      for (j in 2:length(N_rows)) {
        if (N_rows[j] - 1 != N_rows[j-1])  {
          tb_Ys_tmp <- data[(N_rows[j-1] + 1):(N_rows[j] - 1),]
          
          if (unique(pull(tb_Ys_tmp, value_col)) != symbol) {
            tb_Ys_list[[length(tb_Ys_list)+1]] <- tb_Ys_tmp %>% 
              mutate(CCS = seq(1,nrow(.)), abs_CCS = rep(nrow(.), times = nrow(.)),
                     start = min(index) - 1, end = max(index),
                     gt_threshold = ifelse(abs_CCS >= threshold, "Y", "N"))
          } else if (unique(pull(tb_Ys_tmp, value_col)) == symbol) {
            tb_Ys_list[[length(tb_Ys_list)+1]] <- tb_Ys_tmp %>% 
              mutate(CCS = seq(-1,-nrow(.),-1), abs_CCS = rep(nrow(.), times = nrow(.)),
                     start = min(index) - 1, end = max(index),
                     gt_threshold = ifelse(abs_CCS >= threshold, "Y", "N"))
          }
        }
      }
    }
    
    if (N_rows[length(N_rows)] != nrow(data))  {
      tb_Ys_tmp <- data[(N_rows[length(N_rows)] + 1):nrow(data),]
      if (unique(pull(tb_Ys_tmp, value_col)) != symbol) {
        tb_Ys_list[[length(tb_Ys_list)+1]] <- tb_Ys_tmp %>% 
          mutate(CCS = seq(1,nrow(.)), abs_CCS = rep(nrow(.), times = nrow(.)),
                 start = min(index) - 1, end = max(index),
                 gt_threshold = ifelse(abs_CCS >= threshold, "Y", "N"))
      } else if (unique(pull(tb_Ys_tmp, value_col)) == symbol) {
        tb_Ys_list[[length(tb_Ys_list)+1]] <- tb_Ys_tmp %>% 
          mutate(CCS = seq(-1,-nrow(.),-1), abs_CCS = rep(nrow(.), times = nrow(.)),
                 start = min(index) - 1, end = max(index),
                 gt_threshold = ifelse(abs_CCS >= threshold, "Y", "N"))
      }
    }
    
    tb_gt_threshold <- list()
    
    for (i in 1:length(tb_Ys_list)) {
      if (unique(tb_Ys_list[[i]]$gt_threshold) == "Y") {
        tb_gt_threshold[[length(tb_gt_threshold)+1]] <- tb_Ys_list[[i]]
      }
    }
    
    if(length(tb_gt_threshold) > 1) {
      cps <- list()
      for (i in 1:(length(tb_gt_threshold)-1)) {
        if (tb_gt_threshold[[i]]$CCS[1]*tb_gt_threshold[[i+1]]$CCS[1] < 0) {
          cps[[length(cps)+1]] <- tibble(index = i, start_row = unique(tb_gt_threshold[[i]]$end), end_row = unique(tb_gt_threshold[[i+1]]$start))
        }
      }
      
      if (length(cps) > 0) {
        cps <- reduce(cps, bind_rows) %>% 
          mutate(id = seq(1, nrow(.))) %>% 
          select(-index) %>% 
          select(id, start_row, end_row)
      } else {
        cps <- tibble(id = as.numeric(), start_row = as.numeric(), end_row = as.numeric())
      }
    } else {
      cps <- tibble(id = as.numeric(), start_row = as.numeric(), end_row = as.numeric())
    }
    
    if (full_result) {
      res <- reduce(tb_Ys_list, bind_rows) %>%
        bind_rows(tb_Ns) %>%
        arrange(index) %>% 
        select(-index, -similarity, -abs_CCS, -start, -end, -gt_threshold, -paste0({{value_col}}, '_+_1'))
      
      out <- list(full_result = res, changepoint = cps)
    } else {
      out <- cps
    }
  } else {
    cat('Message: the input has 0 rows.')
  }
  
  return(out)
}
