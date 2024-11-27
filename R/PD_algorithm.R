#' Proportional Difference (PD) algorithm
#' @description Proportional Difference (PD) algorithm for locating change points.
#' @usage PD_algorithm(x, value_col = "...", window_size = 1100, threshold = 0.8, full_result = FALSE)
#' @param x a data frame or tibble.
#' @param value_col character; the column in which the change points will be analysed. This column should contain binary data, e.g. 0 and 1.
#' @param window_size numeric; the size of sliding windows in which the proportion of symbol is checked. Default: 1100.
#' @param threshold numeric; the PD value above which a finer step will taken to locate the change points. Default = 0.8.
#' @param full_result logical; it indicates whether the full results will be produced.
#' 
#' @export
PD_algorithm <- function(x, value_col, window_size = 1100, threshold = 0.8, full_result = FALSE) {
  seqle_mod <- function(x) {
    n <- length(x)  
    y <- x[-1L] != x[-n] + 1
    i <- c(which(y|is.na(y)),n) 
    return(x[head(c(0L,i)+1L,-1L)])
  }
  
  find_local_maxima <- function(tb) {
    tb[which(tb$diff == max(tb$diff)),]
  }
  
  get_proportion <- function(input, tb, value_col, symbol) {
    symbols_vector <- pull(tb, value_col)[input$seg_start[1]:input$seg_end[1]]
    proportion <- length(symbols_vector[which(symbols_vector == symbol)]) / length(symbols_vector)
    return(proportion)
  }
  
  finer_running_difference <- function(input, proportion_diff_tb, x, value_col, symbol, radius, threshold) {
    i <- input$row_start[1]
    while (proportion_diff_tb$diff[i] >= threshold) {
      i <- i + 1
    }
    cut_point <- seq(proportion_diff_tb$cut_point[ifelse((input$row_start[1] - 1) == 0, 1, input$row_start[1] - 1)] + 1, proportion_diff_tb$cut_point[i] - 1, 1)
    seg_1_proportion <- tibble(cut_point = cut_point, seg_start = cut_point - radius, seg_end = cut_point - 1) %>% 
      nest(input = !cut_point) %>% 
      mutate(p_1 = map(input, get_proportion, x, value_col, symbol)) %>% 
      unnest(cols = c("input", "p_1"))
    seg_2_proportion <- tibble(cut_point = cut_point, seg_start = cut_point, seg_end = cut_point + radius - 1) %>% 
      nest(input = !cut_point) %>% 
      mutate(p_2 = map(input, get_proportion, x, value_col, symbol)) %>% 
      unnest(cols = c("input", "p_2"))
    out <- inner_join(seg_1_proportion, seg_2_proportion, by = "cut_point") %>% 
      mutate(diff = abs(p_1 - p_2))
    return(out)
  }
  
  radius <- window_size %/% 2
  rows <- nrow(x)
  if ((rows - radius * 2) > 1) {
    symbol <- pull(x, value_col) %>% 
      unique() %>% 
      .[1]
    step <- radius %/% 3
    cut_point <- seq(1, rows - radius * 2, step) + radius
    
    seg_1_proportion <- tibble(cut_point = cut_point, seg_start = cut_point - radius, seg_end = cut_point - 1) %>% 
      nest(input = !cut_point) %>% 
      mutate(p_1 = map(input, get_proportion, x, value_col, symbol)) %>% 
      unnest(cols = c("input", "p_1"))
    seg_2_proportion <- tibble(cut_point = cut_point, seg_start = cut_point, seg_end = cut_point + radius - 1) %>% 
      nest(input = !cut_point) %>% 
      mutate(p_2 = map(input, get_proportion, x, value_col, symbol)) %>% 
      unnest(cols = c("input", "p_2"))
    segment_proportion_diff <- inner_join(seg_1_proportion, seg_2_proportion, by = "cut_point") %>% 
      mutate(diff = abs(p_1 - p_2))
  } else {
    segment_proportion_diff <- tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric())
  }
  
  rows_above_thrshd <- seqle_mod(which(segment_proportion_diff$diff >= threshold)) %>% 
    discard(is.na)
  if (length(rows_above_thrshd) > 0) {
    segment_proportion_finer <- tibble(row_start = rows_above_thrshd) %>% 
      mutate(id = seq(1, nrow(.))) %>% 
      nest(input = !id) %>% 
      mutate(output = map(input, finer_running_difference, segment_proportion_diff, x, value_col, symbol, radius, threshold)) %>% 
      mutate(local_maxima = map(output, find_local_maxima)) %>% 
      select(output, local_maxima)
  } else {
    segment_proportion_finer <- tibble(output = tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric()),
                                           local_maxima = tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric()))
  }
  
  cps <- segment_proportion_finer %>% 
    unnest(cols = "local_maxima") %>% 
    pull(cut_point)
  
  if (full_result) {
    out <- list(PD = segment_proportion_diff, PD_local = segment_goo_proportion_finer, changepoint = cps)
  } else {
    out <- cps
  }
  
  return(out)
}