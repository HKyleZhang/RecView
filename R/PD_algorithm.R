#' Proportional Difference (PD) algorithm
#' @description Proportional Difference (PD) algorithm for locating change points.
#' @usage PD_algorithm(x, value_col = "...", window_size = 400, min_threshold = 0.5, full_result = FALSE)
#' @param x a data frame or tibble.
#' @param value_col character; the column in which the change points will be analysed. This column should contain binary data, e.g. 0 and 1.
#' @param window_size numeric; the size of sliding windows in which the proportion of symbol is checked. Default: 400.
#' @param min_threshold numeric; minimal threshold to detect optimal threshold for locating the change points. Default: 0.5.
#' @param full_result logical; it indicates whether the full results will be produced.
#' 
#' @export
PD_algorithm <- function(x, value_col, window_size = 400, min_threshold = 0.5, full_result = FALSE) {
  seqle_mod <- function(x) {
    n <- length(x)  
    y <- x[-1L] != x[-n] + 1
    i <- c(which(y|is.na(y)),n) 
    return(x[head(c(0L,i)+1L,-1L)])
  }
  
  find_local_maxima <- function(tb) {
    tb_max <- tb[which(tb$diff == max(tb$diff)),]
    if (nrow(tb_max) > 1) {
      tb_max <- tb_max[nrow(tb_max) %/% 2,]
    }
    return(tb_max)
  }
  
  get_proportion <- function(input, tb, value_col, symbol) {
    symbols_vector <- pull(tb, value_col)[input$seg_start[1]:input$seg_end[1]]
    proportion <- length(symbols_vector[which(symbols_vector == symbol)]) / length(symbols_vector)
    return(proportion)
  }
  
  finer_running_difference <- function(input, proportion_diff_tb, x, value_col, symbol, radius, threshold) {
    i <- input$row_start[1]
    while (proportion_diff_tb$diff[i] >= threshold && i < nrow(proportion_diff_tb)) {
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
  
  jonas_elbow_finder <- function(x_values, y_values) {
    # Max values to create line
    max_x_x <- max(x_values)
    max_x_y <- y_values[which.max(x_values)]
    max_y_y <- max(y_values)
    max_y_x <- x_values[which.max(y_values)]
    max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
    
    # Creating straight line between the max values
    fit <- lm(max_df$y ~ max_df$x)
    
    # Distance from point to line
    distances <- c()
    for(i in 1:length(x_values)) {
      distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
    }
    
    # Max distance point
    x_max_dist <- unname(which.max(distances))
    y_max_dist <- unname(which.max(distances))
    
    return(c(x_max_dist, y_max_dist))
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
  
  num_peak_tb <- tibble(threshold = seq(0.95, min_threshold, -0.05), num_peak = 0)
  for (i in 1:nrow(num_peak_tb)) {
    rows_above_thrshd <- seqle_mod(which(segment_proportion_diff$diff >= num_peak_tb$threshold[i])) %>% 
      discard(is.na)
    if (length(rows_above_thrshd) > 0) {
      segment_proportion_finer <- tibble(row_start = rows_above_thrshd) %>% 
        mutate(id = seq(1, nrow(.))) %>% 
        nest(input = !id) %>% 
        mutate(output = map(input, finer_running_difference, segment_proportion_diff, x, value_col, symbol, radius, num_peak_tb$threshold[i])) %>% 
        mutate(local_maxima = map(output, find_local_maxima)) %>% 
        select(output, local_maxima)
    } else {
      segment_proportion_finer <- tibble(output = tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric()),
                                         local_maxima = tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric()))
    }
    
    num_peak_tb$num_peak[i] <- segment_proportion_finer %>% 
      unnest(cols = "local_maxima") %>% 
      pull(cut_point) %>% 
      length()
  }
  
  uniq_num_peak <- unique(num_peak_tb$num_peak)
  if (length(uniq_num_peak) == 1) {
    optimal_threshold <- num_peak_tb$threshold[1]
  } else if (length(uniq_num_peak) == 2) {
    optimal_threshold <- num_peak_tb$threshold[which.max(diff(num_peak_tb$num_peak))+1]
  } else {
    elbow_index <- jonas_elbow_finder(num_peak_tb$threshold, num_peak_tb$num_peak)[1]
    if (is.na(elbow_index)) elbow_index <- 0
    optimal_threshold <- num_peak_tb$threshold[elbow_index + 1]
  }
  
  rows_above_thrshd <- seqle_mod(which(segment_proportion_diff$diff >= optimal_threshold)) %>% 
    discard(is.na)
  if (length(rows_above_thrshd) > 0) {
    segment_proportion_finer <- tibble(row_start = rows_above_thrshd) %>% 
      mutate(id = seq(1, nrow(.))) %>% 
      nest(input = !id) %>% 
      mutate(output = map(input, finer_running_difference, segment_proportion_diff, x, value_col, symbol, radius, optimal_threshold)) %>% 
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
    out <- list(PD = segment_proportion_diff, PD_local = segment_proportion_finer, changepoint = cps)
  } else {
    out <- cps
  }
  
  return(out)
}
