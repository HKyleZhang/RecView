#' Locate recombination with two-generation pedigree
#' @description Locate recombination with two-generation pedigree that consists of father, mother and at least two offspring.
#' @usage rec_2gen(data, sc_order, chromosome = 1, offspring = c("off1", "off2", "off3"), grouping_distance = 1e4, method = "changepoint")
#' @param data a data frame or tibble of '012'-formatted genotype, which can be made using make_012gt() or make_012gt_from_vcf().
#' @param sc_order a data frame or tibble of scaffolds' order and orientation.
#' @param chromosome chromosome.
#' @param offspring the offspring to be included in the analysis.
#' @param grouping_distance the maximal base pairs in which the detected recombination positions will be grouped.
#' @param method the algorithm to identify the change point. "CCS", "PD", or "changepoint".
#' @param method.args more arguments passed down to the algorithm that locates recombination.
#' @note This function produces a S3 class "RecView".
#' 
#' @export
rec_2gen <- function(data, sc_order, chromosome, offspring, grouping_distance = 1e4, method = "changepoint", method.args = NULL) {
  if(is.null(method.args)) {
    if (method == "CCS") method.args <- list(threshold = 50)
    if (method == "PD") method.args <- list(window_size = 1100, threshold = 0.8)
  } else {
    if (method == "PD") {
      if (is.null(method.args$window_size)) method.args <- list(window_size = 1100, threshold = as.vector(unlist(method.args)))
      if (is.null(method.args$threshold)) method.args <- list(threshold = 0.8, window_size = as.vector(unlist(method.args)))
      method.args <- list(window_size = method.args$window_size, threshold = method.args$threshold)
    }
  }
  
  # functions ----
  calc_distm <- function(x) {
    x <- str_split_1(x, "_")
    dm <- suppressWarnings(as.matrix(dist(x)))
    return(dm)
  }
  
  conditional_overlaps <- function(x, y, threshold = 0.5) {
    x <- sort(x)
    y <- sort(y)
    x_seq <- seq(x[1], x[2])
    y_seq <- seq(y[1], y[2])
    olap <- length(which(x_seq %in% y_seq))
    if (olap / length(x_seq) >= threshold || olap / length(y_seq) >= threshold) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  get_cp <- function(input, data_in, offspring_in, grouping_distance, method, method.args) {
    dms <- data_in$dm
    index <- input$index[1]
    off_index <- which(seq(1,length(offspring_in)) != index)
    
    pairwise_tb <- tibble(group = as.character())
    for (i in off_index) {
      if (i > index) {
        pairwise_tb <- add_row(pairwise_tb, group = paste(offspring_in[index], offspring_in[i], sep = "_<>_"))
      } else {
        pairwise_tb <- add_row(pairwise_tb, group = paste(offspring_in[i], offspring_in[index], sep = "_<>_"))
      }
    }
    
    for(i in 1:length(dms)) {
      column <- dms[[i]][,index] %>% 
        .[seq(1,length(.)) != index] %>% 
        as_tibble() %>% 
        `colnames<-`(., paste0("m", i))
      pairwise_tb <- bind_cols(pairwise_tb, column)
    }
    
    pairwise_tb <- pairwise_tb %>% 
      column_to_rownames(var = "group") %>% 
      t() %>% 
      as_tibble()
    
    data_in_mod <- bind_cols(select(data_in, input), pairwise_tb) %>% 
      unnest(cols = "input") %>% 
      arrange(POS_chr)
    
    cps <- list()
    if (method == "changepoint") {
      for (i in colnames(pairwise_tb)) {
        data_in_mod_tmp <- data_in_mod %>% filter(is.na(get(i,.)) == FALSE)
        cpt <- changepoint::cpt.mean(pull(data_in_mod_tmp, i), minseglen = 1)
        cps[[length(cps)+1]] <- tibble(cps_id = length(cps)+1,
                                       start_POS = data_in_mod_tmp$POS_chr[changepoint::cpts(cpt)] - grouping_distance,
                                       end_POS = data_in_mod_tmp$POS_chr[changepoint::cpts(cpt)] + grouping_distance)
      }
    } else if (method == "CCS") {
      for (i in colnames(pairwise_tb)) {
        data_in_mod_tmp <- data_in_mod %>% filter(is.na(get(i,.)) == FALSE)
        cps[[length(cps)+1]] <- CCS_algorithm(x = data_in_mod_tmp, value_col = i, threshold = ifelse(is.null(method.args$threshold), 50, method.args$threshold), full_result = FALSE) %>% 
          mutate(cps_id = length(cps)+1, start_POS = data_in_mod_tmp$POS_chr[start_row], end_POS = data_in_mod_tmp$POS_chr[end_row]) %>% 
          select(cps_id, start_POS, end_POS)
      }
    } else if (method == "PD") {
      for (i in colnames(pairwise_tb)) {
        data_in_mod_tmp <- data_in_mod %>% filter(is.na(get(i,.)) == FALSE)
        cpt <- PD_algorithm(x = data_in_mod_tmp, value_col = i, window_size = ifelse(is.null(method.args$window_size), 1100, method.args$window_size), threshold = ifelse(is.null(method.args$threshold), 0.8, method.args$threshold), full_result = FALSE)
        cps[[length(cps)+1]] <- tibble(cps_id = length(cps)+1,
                                       start_POS = data_in_mod_tmp$POS_chr[cpt] - grouping_distance,
                                       end_POS = data_in_mod_tmp$POS_chr[cpt] + grouping_distance)
      }
    }
    
    
    position_result <- reduce(cps, bind_rows) %>% 
      mutate(ID = seq(1, nrow(.)), K = 0)
    
    if (nrow(position_result) == 1) {
      position_result <- position_result %>% 
        mutate(Mean_bp = (start_POS + end_POS)/2, SE_bp = as.numeric(NA)) %>% 
        select(ID, Mean_bp, SE_bp)
    } else if (nrow(position_result) > 1) {
      overlaps_matrix <- matrix(rep("", times = nrow(position_result)^2), nrow = nrow(position_result))
      for (i in 1:nrow(position_result)) {
        for (j in i:nrow(position_result)) {
          rangeA <- c(position_result$start_POS[i], position_result$end_POS[i])
          rangeB <- c(position_result$start_POS[j], position_result$end_POS[j])
          A_threshold <- ifelse(grouping_distance/abs(rangeA[2]-rangeA[1]) <= 1, grouping_distance/abs(rangeA[2]-rangeA[1]), abs(rangeA[2]-rangeA[1])/grouping_distance)
          B_threshold <- ifelse(grouping_distance/abs(rangeB[2]-rangeB[1]) <= 1, grouping_distance/abs(rangeB[2]-rangeB[1]), abs(rangeB[2]-rangeB[1])/grouping_distance)
          overlaps_matrix[j,i] <- conditional_overlaps(x = rangeA, y = rangeB, threshold = min(A_threshold, B_threshold))
        }
      }
      
      all_i <- vector()
      for (i in 1:nrow(position_result)) {
        if (!( i %in% all_i)) {
          position_result$K[which(overlaps_matrix[i:nrow(position_result),i] == "TRUE")+i-1] <- max(position_result$K)+1
        }
        all_i <- unique(c(all_i, which(overlaps_matrix[i:nrow(position_result),i] == "TRUE")+i-1))
      }
      position_result <- position_result %>% 
        select(cps_id, K, start_POS, end_POS) %>% 
        `colnames<-`(.,c("cps_id", "ID", "start_POS", "end_POS"))
      
      if (length(cps) == 1) {
        position_result <- position_result %>% 
          mutate(POS_chr = (start_POS + end_POS)/2) %>% 
          summarise(Mean_bp = mean(POS_chr), SE_bp = sd(POS_chr)/sqrt(length(POS_chr)) , .by = "ID") %>% 
          arrange(ID)
      } else if (length(cps) > 1) {
        position_result <- position_result %>% 
          summarise(start_POS = min(start_POS), end_POS = max(end_POS), .by = c("cps_id", "ID")) %>% 
          mutate(POS_chr = (start_POS + end_POS)/2) %>% 
          summarise(Mean_bp = mean(POS_chr), SE_bp = sd(POS_chr) / sqrt(length(POS_chr)) , n = length(POS_chr), .by = "ID") %>% 
          filter(n > 1) %>% 
          select(-n)
        
        if (nrow(position_result) > 0) {
          position_result <- position_result %>% mutate(ID = seq(1, nrow(.)))
        } else {
          position_result <- tibble(ID = as.numeric(), Mean_bp = as.numeric(), SE_bp = as.numeric())
        }
      }
    } else {
      position_result <- tibble(ID = as.numeric(), Mean_bp = as.numeric(), SE_bp = as.numeric())
    }
    
    comparison_result <- data_in_mod %>% 
      gather(key = "group", value = "Concord_Discord", colnames(pairwise_tb))
    
    out <- tibble(comparison_result = list(comparison_result), position_result = list(position_result))
    return(out)
  }
  
  # Chromosomal position ----
  start_time <- Sys.time()
  chr <- sc_order %>%
    filter(CHR == !!chromosome) %>%
    arrange(order) %>%
    mutate(accumulated_size = 0)
  
  if (nrow(chr) == 1) {
    chr$accumulated_size[1] <- 0
  } else if (nrow(chr) > 1) {
    for (i in 2:nrow(chr)) {
      chr$accumulated_size[i] <- sum(chr$size[1:(i-1)])
    }
  }
  
  data_chr <- semi_join(data, chr, by = c("CHROM" = "scaffold")) %>%
    left_join(chr %>% select(scaffold, size, orientation, accumulated_size), by = c("CHROM" = "scaffold")) %>%
    mutate(scaffold_orientation = paste0(CHROM, ' ', orientation))
  data_chr <- data_chr %>% arrange(factor(CHROM, levels = chr$scaffold))
  
  data_chr_minus <- data_chr %>%
    filter(orientation == "-") %>%
    mutate(POS_chr = accumulated_size + size - POS) %>%
    select(POS_chr, everything())
  
  data_chr_plus <- data_chr %>%
    filter(orientation == "+") %>%
    mutate(POS_chr = accumulated_size + POS) %>%
    select(POS_chr, everything())
  
  data_chr <- bind_rows(data_chr_minus, data_chr_plus) %>% 
    arrange(POS_chr)
  
  end_time <- Sys.time()
  msg <- paste0('Step 1/4: Chromosomal position conversion completed. Used ', round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "secs")), 3), " sec.")
  cat(msg)
  
  # Paternal ----
  start_time <- Sys.time()
  data_chr_mod <- data_chr %>% 
    select(id, CHROM, orientation, POS, POS_chr, AB, CD, all_of(offspring)) %>% 
    mutate(diffAB_CD = abs(AB - CD), sumAB_CD = AB + CD) %>% 
    filter(diffAB_CD == 1, AB == 1) %>% 
    rowwise() %>% 
    mutate(minAB_CD = min(AB, CD)) %>% 
    ungroup()
  
  for (i in offspring) {
    data_chr_tmp <- data_chr_mod %>% 
      mutate(condition_1 = get(i,.) <= sumAB_CD & get(i,.) >= minAB_CD) 
    
    offspring_gts <- get(i, data_chr_tmp)
    offspring_gts[which(get("condition_1", data_chr_tmp) == FALSE)] <- NA
    
    data_chr_mod <- data_chr_tmp %>% 
      mutate({{i}} := offspring_gts) %>% 
      mutate(condition_2 = is.na(get(i,.))) %>% 
      mutate(keep = condition_1 | condition_2) %>% 
      filter(keep == TRUE) %>% 
      select(-condition_1, -condition_2, -keep)
  }
  
  data_chr_mod <- data_chr_mod %>% 
    select(-diffAB_CD, -sumAB_CD, -minAB_CD) %>% 
    select(id, CHROM, orientation, POS, POS_chr, all_of(offspring)) %>% 
    unite(col = "GTS", all_of(offspring), sep = "_", remove = FALSE) %>% 
    nest(input = !GTS) %>% 
    rowwise() %>% 
    mutate(dm = list(calc_distm(GTS)))
  
  res_pat <- tibble(Offspring = offspring) %>% 
    mutate(index = seq(1, nrow(.))) %>% 
    nest(input = !Offspring) %>% 
    mutate(output = map(input, get_cp, data_chr_mod, offspring, grouping_distance, method, method.args)) %>% 
    unnest(cols = "output") %>% 
    select(-input) %>% 
    mutate(Note = '-', Side = "Paternal")
  
  if (nrow(res_pat) == 2) {
    position_result_pat <- res_pat %>% 
      slice(1) %>% 
      select(-comparison_result, -Offspring) %>% 
      unnest(cols = "position_result") %>% 
      select(Side, everything())
  } else {
    position_result_pat <- res_pat %>% 
      select(-comparison_result) %>% 
      unnest(cols = "position_result") %>% 
      select(Side, Offspring, everything())
    
    dup <- position_result_pat %>% 
      group_by(Mean_bp) %>% 
      tally() %>% 
      filter(n > 1) %>% 
      pull(Mean_bp)
    
    position_result_pat$Note[position_result_pat$Mean_bp %in% dup] <- "Likely_Problematic"
  }
  
  comparison_result_pat <- res_pat %>% 
    select(comparison_result, Side) %>% 
    unnest(cols = "comparison_result")
  
  end_time <- Sys.time()
  msg <- paste0('\nStep 2/4: Paternal chromosome analysis completed. Used ', round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "secs")), 3), " sec.")
  cat(msg)
  
  # Maternal ----
  start_time <- Sys.time()
  data_chr_mod <- data_chr %>% 
    select(id, CHROM, orientation, POS, POS_chr, AB, CD, all_of(offspring)) %>% 
    mutate(diffAB_CD = abs(AB - CD), sumAB_CD = AB + CD) %>% 
    filter(diffAB_CD == 1, CD == 1) %>% 
    rowwise() %>% 
    mutate(minAB_CD = min(AB, CD)) %>% 
    ungroup()
  
  for (i in offspring) {
    data_chr_tmp <- data_chr_mod %>% 
      mutate(condition_1 = get(i,.) <= sumAB_CD & get(i,.) >= minAB_CD) 
    
    offspring_gts <- get(i, data_chr_tmp)
    offspring_gts[which(get("condition_1", data_chr_tmp) == FALSE)] <- NA
    
    data_chr_mod <- data_chr_tmp %>% 
      mutate({{i}} := offspring_gts) %>% 
      mutate(condition_2 = is.na(get(i,.))) %>% 
      mutate(keep = condition_1 | condition_2) %>% 
      filter(keep == TRUE) %>% 
      select(-condition_1, -condition_2, -keep)
  }
  
  data_chr_mod <- data_chr_mod %>% 
    select(-diffAB_CD, -sumAB_CD, -minAB_CD) %>% 
    select(id, CHROM, orientation, POS, POS_chr, all_of(offspring)) %>% 
    unite(col = "GTS", all_of(offspring), sep = "_", remove = FALSE) %>% 
    nest(input = !GTS) %>% 
    rowwise() %>% 
    mutate(dm = list(calc_distm(GTS)))
  
  res_mat <- tibble(Offspring = offspring) %>% 
    mutate(index = seq(1, nrow(.))) %>% 
    nest(input = !Offspring) %>% 
    mutate(output = map(input, get_cp, data_chr_mod, offspring, grouping_distance, method, method.args)) %>% 
    unnest(cols = "output") %>% 
    select(-input) %>% 
    mutate(Note = '-', Side = "Maternal")
  
  if (nrow(res_mat) == 2) {
    position_result_mat <- res_mat %>% 
      slice(1) %>% 
      select(-comparison_result, -Offspring) %>% 
      unnest(cols = "position_result") %>% 
      select(Side, everything())
  } else {
    position_result_mat <- res_mat %>% 
      select(-comparison_result) %>% 
      unnest(cols = "position_result") %>% 
      select(Side, Offspring, everything())
    
    dup <- position_result_mat %>% 
      group_by(Mean_bp) %>% 
      tally() %>% 
      filter(n > 1) %>% 
      pull(Mean_bp)
    
    position_result_mat$Note[position_result_mat$Mean_bp %in% dup] <- "Likely_Problematic"
  }
  
  comparison_result_mat <- res_mat %>% 
    select(comparison_result, Side) %>% 
    unnest(cols = "comparison_result")
  
  end_time <- Sys.time()
  msg <- paste0('\nStep 3/4: Maternal chromosome analysis completed. Used ', round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "secs")), 3), " sec.")
  cat(msg)
  
  # Results collecting ----
  start_time <- Sys.time()
  position_result <- bind_rows(position_result_pat, position_result_mat)
  comparison_result <- bind_rows(comparison_result_pat, comparison_result_mat) 
  
  end_time <- Sys.time()
  msg <- paste0('\nStep 4/4: Results collected. Used ', round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "secs")), 3), " sec.")
  cat(msg)
  
  out <- list(chromosome = chromosome,
              offspring = offspring, 
              pairwise_comparison = comparison_result, 
              recombination_position = position_result,
              grouping_distance = grouping_distance,
              method = method,
              method.args = unlist(method.args))
  
  cat('\nDone!')
  
  attr(out, "class") <- "RecView"
  return(out)
}
