#' Locate recombination with two-generation pedigree
#' @description Locate recombination with two-generation pedigree that consists of father, mother and at least two offspring.
#' @usage rec_2gen(data, scaffold_info, chromosome = 1, offspring = c("off1", "off2", "off3"), method = "changepoint")
#' @param data a data frame or tibble of '012'-formatted genotype, which can be made using make_012gt() or make_012gt_from_vcf().
#' @param scaffold_info a data frame or tibble of scaffolds' order and orientation.
#' @param chromosome chromosome.
#' @param offspring the offspring to be included in the analysis.
#' @param method the algorithm to identify the change point. "CCS", "PD", or "changepoint".
#' @param method_args more arguments passed down to the algorithm that locates recombination.
#' @note This function produces a S3 class "RecView".
#' 
#' @export
rec_2gen <- function(data, scaffold_info, chromosome, offspring, method = "changepoint", method_args = NULL) {
  if(is.null(method_args)) {
    if (method == "CCS") method_args <- list(threshold = 50)
    if (method == "PD") method_args <- list(window_size = 1100, threshold = 0.8)
  } else {
    if (method == "PD") {
      if (is.null(method_args$window_size)) method_args <- list(window_size = 1100, threshold = as.vector(unlist(method_args)))
      if (is.null(method_args$threshold)) method_args <- list(threshold = 0.8, window_size = as.vector(unlist(method_args)))
      method_args <- list(window_size = method_args$window_size, threshold = method_args$threshold)
    }
  }
  
  # functions ----
  calc_distm <- function(x) {
    x <- str_split_1(x, "_")
    dm <- suppressWarnings(as.matrix(dist(x)))
    return(dm)
  }
  
  check_overlaps <- function(x, y) {
    x <- sort(x)
    y <- sort(y)
    x_tb <- tibble(seq = as.integer(seq(x[1], x[2])))
    y_tb <- tibble(seq = as.integer(seq(y[1], y[2])))
    olap <- nrow(inner_join(x_tb, y_tb, by = "seq"))
    if (olap > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  get_cp <- function(input, data_in, offspring_in, method, method_args) {
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
                                       start_POS = data_in_mod_tmp$POS_chr[changepoint::cpts(cpt)],
                                       end_POS = data_in_mod_tmp$POS_chr[changepoint::cpts(cpt)])
      }
    } else if (method == "CCS") {
      for (i in colnames(pairwise_tb)) {
        data_in_mod_tmp <- data_in_mod %>% filter(is.na(get(i,.)) == FALSE)
        cps[[length(cps)+1]] <- CCS_algorithm(x = data_in_mod_tmp, value_col = i, threshold = ifelse(is.null(method_args$threshold), 50, method_args$threshold), full_result = FALSE) %>% 
          mutate(cps_id = length(cps)+1, start_POS = data_in_mod_tmp$POS_chr[start_row], end_POS = data_in_mod_tmp$POS_chr[end_row]) %>% 
          select(cps_id, start_POS, end_POS)
      }
    } else if (method == "PD") {
      for (i in colnames(pairwise_tb)) {
        data_in_mod_tmp <- data_in_mod %>% filter(is.na(get(i,.)) == FALSE)
        cpt <- PD_algorithm(x = data_in_mod_tmp, value_col = i, window_size = ifelse(is.null(method_args$window_size), 1100, method_args$window_size), threshold = ifelse(is.null(method_args$threshold), 0.8, method_args$threshold), full_result = FALSE)
        cps[[length(cps)+1]] <- tibble(cps_id = length(cps)+1,
                                       start_POS = data_in_mod_tmp$POS_chr[cpt],
                                       end_POS = data_in_mod_tmp$POS_chr[cpt])
      }
    }
    
    
    position_result <- reduce(cps, bind_rows) %>% 
      mutate(mid_POS = (start_POS + end_POS) / 2)
    
    if (nrow(position_result) == 1) {
      position_result <- position_result %>% 
        mutate(Mean_bp = (start_POS + end_POS)/2, SE_bp = as.numeric(NA)) %>% 
        select(ID, Mean_bp, SE_bp)
    } else if (nrow(position_result) > 1) {
      range_seq <- list()
      for (i in 1:nrow(position_result)) {
        range_seq[[length(range_seq)+1]] <- seq(position_result$start_POS[i], position_result$end_POS[i])
      }
      
      overlaps <- 1
      for (i in 2:length(range_seq)) {
        if (length(which((range_seq[[1]] %in% range_seq[[i]]) == TRUE)) > 0) overlaps <- overlaps + 1
      }
      
      if (overlaps != nrow(position_result) && nrow(position_result) > 2) {
        kmeans_res <- kmeans(x = position_result$mid_POS, centers = 2, nstart = 2)
        ss_max <- round(kmeans_res$betweenss / kmeans_res$totss, 3)
        continue <- TRUE
        centers <- 2
        while (continue) {
          centers <- centers + 1
          if (centers <= length(unique(position_result$mid_POS)) && nrow(position_result) - centers >= 1) {
            kmeans_res <- kmeans(x = position_result$mid_POS, centers = centers, nstart = centers, iter.max = 25)
            ss_prop <- round(kmeans_res$betweenss / kmeans_res$totss, 3)
            if (ss_prop - ss_max > 0.001) {
              continue <- TRUE
              ss_max <- ss_prop
            } else {
              continue <- FALSE
              kmeans_res <- kmeans(x = position_result$mid_POS, centers = centers-1, nstart = centers-1)
              position_result <- position_result %>% 
                mutate(ID = kmeans_res$cluster)
            }
          } else {
            continue <- FALSE
            kmeans_res <- kmeans(x = position_result$mid_POS, centers = centers-1, nstart = centers-1)
            position_result <- position_result %>% 
              mutate(ID = kmeans_res$cluster)
          }
        }
      } else {
        position_result <- position_result %>% 
          mutate(ID = seq(1,nrow(.)))
      }

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
  chr <- scaffold_info %>%
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
    mutate(output = map(input, get_cp, data_chr_mod, offspring, method, method_args)) %>% 
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
    
    overlaps_matrix <- matrix(rep("", times = nrow(position_result_pat)^2), nrow = nrow(position_result_pat))
    for (i in 1:nrow(position_result_pat)) {
      for (j in i:nrow(position_result_pat)) {
        rangeA <- c(position_result_pat$Mean_bp[i] - position_result_pat$SE_bp[i], position_result_pat$Mean_bp[i] + position_result_pat$SE_bp[i])
        rangeB <- c(position_result_pat$Mean_bp[j] - position_result_pat$SE_bp[j], position_result_pat$Mean_bp[j] + position_result_pat$SE_bp[j])
        overlaps_matrix[j,i] <- check_overlaps(x = rangeA, y = rangeB)
      }
    }
    
    all_i <- vector()
    remove_i <- vector()
    for (i in 1:(nrow(position_result_pat)-1)) {
      j <- seq(i+1, nrow(position_result_pat))
      if (!( i %in% all_i)) {
        remove_i <- unique(c(remove_i, j[which(overlaps_matrix[j,i] == "TRUE")]))
      }
      all_i <- unique(c(all_i, which(overlaps_matrix[i:nrow(position_result_pat),i] == "TRUE")+i-1))
    }
    
    position_result_pat <- position_result_pat %>% 
      slice(which((seq(1,nrow(position_result_pat)) %in% remove_i) == FALSE))
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
    mutate(output = map(input, get_cp, data_chr_mod, offspring, method, method_args)) %>% 
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
    
    overlaps_matrix <- matrix(rep("", times = nrow(position_result_mat)^2), nrow = nrow(position_result_mat))
    for (i in 1:nrow(position_result_mat)) {
      for (j in i:nrow(position_result_mat)) {
        rangeA <- c(position_result_mat$Mean_bp[i] - position_result_mat$SE_bp[i], position_result_mat$Mean_bp[i] + position_result_mat$SE_bp[i])
        rangeB <- c(position_result_mat$Mean_bp[j] - position_result_mat$SE_bp[j], position_result_mat$Mean_bp[j] + position_result_mat$SE_bp[j])
        overlaps_matrix[j,i] <- check_overlaps(x = rangeA, y = rangeB)
      }
    }
    
    all_i <- vector()
    remove_i <- vector()
    for (i in 1:(nrow(position_result_mat)-1)) {
      j <- seq(i+1, nrow(position_result_mat))
      if (!( i %in% all_i)) {
        remove_i <- unique(c(remove_i, j[which(overlaps_matrix[j,i] == "TRUE")]))
      }
      all_i <- unique(c(all_i, which(overlaps_matrix[i:nrow(position_result_mat),i] == "TRUE")+i-1))
    }
    
    position_result_mat <- position_result_mat %>% 
      slice(which((seq(1,nrow(position_result_mat)) %in% remove_i) == FALSE))
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
              method = method,
              method_args = unlist(method_args))
  
  cat('\nDone!')
  
  attr(out, "class") <- "RecView"
  return(out)
}
