#' Locate recombination with two-generation pedigree
#' @description Locate recombination with two-generation pedigree that consists of father, mother and at least two offspring.
#' @usage rec2gen(data, scaffold_info, chromosome = 1, offspring = c("off1", "off2", "off3"), method = "changepoint")
#' @param data a data frame or tibble of '012'-formatted genotype, which can be made using make_012gt() or make_012gt_from_vcf().
#' @param scaffold_info a data frame or tibble of scaffolds' order and orientation.
#' @param chromosome chromosome.
#' @param offspring the offspring to be included in the analysis.
#' @param method the algorithm to identify the change point. "CCS", "PD", or "changepoint".
#' @param method_args more arguments passed down to the algorithm that locates recombination.
#' @note This function produces a S3 class "RecView".
#' @export
rec2gen <- function(data, scaffold_info, chromosome, offspring, method = "changepoint", method_args = NULL) {
  if(is.null(method_args)) {
    if (method == "CCS") method_args <- list(threshold = 50, se = FALSE)
    if (method == "PD") method_args <- list(window_size = 400, threshold = 0.8)
  } else {
    if (method == "CCS") {
      if (is.null(method_args$threshold)) method_args <- list(threshold = 50, se = as.logical(unlist(method_args)))
      if (is.null(method_args$se)) method_args <- list(threshold = as.vector(unlist(method_args)), se = FALSE)
    }
    
    if (method == "PD") {
      if (is.null(method_args$window_size)) method_args <- list(window_size = 400, threshold = as.vector(unlist(method_args)))
      if (is.null(method_args$threshold)) method_args <- list(threshold = 0.8, window_size = as.vector(unlist(method_args)))
      method_args <- list(window_size = method_args$window_size, threshold = method_args$threshold)
    }
  }
  pos_col <- ifelse(method == "CCS", ifelse(method_args$se, "Mean_bp", "Middle_bp"), "Position_bp")
  
  # functions ----
  calc_distm <- function(x) {
    x <- str_split_1(x, "_")
    dm <- suppressWarnings(as.matrix(dist(x)))
    return(dm)
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
  
  # for (i in offspring) {
  #   data_chr_mod <- data_chr_mod %>%
  #     filter(get(i, .) <= sumAB_CD, get(i, .) >= minAB_CD)
  # }
  
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
  
  infosite_tb <- tibble(Offspring = offspring) %>% 
    rowwise() %>% 
    mutate(Paternal = length(which(!is.na(pull(data_chr_mod, Offspring)))))
  
  data_chr_mod <- data_chr_mod %>% 
    select(-diffAB_CD, -sumAB_CD, -minAB_CD) %>% 
    select(id, CHROM, orientation, POS, POS_chr, all_of(offspring)) %>% 
    unite(col = "GTS", all_of(offspring), sep = "_", remove = FALSE) %>% 
    nest(input = !GTS) %>% 
    rowwise() %>% 
    mutate(dm = list(calc_distm(GTS)))
  
  if (length(offspring) > 2) {
    dms <- data_chr_mod$dm
    pairwise_tbs <- list()
    for (index in 1:(length(offspring)-1)) {
      off_index <- seq(index+1, length(offspring))
      
      pairwise_tb <- tibble(group = as.character())
      for (i in off_index) {
        if (i > index) {
          pairwise_tb <- add_row(pairwise_tb, group = paste(index, i, sep = "_<>_"))
        } else {
          pairwise_tb <- add_row(pairwise_tb, group = paste(i, index, sep = "_<>_"))
        }
      }
      
      for(i in 1:length(dms)) {
        column <- dms[[i]][,index] %>% 
          .[off_index] %>% 
          as_tibble() %>% 
          `colnames<-`(., paste0("m", i))
        pairwise_tb <- bind_cols(pairwise_tb, column)
      }
      
      pairwise_tbs[[length(pairwise_tbs)+1]] <- pairwise_tb %>% 
        column_to_rownames(var = "group") %>% 
        t() %>% 
        as_tibble()
    }
    pairwise_tbs <- reduce(pairwise_tbs, bind_cols)
    
    pairwise_tb <- bind_cols(select(data_chr_mod, input), pairwise_tbs) %>% 
      unnest(cols = "input") %>% 
      arrange(POS_chr)
    
    rec_num_Q <- tibble(comparison = colnames(pairwise_tbs))
    if (method == "changepoint") {
      rec_num_Q <- rec_num_Q %>% 
        rowwise() %>% 
        mutate(n = length(changepoint::cpts(changepoint::cpt.mean(pull(pairwise_tb %>% filter(!is.na(pull(pairwise_tb, comparison))), comparison), minseglen = 1))))
    } else if (method == "PD") {
      rec_num_Q <- rec_num_Q %>% 
        mutate(id = seq(1,nrow(.))) %>% 
        nest(input = !id) %>% 
        mutate(n = map(input, function(tb) length(PD_algorithm(x = pairwise_tb %>% filter(!is.na(pull(pairwise_tb, tb$comparison))), value_col = tb$comparison, window_size = ifelse(is.null(method_args$window_size), 400, method_args$window_size), threshold = ifelse(is.null(method_args$threshold), 0.8, method_args$threshold), full_result = FALSE)))) %>% 
        select(input, n) %>% 
        unnest(cols = c("input", "n"))
    } else if (method == "CCS") {
      rec_num_Q <- rec_num_Q %>% 
        mutate(id = seq(1,nrow(.))) %>% 
        nest(input = !id) %>% 
        mutate(n = map(input, function(tb) nrow(CCS_algorithm(x = pairwise_tb %>% filter(!is.na(pull(pairwise_tb, tb$comparison))), value_col = tb$comparison, threshold = ifelse(is.null(method_args$threshold), 50, method_args$threshold), full_result = FALSE)))) %>% 
        select(input, n) %>% 
        unnest(cols = c("input", "n"))
    }
    
    off_no_rec <- tibble(Offspring = offspring, n = 0)
    for (i in 1:(length(offspring)-2)) {
      for (j in (i+1):(length(offspring)-1)) {
        for (k in (j+1):length(offspring)) {
          D <-  c(rec_num_Q %>% filter(comparison == paste(i, j, sep = "_<>_")) %>% pull(n), 
                  rec_num_Q %>% filter(comparison == paste(i, k, sep = "_<>_")) %>% pull(n),
                  rec_num_Q %>% filter(comparison == paste(j, k, sep = "_<>_")) %>% pull(n))
          
          comp_index <- list(c(i,j), c(i,k), c(j,k))
          Dmax_index <- which(D == max(D))
          if (length(Dmax_index) == 1) {
            if(max(D) == sum(D[which(D != max(D))])) {
              off_no_rec$n[off_no_rec$Offspring == offspring[c(i,j,k)[which(!(c(i,j,k) %in% comp_index[[Dmax_index]]))]]] <- off_no_rec$n[off_no_rec$Offspring == offspring[c(i,j,k)[which(!(c(i,j,k) %in% comp_index[[Dmax_index]]))]]] + 1
              # off_no_rec <- unique(c(off_no_rec, offspring[c(i,j,k)[which(!(c(i,j,k) %in% comp_index[[Dmax_index]]))]]))
            } else if (length(which(D == 0)) > 0) {
              off_no_rec$n[off_no_rec$Offspring == offspring[unique(unlist(comp_index[which(D == 0)]))]] <- off_no_rec$n[off_no_rec$Offspring == offspring[unique(unlist(comp_index[which(D == 0)]))]] + 1
              # off_no_rec <- unique(c(off_no_rec, offspring[unique(unlist(comp_index[which(D == 0)]))]))
            }
          }
        }
      }
    }
  }
  off_no_rec_pat <- off_no_rec %>% 
    filter(n == max(n), n > 0) %>% 
    # filter(n == (length(offspring)-1)*(length(offspring)-2)/2) %>% 
    pull(Offspring)
  
  if (length(off_no_rec_pat) != 0) {
    cat(paste0("\n- Note: Found ", length(off_no_rec_pat) ," individual(s) without recombination on paternal chromosome."))
    
    res_pat <- tibble(Offspring = offspring) %>% 
      filter(!(Offspring %in% off_no_rec_pat)) %>% 
      rowwise() %>% 
      mutate(output = list(RecView:::rec2gen_internal(data = data, scaffold_info = scaffold_info, chromosome = chromosome, offspring = c(Offspring, off_no_rec_pat), side = "Paternal", method = method, method_args = method_args, verbose = FALSE)))
    
    position_result_pat <- list()
    comparison_result_pat <- list()
    for (i in 1:nrow(res_pat)) {
      position_result_pat[[length(position_result_pat)+1]] <- res_pat$output[[i]]$recombination_position %>% mutate(Offspring = res_pat$Offspring[i])
      comparison_result_pat[[length(comparison_result_pat)+1]] <- res_pat$output[[i]]$pairwise_comparison
    }
    position_result_pat <- reduce(position_result_pat, bind_rows)
    comparison_result_pat <- reduce(comparison_result_pat, bind_rows) %>% 
      distinct(id, CHROM, orientation, POS, POS_chr, group, Concord_Discord, Chromosome_origin)
  } else {
    res_pat <- RecView:::rec2gen_internal(data = data, scaffold_info = scaffold_info, chromosome = chromosome, offspring = offspring, side = "Paternal", method = method, method_args = method_args, verbose = FALSE)
    position_result_pat <- res_pat$recombination_position
    comparison_result_pat <- res_pat$pairwise_comparison
  }
  
  dup <- position_result_pat %>% 
      group_by_at(pos_col) %>% 
      tally() %>% 
      filter(n > 1) %>% 
      get(pos_col, .)
  
  position_result_pat <- position_result_pat %>% filter(!(get(pos_col, .) %in% dup))
  
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
  
  # for (i in offspring) {
  #   data_chr_mod <- data_chr_mod %>%
  #     filter(get(i, .) <= sumAB_CD, get(i, .) >= minAB_CD)
  # }
  
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
  
  infosite_tb <- infosite_tb %>% 
    rowwise() %>% 
    mutate(Maternal = length(which(!is.na(pull(data_chr_mod, Offspring))))) %>% 
    ungroup()
  
  data_chr_mod <- data_chr_mod %>% 
    select(-diffAB_CD, -sumAB_CD, -minAB_CD) %>% 
    select(id, CHROM, orientation, POS, POS_chr, all_of(offspring)) %>% 
    unite(col = "GTS", all_of(offspring), sep = "_", remove = FALSE) %>% 
    nest(input = !GTS) %>% 
    rowwise() %>% 
    mutate(dm = list(calc_distm(GTS)))
  
  if (length(offspring) > 2) {
    dms <- data_chr_mod$dm
    pairwise_tbs <- list()
    for (index in 1:(length(offspring)-1)) {
      off_index <- seq(index+1, length(offspring))
      
      pairwise_tb <- tibble(group = as.character())
      for (i in off_index) {
        if (i > index) {
          pairwise_tb <- add_row(pairwise_tb, group = paste(index, i, sep = "_<>_"))
        } else {
          pairwise_tb <- add_row(pairwise_tb, group = paste(i, index, sep = "_<>_"))
        }
      }
      
      for(i in 1:length(dms)) {
        column <- dms[[i]][,index] %>% 
          .[off_index] %>% 
          as_tibble() %>% 
          `colnames<-`(., paste0("m", i))
        pairwise_tb <- bind_cols(pairwise_tb, column)
      }
      
      pairwise_tbs[[length(pairwise_tbs)+1]] <- pairwise_tb %>% 
        column_to_rownames(var = "group") %>% 
        t() %>% 
        as_tibble()
    }
    pairwise_tbs <- reduce(pairwise_tbs, bind_cols)
    
    pairwise_tb <- bind_cols(select(data_chr_mod, input), pairwise_tbs) %>% 
      unnest(cols = "input") %>% 
      arrange(POS_chr)
    
    rec_num_Q <- tibble(comparison = colnames(pairwise_tbs))
    if (method == "changepoint") {
      rec_num_Q <- rec_num_Q %>% 
        rowwise() %>% 
        mutate(n = length(changepoint::cpts(changepoint::cpt.mean(pull(pairwise_tb %>% filter(!is.na(pull(pairwise_tb, comparison))), comparison), minseglen = 1))))
    } else if (method == "PD") {
      rec_num_Q <- rec_num_Q %>% 
        mutate(id = seq(1,nrow(.))) %>% 
        nest(input = !id) %>% 
        mutate(n = map(input, function(tb) length(PD_algorithm(x = pairwise_tb %>% filter(!is.na(pull(pairwise_tb, tb$comparison))), value_col = tb$comparison, window_size = ifelse(is.null(method_args$window_size), 400, method_args$window_size), threshold = ifelse(is.null(method_args$threshold), 0.8, method_args$threshold), full_result = FALSE)))) %>% 
        select(input, n) %>% 
        unnest(cols = c("input", "n"))
    } else if (method == "CCS") {
      rec_num_Q <- rec_num_Q %>% 
        mutate(id = seq(1,nrow(.))) %>% 
        nest(input = !id) %>% 
        mutate(n = map(input, function(tb) nrow(CCS_algorithm(x = pairwise_tb %>% filter(!is.na(pull(pairwise_tb, tb$comparison))), value_col = tb$comparison, threshold = ifelse(is.null(method_args$threshold), 50, method_args$threshold), full_result = FALSE)))) %>% 
        select(input, n) %>% 
        unnest(cols = c("input", "n"))
    }
    
    off_no_rec <- tibble(Offspring = offspring, n = 0)
    for (i in 1:(length(offspring)-2)) {
      for (j in (i+1):(length(offspring)-1)) {
        for (k in (j+1):length(offspring)) {
          D <-  c(rec_num_Q %>% filter(comparison == paste(i, j, sep = "_<>_")) %>% pull(n), 
                  rec_num_Q %>% filter(comparison == paste(i, k, sep = "_<>_")) %>% pull(n),
                  rec_num_Q %>% filter(comparison == paste(j, k, sep = "_<>_")) %>% pull(n))
          
          comp_index <- list(c(i,j), c(i,k), c(j,k))
          Dmax_index <- which(D == max(D))
          if (length(Dmax_index) == 1) {
            if(max(D) == sum(D[which(D != max(D))])) {
              off_no_rec$n[off_no_rec$Offspring == offspring[c(i,j,k)[which(!(c(i,j,k) %in% comp_index[[Dmax_index]]))]]] <- off_no_rec$n[off_no_rec$Offspring == offspring[c(i,j,k)[which(!(c(i,j,k) %in% comp_index[[Dmax_index]]))]]] + 1
              # off_no_rec <- unique(c(off_no_rec, offspring[c(i,j,k)[which(!(c(i,j,k) %in% comp_index[[Dmax_index]]))]]))
            } else if (length(which(D == 0)) > 0) {
              off_no_rec$n[off_no_rec$Offspring == offspring[unique(unlist(comp_index[which(D == 0)]))]] <- off_no_rec$n[off_no_rec$Offspring == offspring[unique(unlist(comp_index[which(D == 0)]))]] + 1
              # off_no_rec <- unique(c(off_no_rec, offspring[unique(unlist(comp_index[which(D == 0)]))]))
            }
          }
        }
      }
    }
  }
  off_no_rec_mat <- off_no_rec %>% 
    filter(n == max(n), n > 0) %>% 
    # filter(n == (length(offspring)-1)*(length(offspring)-2)/2) %>% 
    pull(Offspring)  
  
  if (length(off_no_rec_mat) != 0) {
    cat(paste0("\n- Note: Found ", length(off_no_rec_mat) ," individuals without recombination on maternal chromosome."))
    
    res_mat <- tibble(Offspring = offspring) %>% 
      filter(!(Offspring %in% off_no_rec_mat)) %>% 
      rowwise() %>% 
      mutate(output = list(RecView:::rec2gen_internal(data = data, scaffold_info = scaffold_info, chromosome = chromosome, offspring = c(Offspring, off_no_rec_mat), side = "Maternal", method = method, method_args = method_args, verbose = FALSE)))
    
    position_result_mat <- list()
    comparison_result_mat <- list()
    for (i in 1:nrow(res_mat)) {
      position_result_mat[[length(position_result_mat)+1]] <- res_mat$output[[i]]$recombination_position %>% mutate(Offspring = res_mat$Offspring[i])
      comparison_result_mat[[length(comparison_result_mat)+1]] <- res_mat$output[[i]]$pairwise_comparison
    }
    position_result_mat <- reduce(position_result_mat, bind_rows)
    comparison_result_mat <- reduce(comparison_result_mat, bind_rows) %>% 
      distinct(id, CHROM, orientation, POS, POS_chr, group, Concord_Discord, Chromosome_origin)
  } else {
    res_mat <- RecView:::rec2gen_internal(data = data, scaffold_info = scaffold_info, chromosome = chromosome, offspring = offspring, side = "Maternal", method = method, method_args = method_args, verbose = FALSE)
    position_result_mat <- res_mat$recombination_position
    comparison_result_mat <- res_mat$pairwise_comparison
  }
  
  dup <- position_result_mat %>% 
    group_by_at(pos_col) %>% 
    tally() %>% 
    filter(n > 1) %>% 
    get(pos_col, .)
  
  position_result_mat <- position_result_mat %>% filter(!(get(pos_col, .) %in% dup))
  
  end_time <- Sys.time()
  msg <- paste0('\nStep 3/4: Maternal chromosome analysis completed. Used ', round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "secs")), 3), " sec.")
  cat(msg)
  
  # Results collecting ----
  start_time <- Sys.time()
  position_result <- bind_rows(position_result_pat, position_result_mat)
  position_result$Chromosome_origin <- factor(position_result$Chromosome_origin, levels = c("Paternal", "Maternal"))
  position_result <- position_result %>% 
    select(-ID) %>% 
    arrange_at(c(if (length(offspring) > 2) {"Offspring"}, "Chromosome_origin", pos_col)) %>% 
    # nest(input = all_of(pos_col)) %>% 
    # mutate(output = map(input, function(tb) tibble(ID = seq(1,nrow(tb))))) %>% 
    # unnest(cols = c("input", "output")) %>% 
    select(all_of(if (length(offspring) > 2) {c("Offspring", "Chromosome_origin")} else {"Chromosome_origin"}), all_of(if (pos_col == "Mean_bp") {c(pos_col, if (length(offspring) > 2) "SE_bp")} else {pos_col}), Precision, Precision_index)
  
  comparison_result <- bind_rows(comparison_result_pat, comparison_result_mat) %>% 
    distinct(id, CHROM, orientation, POS, POS_chr, group, Concord_Discord, Chromosome_origin)
  
  end_time <- Sys.time()
  msg <- paste0('\nStep 4/4: Results collected. Used ', round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "secs")), 3), " sec.")
  cat(msg)
  
  out <- list(chromosome = chromosome,
              method = method,
              method_args = method_args,
              offspring = left_join(tibble(Offspring = offspring), tibble(Offspring = off_no_rec_pat, Paternal_chromosome = "No_recombination"), by = "Offspring") %>% 
                left_join(tibble(Offspring = off_no_rec_mat, Maternal_chromosome = "No_recombination"), by = "Offspring") %>% 
                replace_na(list(Paternal_chromosome = "", Maternal_chromosome = "")), 
              informative_site = infosite_tb,
              pairwise_comparison = comparison_result, 
              recombination_position = position_result)
  
  cat('\nDone!')
  
  attr(out, "class") <- "RecView"
  return(out)
}
#' 
#' @noRd
rec2gen_internal <- function(data, scaffold_info, chromosome, offspring, side = c("Paternal", "Maternal"), method = "changepoint", method_args = NULL, verbose = TRUE) {
  if(is.null(method_args)) {
    if (method == "CCS") method_args <- list(threshold = 50, se = FALSE)
    if (method == "PD") method_args <- list(window_size = 400, threshold = 0.8)
  } else {
    if (method == "CCS") {
      if (is.null(method_args$threshold)) method_args <- list(threshold = 50, se = as.logical(unlist(method_args)))
      if (is.null(method_args$se)) method_args <- list(threshold = as.vector(unlist(method_args)), se = FALSE)
    }
    
    if (method == "PD") {
      if (is.null(method_args$window_size)) method_args <- list(window_size = 400, threshold = as.vector(unlist(method_args)))
      if (is.null(method_args$threshold)) method_args <- list(threshold = 0.8, window_size = as.vector(unlist(method_args)))
      method_args <- list(window_size = method_args$window_size, threshold = method_args$threshold)
    }
  }
  pos_col <- ifelse(method == "CCS", ifelse(method_args$se, "Mean_bp", "Middle_bp"), "Position_bp")
  
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
  
  get_cp <- function(input, data_in, offspring_in, method, method_args, pos_col) {
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
        cpt <- PD_algorithm(x = data_in_mod_tmp, value_col = i, window_size = ifelse(is.null(method_args$window_size), 400, method_args$window_size), threshold = ifelse(is.null(method_args$threshold), 0.8, method_args$threshold), full_result = FALSE)
        cps[[length(cps)+1]] <- tibble(cps_id = length(cps)+1,
                                       start_POS = data_in_mod_tmp$POS_chr[cpt],
                                       end_POS = data_in_mod_tmp$POS_chr[cpt])
      }
    }
    
    
    position_result <- reduce(cps, bind_rows)
    
    if (nrow(position_result) == 0) {
      position_result_pass <- tibble(ID = as.numeric(), !!pos_col := as.numeric())
    } else {
      position_result <- position_result %>% mutate(ID = seq(1, nrow(.)))
      if (length(cps) == 1) {
        position_result_pass <- position_result %>% 
          mutate(ID = seq(1,nrow(.)), !!pos_col := as.integer((start_POS + end_POS)/2)) %>% 
          select(ID, all_of(pos_col)) %>% 
          arrange_at(c("ID", pos_col))
        
        interval <- 1e4
        if (max(data_in_mod$POS_chr) / interval < 100) interval <- max(data_in_mod$POS_chr)/100
        digits <- floor(log10(interval)) + 1
        position_density_index <- which(hist(get(pos_col, position_result_pass), breaks = seq(0, ceiling(max(data_in_mod$POS_chr/10^digits)) * 10^digits, interval), plot = FALSE)$count != 0)
        infosite_density <- hist(data_in_mod$POS_chr, breaks = seq(0, ceiling(max(data_in_mod$POS_chr/10^digits)) * 10^digits, interval), plot = FALSE)$count
        position_result_pass <- position_result_pass %>% 
          mutate(Precision = round(interval / infosite_density[position_density_index], digits = 3), 
                 Precision_index = round((max(data_in_mod$POS_chr) / sum(infosite_density) - interval / infosite_density[position_density_index]) / (max(data_in_mod$POS_chr) / sum(infosite_density)), digits = 1))
      } else {
        position_result <- position_result %>% mutate(n = 0, index = "")
        
        overlaps_matrix <- matrix(rep("", times = nrow(position_result)^2), nrow = nrow(position_result))
        for (i in 1:nrow(position_result)) {
          for (j in i:nrow(position_result)) {
            rangeA <- c(position_result$start_POS[i], position_result$end_POS[i])
            rangeB <- c(position_result$start_POS[j], position_result$end_POS[j])
            overlaps_matrix[j,i] <- check_overlaps(x = rangeA, y = rangeB)
          }
        }
        
        all_i <- vector()
        for (i in 1:nrow(position_result)) {
          if (!( i %in% all_i)) {
            position_result$n[i] <- length(which(overlaps_matrix[i:nrow(position_result),i] == TRUE))
            position_result$index[i] <- str_c(which(overlaps_matrix[i:nrow(position_result),i] == TRUE)+i-1, collapse = ",")
          }
          all_i <- unique(c(all_i, which(overlaps_matrix[i:nrow(position_result),i] == "TRUE")+i-1))
        }
        n_threshold <- length(cps) / 2
        if (!is.integer(n_threshold)) n_threshold <- ceiling(length(cps)/2)
        position_result_pass <- position_result %>% filter(n > n_threshold)  # occurrence must satisfy a certain condition
        
        if (nrow(position_result_pass) > 0) {
          if (pos_col == "Mean_bp") {
            position_result_pass <- position_result_pass %>% 
              mutate(Mean_bp = as.numeric(NA), SE_bp = as.numeric(NA))
            
            for (i in 1:nrow(position_result_pass)) {
              index <- as.numeric(str_split_1(position_result_pass$index[i], pattern = ","))
              mid_POS <- vector()
              for (j in 1:1e3) {
                index_sample <- sample(index, size = length(index), replace = TRUE)
                start_POS <- max(position_result$start_POS[index_sample])
                end_POS <- min(position_result$end_POS[index_sample])
                mid_POS[length(mid_POS)+1] <- (start_POS + end_POS)/2
              }
              position_result_pass$Mean_bp[i] <- as.integer(mean(mid_POS))
              position_result_pass$SE_bp[i] <- round(sd(mid_POS) / sqrt(1e3), digits = 3)
            }
            position_result_pass <- position_result_pass %>% 
              arrange(Mean_bp) %>% 
              mutate(ID = seq(1, nrow(.))) %>% 
              select(ID, Mean_bp, SE_bp)
            
            interval <- 1e4
            if (max(data_in_mod$POS_chr) / interval < 100) interval <- max(data_in_mod$POS_chr)/100
            digits <- floor(log10(interval)) + 1
            position_density_index <- which(hist(get(pos_col, position_result_pass), breaks = seq(0, ceiling(max(data_in_mod$POS_chr/10^digits)) * 10^digits, interval), plot = FALSE)$count != 0)
            infosite_density <- hist(data_in_mod$POS_chr, breaks = seq(0, ceiling(max(data_in_mod$POS_chr/10^digits)) * 10^digits, interval), plot = FALSE)$count
            position_result_pass <- position_result_pass %>% 
              mutate(Precision = round(interval / infosite_density[position_density_index], digits = 3), 
                     Precision_index = round((max(data_in_mod$POS_chr) / sum(infosite_density) - interval / infosite_density[position_density_index]) / (max(data_in_mod$POS_chr) / sum(infosite_density)), digits = 1))
          } else {
            for (i in 1:nrow(position_result_pass)) {
              index <- as.numeric(str_split_1(position_result_pass$index[i], pattern = ","))
              position_result_pass$start_POS[i] <- max(position_result$start_POS[index])
              position_result_pass$end_POS[i] <- min(position_result$end_POS[index])
            }
            
            position_result_pass <- position_result_pass %>% 
              mutate(!!pos_col := as.integer((start_POS + end_POS)/2)) %>% 
              arrange_at(pos_col) %>% 
              mutate(ID = seq(1, nrow(.))) %>% 
              select(ID, all_of(pos_col))
            
            interval <- 1e4
            if (max(data_in_mod$POS_chr) / interval < 100) interval <- max(data_in_mod$POS_chr)/100
            digits <- floor(log10(interval)) + 1
            position_density_index <- vector()
            for (i in 1:nrow(position_result_pass)) {
              position_density_index[length(position_density_index)+1] <- which(hist(get(pos_col, position_result_pass[i,]), breaks = seq(0, ceiling(max(data_in_mod$POS_chr/10^digits)) * 10^digits, interval), plot = FALSE)$count != 0)
            }
            infosite_density <- hist(data_in_mod$POS_chr, breaks = seq(0, ceiling(max(data_in_mod$POS_chr/10^digits)) * 10^digits, interval), plot = FALSE)$count
            position_result_pass <- position_result_pass %>% 
              mutate(Precision = round(interval / infosite_density[position_density_index], digits = 3), 
                     Precision_index = round((max(data_in_mod$POS_chr) / sum(infosite_density) - interval / infosite_density[position_density_index]) / (max(data_in_mod$POS_chr) / sum(infosite_density)), digits = 1))
          }
        } else {
          position_result_pass <- tibble(ID = as.numeric(), !!pos_col := as.numeric(), Precision = as.numeric(), Precision_index = as.numeric())
        }
      }
    }
    
    comparison_result <- data_in_mod %>% 
      gather(key = "group", value = "Concord_Discord", colnames(pairwise_tb))
    
    out <- tibble(comparison_result = list(comparison_result), position_result = list(position_result_pass))
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
  if (verbose) cat(msg)
  
  # Paternal ----
  start_time <- Sys.time()
  data_chr_mod <- data_chr %>% 
    select(id, CHROM, orientation, POS, POS_chr, AB, CD, all_of(offspring)) %>% 
    mutate(diffAB_CD = abs(AB - CD), sumAB_CD = AB + CD) %>% 
    filter(diffAB_CD == 1, AB == 1) %>% 
    rowwise() %>% 
    mutate(minAB_CD = min(AB, CD)) %>% 
    ungroup()
  
  # for (i in offspring) {
  #   data_chr_mod <- data_chr_mod %>%
  #     filter(get(i, .) <= sumAB_CD, get(i, .) >= minAB_CD)
  # }
  
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
  
  if ("Paternal" %in% side) {
    res_pat <- tibble(Offspring = offspring) %>% 
      mutate(index = seq(1, nrow(.))) %>% 
      nest(input = !Offspring) %>% 
      mutate(output = map(input, get_cp, data_chr_mod, offspring, method, method_args, pos_col)) %>% 
      unnest(cols = "output") %>% 
      select(-input) %>% 
      mutate(Note = '-', Chromosome_origin = "Paternal")
    
    if (nrow(res_pat) == 2) {
      position_result_pat <- res_pat %>% 
        slice(1) %>% 
        select(-comparison_result, -Offspring) %>% 
        unnest(cols = "position_result") %>% 
        select(Chromosome_origin, everything())
    } else {
      position_result_pat <- res_pat %>% 
        select(-comparison_result) %>% 
        unnest(cols = "position_result") %>% 
        select(Chromosome_origin, Offspring, everything())
      
      dup <- position_result_pat %>% 
        group_by_at(pos_col) %>% 
        tally() %>% 
        filter(n > 1) %>% 
        get(pos_col, .)
      
      position_result_pat$Note[which((get(pos_col, position_result_pat) %in% dup) == TRUE)] <- "Duplication"
    }
    
    comparison_result_pat <- res_pat %>% 
      select(comparison_result, Chromosome_origin) %>% 
      unnest(cols = "comparison_result")
    
    end_time <- Sys.time()
    msg <- paste0('\nStep 2/4: Paternal chromosome analysis completed. Used ', round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "secs")), 3), " sec.")
    if (verbose) cat(msg)
  } else { 
    position_result_pat <- tibble(Chromosome_origin = as.character(),
                                  Offspring = as.character(),
                                  ID = as.numeric(),
                                  !!pos_col := as.numeric(),
                                  Precision = as.numeric(),
                                  Precision_index = as.numeric(),
                                  Note = as.character())
    
    comparison_result_pat <- data_chr_mod$input[[1]] %>% 
      slice(1) %>% 
      mutate(Concord_Discord = 1, Chromosome_origin = "Paternal") %>% 
      filter(id == "None")
  }
  
  # Maternal ----
  start_time <- Sys.time()
  data_chr_mod <- data_chr %>% 
    select(id, CHROM, orientation, POS, POS_chr, AB, CD, all_of(offspring)) %>% 
    mutate(diffAB_CD = abs(AB - CD), sumAB_CD = AB + CD) %>% 
    filter(diffAB_CD == 1, CD == 1) %>% 
    rowwise() %>% 
    mutate(minAB_CD = min(AB, CD)) %>% 
    ungroup()
  
  # for (i in offspring) {
  #   data_chr_mod <- data_chr_mod %>%
  #     filter(get(i, .) <= sumAB_CD, get(i, .) >= minAB_CD)
  # }
  
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
  
  if ("Maternal" %in% side) {
    res_mat <- tibble(Offspring = offspring) %>% 
      mutate(index = seq(1, nrow(.))) %>% 
      nest(input = !Offspring) %>% 
      mutate(output = map(input, get_cp, data_chr_mod, offspring, method, method_args, pos_col)) %>% 
      unnest(cols = "output") %>% 
      select(-input) %>% 
      mutate(Note = '-', Chromosome_origin = "Maternal")
    
    if (nrow(res_mat) == 2) {
      position_result_mat <- res_mat %>% 
        slice(1) %>% 
        select(-comparison_result, -Offspring) %>% 
        unnest(cols = "position_result") %>% 
        select(Chromosome_origin, everything())
    } else {
      position_result_mat <- res_mat %>% 
        select(-comparison_result) %>% 
        unnest(cols = "position_result") %>% 
        select(Chromosome_origin, Offspring, everything())
      
      dup <- position_result_mat %>% 
        group_by_at(pos_col) %>% 
        tally() %>% 
        filter(n > 1) %>% 
        get(pos_col, .)
      
      position_result_mat$Note[which((get(pos_col, position_result_mat) %in% dup) == TRUE)] <- "Duplication"
    }
    
    comparison_result_mat <- res_mat %>% 
      select(comparison_result, Chromosome_origin) %>% 
      unnest(cols = "comparison_result")
    
    end_time <- Sys.time()
    msg <- paste0('\nStep 3/4: Maternal chromosome analysis completed. Used ', round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "secs")), 3), " sec.")
    if (verbose) cat(msg)
  } else { 
    position_result_mat <- tibble(Chromosome_origin = as.character(),
                                  Offspring = as.character(),
                                  ID = as.numeric(),
                                  !!pos_col := as.numeric(),
                                  Precision = as.numeric(),
                                  Precision_index = as.numeric(),
                                  Note = as.character())
    
    comparison_result_mat <- data_chr_mod$input[[1]] %>% 
      slice(1) %>% 
      mutate(Concord_Discord = 1, Chromosome_origin = "Maternal") %>% 
      filter(id == "None")
  }
  
  # Results collecting ----
  start_time <- Sys.time()
  position_result <- bind_rows(position_result_pat, position_result_mat)
  position_result$Chromosome_origin <- factor(position_result$Chromosome_origin, levels = c("Paternal", "Maternal"))
  if (nrow(position_result) > 0) {
    position_result <- position_result %>% 
      filter(Note != "Duplication") %>% 
      select(-Note, -ID) %>% 
      arrange_at(c(if (length(offspring) > 2) {"Offspring"}, "Chromosome_origin", pos_col)) %>% 
      nest(input = all_of(pos_col)) %>% 
      mutate(output = map(input, function(tb) tibble(ID = seq(1,nrow(tb))))) %>% 
      unnest(cols = c("input", "output")) %>% 
      select(all_of(if (length(offspring) > 2) {c("Offspring", "Chromosome_origin")} else {"Chromosome_origin"}), ID, all_of(if (pos_col == "Mean_bp") {c(pos_col, if (length(offspring) > 2) "SE_bp")} else {pos_col}), Precision, Precision_index)
  } else {
    position_result <- position_result %>% 
      select(-Note) %>% 
      select(all_of(if (length(offspring) > 2) {c("Offspring", "Chromosome_origin")} else {"Chromosome_origin"}), ID, all_of(if (pos_col == "Mean_bp") {c(pos_col, if (length(offspring) > 2) "SE_bp")} else {pos_col}), Precision, Precision_index)
  }
  comparison_result <- bind_rows(comparison_result_pat, comparison_result_mat) 
  
  end_time <- Sys.time()
  msg <- paste0('\nStep 4/4: Results collected. Used ', round(as.numeric(difftime(time1 = end_time, time2 = start_time, units = "secs")), 3), " sec.")
  if (verbose) cat(msg)
  
  out <- list(chromosome = chromosome,
              offspring = offspring, 
              method = method,
              method_args = method_args,
              pairwise_comparison = comparison_result, 
              recombination_position = position_result)
  
  if (verbose) cat('\nDone!')
  
  attr(out, "class") <- "RecView"
  return(out)
}
