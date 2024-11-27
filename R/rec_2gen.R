#' Locate recombination with two-generation pedigree
#' @description Locate recombination with two-generation pedigree that consists of father, mother and at least two offspring
#' @usage rec_2gen(data, sc_order, chromosome, offspring, precision)
#' @param data '012'-formatted genotype file, which is the output from make_012gt() or make_012gt_from_vcf().
#' @param sc_order scaffold file.
#' @param chromosome chromosome.
#' @param offspring the offspring to be included in the analysis.
#' @param precision the maximal base pairs in which the detected recombination positions will be pooled.
#' @param method the algorithm to identify the change point. "CCS", "PD", or "changepoint".
#' @note This function produces a S3 class "RecView".
#' 
#' @export
rec_2gen <- function(data, sc_order, chromosome, offspring, precision = 1e4, method = "changepoint", ...) {
  # functions ----
  calc_distm <- function(x) {
    x <- str_split_1(x, "_")
    dm <- as.matrix(dist(x))
    return(dm)
  }
  
  get_cp <- function(input, data_in, offspring_in, precision, method, ...) {
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
        cpt <- changepoint::cpt.mean(pull(data_in_mod, i), minseglen = 1, ...)
        cps[[length(cps)+1]] <- data_in_mod$POS_chr[changepoint::cpts(cpt)]
      }
    } else if (method == "CCS") {
      for (i in colnames(pairwise_tb)) {
        cps[[length(cps)+1]] <- CCS_algorithm(x = data_in_mod, value_col = i, threshold = 50, full_result = FALSE) %>% 
          mutate(start_POS = data_in_mod$POS_chr[start_row], end_POS = data_in_mod$POS_chr[end_row]) %>% 
          mutate(middle_POS = (start_POS + end_POS) / 2) %>% 
          pull(middle_POS)
      }
    } else if (method == "PD") {
      for (i in colnames(pairwise_tb)) {
        cpt <- PD_algorithm(x = data_in_mod, value_col = i, window_size = 1100, threshold = 0.8, full_result = FALSE)
        cps[[length(cps)+1]] <- data_in_mod$POS_chr[cpt]
      }
    }
    
    if (length(cps) == 1) {
      cps <- unlist(cps)
      position_result <- tibble(POS_chr = cps) %>% 
        mutate(ID = seq(1, nrow(.))) %>%
        arrange(POS_chr)
      if (length(cps) > 1) {
        for (i in 1:(nrow(position_result)-1)) {
          for (j in (i+1):nrow(position_result)) {
            if (abs(position_result$POS_chr[j] - position_result$POS_chr[i]) <= precision) {
              position_result$ID[j] <- position_result$ID[i]
            } else {
              position_result$ID[j] <- position_result$ID[i] + 1
            }
          }
        }
        
        position_result <- position_result %>% 
          summarise(Mean_bp = mean(POS_chr), SD_bp = sd(POS_chr), .by = "ID")
      } else {
        position_result <- position_result %>% 
          mutate(SD_bp = as.numeric(NA)) %>% 
          select(ID, POS_chr, SD_bp) %>% 
          `colnames<-`(.,c("ID", "Mean_bp", "SD_bp"))
        
      }
    } else if (length(cps) > 1) {
      cps <- reduce(cps, c)
      if (length(unique(cps[duplicated(cps)])) == 1) {
        position_result <- tibble(POS_chr = unique(cps[duplicated(cps)])) %>% 
          mutate(ID = 1, Mean_bp = POS_chr, SD_bp = as.numeric(NA)) %>% 
          select(ID, Mean_bp, SD_bp)
        
      } else if (length(unique(cps[duplicated(cps)])) > 1) {
        cps_mod <- tibble(POS_chr = unique(cps[duplicated(cps)])) %>% 
          arrange(POS_chr) %>% 
          mutate(ID = c(1, rep(0, times = nrow(.) - 1)))
        
        for (i in 1:(nrow(cps_mod)-1)) {
          for (j in (i+1):nrow(cps_mod)) {
            if (abs(cps_mod$POS_chr[j] - cps_mod$POS_chr[i]) <= precision) {
              cps_mod$ID[j] <- cps_mod$ID[i]
            } else {
              cps_mod$ID[j] <- cps_mod$ID[i] + 1
            }
          }
        }
        
        position_result <- cps_mod %>% 
          summarise(Mean_bp = mean(POS_chr), SD_bp = sd(POS_chr), .by = "ID")
      } else {
        position_result <- tibble(ID = as.numeric(), Mean_bp = as.numeric(), SD_bp = as.numeric())
      }
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
    data_chr_mod <- data_chr_mod %>% 
      filter(get(i, .) <= sumAB_CD, get(i, .) >= minAB_CD)
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
    mutate(output = map(input, get_cp, data_chr_mod, offspring, precision, method)) %>% 
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
    data_chr_mod <- data_chr_mod %>% 
      filter(get(i, .) <= sumAB_CD, get(i, .) >= minAB_CD)
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
    mutate(output = map(input, get_cp, data_chr_mod, offspring, precision, method)) %>% 
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
              precision = precision,
              method = method,
              params = ...)
  
  cat('\nDone!')
  
  attr(out, "class") <- "RecView"
  return(out)
}
