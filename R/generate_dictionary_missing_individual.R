#' Generate <Dictionary of grandparent-of-origin> for cases of missing genotypes of individuals
#'
#' @noRd
generate_dictionary_missing_individual <- function(dictionary, save_filename = "dict.tsv") {
  ind <- c("A", "B", "C", "D", "AB", "CD")
  dictionary_new <- tibble(Missing_ind_num = c(1, 2, 3, 4))
  
  comb_1 <- list() # Missing One Individual
  for (i in 1:length(ind)) {
    comb_1[[length(comb_1)+1]] <- ind[i]
  }
  
  comb_2 <- list() # Missing Two Individual
  for (i in 1:(length(ind)-1)) {
    for (j in (i+1):length(ind)) {
      comb_2[[length(comb_2)+1]] <- c(ind[i], ind[j])
    }
  }
  
  comb_3 <- list() # Missing Three Individual
  for (i in 1:(length(ind)-2)) {
    for (j in (i+1):(length(ind)-1)) {
      for(k in (j+1):(length(ind))) {
        comb_3[[length(comb_3)+1]] <- c(ind[i], ind[j], ind[k])
      }
    }
  }
  
  comb_4 <- list() # Missing Four Individual
  for (i in 1:(length(ind)-3)) {
    for (j in (i+1):(length(ind)-2)) {
      for (k in (j+1):(length(ind)-1)) {
        for(l in (k+1):(length(ind))) {
          comb_4[[length(comb_4)+1]] <- c(ind[i], ind[j], ind[k], ind[l])
        }
      }
    }
  }
  
  dictionary_new <- dictionary_new %>%
    mutate(comb = list(comb_1, comb_2, comb_3, comb_4)) %>%
    mutate(dict_new = map(comb, make_dictionary, dictionary))
  
  dictionary_to_save <- dictionary_new %>% 
    select(Missing_ind_num, dict_new) %>% 
    unnest(cols = "dict_new") %>% 
    select(GT_string, Paternal, Maternal, everything())
  
  write_tsv(dictionary_to_save, file = save_filename)
}
#'
#' Function to identify whether change of GT affects GoO inference
#' @noRd
make_Q <- function(dict_mod, side, combination) {
  Q <- tibble(X1 = c("GT_string", seq(1:length(combination))), X2 = c("Character", seq(1:length(combination)))) %>% 
    column_to_rownames(var = "X1") %>% 
    t() %>% 
    as_tibble() %>% 
    `colnames<-`(.,c("GT_string", lapply(combination, function(v) str_c(v, collapse = "_")) %>% unlist()))
  
  for (i in 1:nrow(dict_mod)) {
    if (i > 1) Q <- list(Q, Q[1,]) %>% reduce(bind_rows)
    Q$GT_string[i] <- dict_mod$GT_string[i]
    
    for (j in 1:length(combination)) {
      gt <- tibble(A=get("A", dict_mod)[i], 
                   B=get("B", dict_mod)[i], 
                   C=get("C", dict_mod)[i], 
                   D=get("D", dict_mod)[i], 
                   AB=get("AB", dict_mod)[i], 
                   CD=get("CD", dict_mod)[i], 
                   off=get("off", dict_mod)[i])
      goo <- get(side, dict_mod)[i]
      
      k_set <- expand.grid(rep(list(c(0,1,2)), times = length(combination[[j]]))) %>% 
        unite(col = "string", remove = F)
      
      gt_comb <- vector()
      for(x in 1:length(combination[[j]])) {
        gt_comb <- c(gt_comb, get(combination[[j]][x], gt))
      }
      gt_comb <- str_c(gt_comb, collapse = "_")
      
      k_set <- k_set[which(k_set$string != gt_comb),] %>% 
        as_tibble() %>% 
        select(-string) %>% 
        `colnames<-`(.,combination[[j]])
      
      gt_mod <- rep(list(gt), times = nrow(k_set)) %>% 
        reduce(bind_rows)
      
      for(col in colnames(k_set)) {
        gt_mod[col] <- k_set[col]
      }
      
      gt_mod <- unite(gt_mod, col = "GT_string_new", c("A", "B", "C", "D", "AB", "CD", "off"), sep = "_", remove = FALSE) %>% 
        mutate(score = 0, goo_ori = goo) %>% 
        left_join(dict_mod %>% 
                    select(GT_string, !!side) %>% 
                    `colnames<-`(.,c("GT_string", "goo_new")), 
                  by = c("GT_string_new" = "GT_string"))
      
      gt_mod$score <- ifelse(gt_mod$goo_new != gt_mod$goo_ori & !is.na(gt_mod$goo_new), 1, 0)
      
      if (sum(gt_mod$score) == 0) {
        Q[colnames(k_set) %>% str_c(collapse = "_")][i,] <- "T"  #T stands for Trivial
      } else {
        Q[colnames(k_set) %>% str_c(collapse = "_")][i,] <- "P"  #P stands for Pivotal
      }
      
    }
  }
  Q <- left_join(Q, dict_mod %>% select(GT_string, !!side, Chromosome_type), by = "GT_string")
  return(Q)
}
#'
#' Function to make dictionary based on chromosome type
#' @noRd
make_dict <- function(dict, chr_type, combination) {
  dict_mod <- dict %>% 
    filter(Chromosome_type == chr_type) %>% 
    mutate(id = seq(1:nrow(.)))
  
  Q_pat <- make_Q(dict_mod = dict_mod, side = "Paternal", combination = combination)
  Q_mat <- make_Q(dict_mod = dict_mod, side = "Maternal", combination = combination)
  
  dict_new <- list()
  for (i in 1:length(combination)) {
    dict_intermediate_pat <- Q_pat %>%
      select(GT_string, str_c(combination[[i]], collapse = "_"), Paternal, Chromosome_type) %>% 
      separate(col = "GT_string", into = c("A_", "B_", "C_", "D_", "AB_", "CD_", "off_"), sep = "_", remove = FALSE)
    dict_intermediate_mat <- Q_mat %>% 
      select(GT_string, str_c(combination[[i]], collapse = "_"), Maternal, Chromosome_type) %>% 
      separate(col = "GT_string", into = c("A_", "B_", "C_", "D_", "AB_", "CD_", "off_"), sep = "_", remove = FALSE)
    
    for(j in 1:length(combination[[i]])) {
      dict_intermediate_pat <- dict_intermediate_pat %>% 
        mutate("{str_c(combination[[i]], collapse = '_')}" := 
                 ifelse(get(str_c(combination[[i]], collapse = "_"), Q_pat) == "T" & get(str_c(combination[[i]], collapse = "_"), Q_mat) == "T", "T", "P")) %>% 
        mutate("{paste0(combination[[i]][j], '_')}" := 
                 ifelse(get(str_c(combination[[i]], collapse = "_"), .) == "T", "NA", get(paste0(combination[[i]][j], "_"), .)))
    }
    dict_intermediate_pat <- dict_intermediate_pat %>% 
      unite(col = "GT_string_new", c("A_", "B_", "C_", "D_", "AB_", "CD_", "off_"), sep = "_") %>% 
      left_join(dict_mod %>% select(id, GT_string), by = "GT_string") %>% 
      filter(get(str_c(combination[[i]], collapse = '_'), .) == "T") %>% 
      select(id, GT_string_new, Paternal, Chromosome_type)
    
    
    for(j in 1:length(combination[[i]])) {
      dict_intermediate_mat <- dict_intermediate_mat %>% 
        mutate("{str_c(combination[[i]], collapse = '_')}" := 
                 ifelse(get(str_c(combination[[i]], collapse = "_"), Q_mat) == "T" & get(str_c(combination[[i]], collapse = "_"), Q_pat) == "T", "T", "P")) %>%  
        mutate("{paste0(combination[[i]][j], '_')}" := 
                 ifelse(get(str_c(combination[[i]], collapse = "_"), .) == "T", "NA", get(paste0(combination[[i]][j], "_"), .)))
    }
    dict_intermediate_mat <- dict_intermediate_mat %>% 
      unite(col = "GT_string_new", c("A_", "B_", "C_", "D_", "AB_", "CD_", "off_"), sep = "_") %>% 
      left_join(dict_mod %>% select(id, GT_string), by = "GT_string") %>% 
      filter(get(str_c(combination[[i]], collapse = '_'), .) == "T") %>% 
      select(id, Maternal)
    
    dict_new[[length(dict_new)+1]] <- full_join(dict_intermediate_pat, dict_intermediate_mat, by = "id") %>% 
      select(GT_string_new, Paternal, Maternal) %>% 
      distinct(GT_string_new, Paternal, Maternal) %>% 
      mutate(Missing_ind = str_c(combination[[i]], collapse = "_"))
    
  }
  dict_new <- reduce(dict_new, bind_rows) %>% 
    mutate(Chromosome_type = chr_type)
  return(dict_new)
}
#'
#' Function to make dictionary
#' @noRd
make_dictionary <- function(combination_input, dictionary_input) {
  dict_ <- list()
  for (i in unique(dictionary_input$Chromosome_type)) {
    dict_[[i]] <- make_dict(dict = dictionary_input, chr_type = i, combination = combination_input)
  }
  dict_ <- reduce(dict_, bind_rows) %>% 
    `colnames<-`(.,c("GT_string", colnames(.)[2:length(colnames(.))]))
  return(dict_)
}
