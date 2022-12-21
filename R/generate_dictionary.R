#' Generate <Dictionary of grandparent-of-origin>
#'
#' @noRd
generate_dictionary <- function(chromosome = "autosome", missing = "none", save_filename = "dict.tsv") {
  pool <- c("0|0", "0|1", "1|0", "1|1")
  pool_012 <- c("0", "1", "1", "2")
  dict <- tibble(A = character(), B = character(), C = character(), D = character(), AB = character(), CD = character(), off = character())
  dict_012 <- tibble(A = character(), B = character(), C = character(), D = character(), AB = character(), CD = character(), off = character())

  for (i1 in 1:length(pool)) {
    for (i2 in 1:length(pool)) {
      for (i3 in 1:length(pool)) {
        for (i4 in 1:length(pool)) {
          for (i5 in 1:length(pool)) {
            for (i6 in 1:length(pool)) {
              for (i7 in 1:length(pool)) {
                dict[nrow(dict) + 1,] <- list(pool[i1], pool[i2], pool[i3], pool[i4], pool[i5], pool[i6], pool[i7])
                dict_012[nrow(dict_012) + 1,] <- list(pool_012[i1], pool_012[i2], pool_012[i3], pool_012[i4], pool_012[i5], pool_012[i6], pool_012[i7])
              }
            }
          }
        }
      }
    }
  }

  dict <- mutate(dict, index = seq(1:nrow(dict))) %>%
    select(index, everything()) %>%
    remove_impossible_cases(chromosome = chromosome, missing = missing)

  dict_012 <- mutate(dict_012, index = seq(1:nrow(dict_012))) %>%
    select(index, everything())

  # Assignment of grandparent-of-origin in the dictionary #
  left_sub <- dict %>% select(index, A, B, C, D, AB, CD, off) %>%
    mutate(off_left = str_sub(off, 1, 1)) %>%
    select(-off)

  right_sub <- dict %>% select(index, A, B, C, D, AB, CD, off) %>%
    mutate(off_right = str_sub(off, 3, 3)) %>%
    select(-off)

  ## Left side autosome
  left_sub$Origin <- rep('N', times = nrow(left_sub))
  if (missing != "pat") {
    if (chromosome != "X") {
      for (i in 1:nrow(left_sub)) {
        if (left_sub$off_left[i] == "0") {
          if (left_sub$AB[i] == "0|1") {
            if (left_sub$A[i] != "1|1" && left_sub$B[i] != "0|0") {
              left_sub$Origin[i] <- "A"
            }
          }
          if (left_sub$AB[i] == "1|0") {
            if (left_sub$A[i] != "0|0" && left_sub$B[i] != "1|1") {
              left_sub$Origin[i] <- "B"
            }
          }
        }

        if (left_sub$off_left[i] == "1") {
          if (left_sub$AB[i] == "0|1") {
            if (left_sub$A[i] != "1|1" && left_sub$B[i] != "0|0") {
              left_sub$Origin[i] <- "B"
            }
          }
          if (left_sub$AB[i] == "1|0") {
            if (left_sub$A[i] != "0|0" && left_sub$B[i] != "1|1") {
              left_sub$Origin[i] <- "A"
            }
          }
        }
      }
    }
  }
  ## X chromosome
  if (chromosome == "X") {
    for (i in 1:nrow(left_sub)) {
      if (left_sub$off_left[i] == "0") {
        if (left_sub$AB[i] == "0|0") {
          if (left_sub$B[i] != "1|1") {
            left_sub$Origin[i] <- "B"
          }
        }
      }
      if (left_sub$off_left[i] == "1") {
        if (left_sub$AB[i] == "1|1") {
          if (left_sub$B[i] != "0|0") {
            left_sub$Origin[i] <- "B"
          }
        }
      }
    }
  }

  ## Right side autosome
  right_sub$Origin <- rep('N', times = nrow(right_sub))
  if (missing != "mat") {
    if (chromosome != "Z") {
      for (i in 1:nrow(right_sub)) {
        if (right_sub$off_right[i] == "0") {
          if (right_sub$CD[i] == "0|1") {
            if (right_sub$C[i] != "1|1" && right_sub$D[i] != "0|0") {
              right_sub$Origin[i] <- "C"
            }
          }
          if (right_sub$CD[i] == "1|0") {
            if (right_sub$C[i] != "0|0" && right_sub$D[i] != "1|1") {
              right_sub$Origin[i] <- "D"
            }
          }
        }

        if (right_sub$off_right[i] == "1") {
          if (right_sub$CD[i] == "0|1") {
            if (right_sub$C[i] != "1|1" && right_sub$D[i] != "0|0") {
              right_sub$Origin[i] <- "D"
            }
          }
          if (right_sub$CD[i] == "1|0") {
            if (right_sub$C[i] != "0|0" && right_sub$D[i] != "1|1") {
              right_sub$Origin[i] <- "C"
            }
          }
        }
      }
    }
  }
  ## Z chromosome
  if (chromosome == "Z") {
    for (i in 1:nrow(right_sub)) {
      if (right_sub$off_right[i] == "0") {
        if (right_sub$CD[i] == "0|0") {
          if (right_sub$C[i] != "1|1") {
            right_sub$Origin[i] <- "C"
          }
        }
      }
      if (right_sub$off_right[i] == "1") {
        if (right_sub$CD[i] == "1|1") {
          if (right_sub$C[i] != "0|0") {
            right_sub$Origin[i] <- "C"
          }
        }
      }
    }
  }

  if (missing == "pat") {
    dict_012 <- dict_012 %>% select(-A, -B, -AB)
    columns <- c("C", "D", "CD", "off")
  } else if (missing == "mat") {
    dict_012 <- dict_012 %>% select(-C, -D, -CD)
    columns <- c("A", "B", "AB", "off")
  } else {
    columns <- c("A", "B", "C", "D", "AB", "CD", "off")
  }

  dict_mod <- dict %>%
    left_join(left_sub %>% select(index, Origin), by = "index") %>%
    left_join(right_sub %>% select(index, Origin), by = "index") %>%
    unite(col = "Origin", c("Origin.x", "Origin.y"), sep = '|') %>%
    select(index, Origin) %>%
    `colnames<-`(., c("index", "Origin")) %>%
    left_join(dict_012, by = "index") %>%
    unite(col = "GT_string", columns, sep = '_', remove = F) %>%
    nest(input = !GT_string) %>%
    mutate(output = map(input, dict_compact)) %>%
    select(-input) %>%
    unnest(cols = "output") %>%
    select(-index) %>%
    separate(col = "Origin", into = c("Paternal", "Maternal")) %>%
    mutate(Chromosome_type = chromosome)

  write_tsv(dict_mod, file = save_filename)
}

#' @noRd
remove_impossible_cases <- function(tb, chromosome, missing = "none") {
  tb <- tb %>% mutate(Possibility = "Y")

  for (i in 1:nrow(tb)) {
    if (missing != "pat") {
      A <- tb$A[i] %>% str_split('[|]') %>% unlist()
      B <- tb$B[i] %>% str_split('[|]') %>% unlist()
      AB <- tb$AB[i] %>% str_split('[|]') %>% unlist()

      if (!(str_sub(tb$AB[i],3,3) %in% B)) tb$Possibility[i] <- "N"
      if (!(str_sub(tb$off[i],1,1) %in% AB)) tb$Possibility[i] <- "N"
      if (chromosome != "X") {
        if (!(str_sub(tb$AB[i],1,1) %in% A)) tb$Possibility[i] <- "N"
      } else if (chromosome == "X") {
        if (sum(as.numeric(A)) == 1) tb$Possibility[i] <- "N"
        if (sum(as.numeric(AB)) == 1) tb$Possibility[i] <- "N"
      }
    }

    if (missing != "mat") {
      C <- tb$C[i] %>% str_split('[|]') %>% unlist()
      D <- tb$D[i] %>% str_split('[|]') %>% unlist()
      CD <- tb$CD[i] %>% str_split('[|]') %>% unlist()

      if (!(str_sub(tb$CD[i],1,1) %in% C)) tb$Possibility[i] <- "N"
      if (!(str_sub(tb$off[i],3,3) %in% CD)) tb$Possibility[i] <- "N"
      if (chromosome != "Z") {
        if (!(str_sub(tb$CD[i],3,3) %in% D)) tb$Possibility[i] <- "N"
      } else if (chromosome == "Z") {
        if (sum(as.numeric(D)) == 1) tb$Possibility[i] <- "N"
        if (sum(as.numeric(CD)) == 1) tb$Possibility[i] <- "N"
      }
    }
  }
  out <- tb %>% filter(Possibility == "Y") %>% select(-Possibility)
  return(out)
}

#' @noRd
dict_compact <- function(tb) {
  summary <- tb %>%
    group_by(Origin) %>%
    summarise(Count = length(index)) %>%
    arrange(Origin) %>%
    separate(col = "Origin", into = c("left", "right"), sep = '[|]')
  left <- summary$left %>% unique()
  right <- summary$right %>% unique()

  if (length(left) != 1) {
    left <- "N"
  }

  if (length(right) != 1) {
    right <- "N"
  }

  out <- tb[1,] %>%
    select(-Origin) %>%
    mutate(Origin = paste(left, right, sep = '|')) %>%
    select(index, Origin, everything())
  return(out)
}

