#' RecView server function
#' @noRd
recview_server <- function(input, output, session) {
  ## UI modification ----
  ### Banner ----
  output$title <- renderImage({
    rnum <- sample(x = c(1, 2, 3), size = 1)
    list(src = list.files(path = system.file("banner", package = "RecView"), pattern = paste0("banner_", rnum, ".png"), full.names = T),
         width = 383,
         height = 70)
  }, deleteFile = FALSE)
  ### Offspring input droplist ----
  output$genoexist <- reactive({
    if(!is.null(input$genofile)) {
      return(1)
    }
  })
  outputOptions(output, 'genoexist', suspendWhenHidden=FALSE)

  output$off <- renderUI({
    if(!is.null(input$genofile)) {
      all_columns <- readLines(input$genofile$datapath, n = 1) %>%
        str_split(pattern = ',') %>%
        unlist() %>%
        str_sort(numeric = T)
      ind_id <- all_columns[!(all_columns %in% c("id", "CHROM", "POS", "A", "B", "C", "D", "AB", "CD"))]
      selectInput(inputId = "off",
                  label = "Choose Offspring(s):",
                  multiple = TRUE,
                  choices = ind_id,
                  width = "95%")
    }
    })
  ### Chromosome input droplist ----
  output$scexist <- reactive({
    if(!is.null(input$scfile)) {
      return(1)
    }
  })
  outputOptions(output, 'scexist', suspendWhenHidden=FALSE)
  
  output$chr <- renderUI({
    if(!is.null(input$scfile)) {
      all_chr <- read_csv(input$scfile$datapath, col_types = cols(CHR = col_character())) %>% 
        .$CHR %>% 
        str_sort(numeric = T)
      selectInput(inputId = "chr",
                  label = "Choose Chromosome(s):",
                  multiple = TRUE,
                  choices = all_chr,
                  width = "95%")
    }
  })
  
  ### Show the result of which chromosome droplist ----
  output$show_chr <- renderUI({
    if(!is.null(input$scfile)) {
      all_chr <- read_csv(input$scfile$datapath, col_types = cols(CHR = col_character())) %>% 
        .$CHR %>% 
        str_sort(numeric = T) %>% 
        unique()
      selectInput(inputId = "show_chr", 
                  label = "Show results of chromosome:", 
                  multiple = FALSE,
                  choices = names(rec_plot_event()),
                  width = "100%",
                  selectize = FALSE)
    }
  })
  
  
  ## 2.1 Read-in and Settings ----
  ### Read-in data file ----
  dd <- reactive({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Loading data...")

    if(!is.null(input$genofile)) {
      read_csv(input$genofile$datapath, col_types = cols())
    } else {
      return(FALSE)
    }
  })

  ### Read-in scaffold order and orientation file ----
  sc <- reactive({
    if(!is.null(input$scfile)) {
      read_csv(input$scfile$datapath,
               col_types = cols(
                 CHR = col_character(),
                 order = col_double(),
                 scaffold = col_character(),
                 size = col_double(),
                 orientation = col_character()))
    } else {
      return(FALSE)
    }
  })
  ### Obtain chromosome setting ----
  ch <- reactive({
    req(input$chr)
  })
  ### Obtain offspring setting ----
  of <- reactive({
    req(input$off)
  })
  ### Obtain resolution setting ----
  rsn <- reactive({
    input$resolution
  })
  ### Locate recombination or not ----
  loc <- reactive({
    input$locate
  })
  ### Obtain algorithm setting ----
  alg <- reactive({
    input$algorithm
  })
  ### Obtain setting for PD----
  rad <- reactive({
    if (input$radius != '') {
      input$radius %>% as.numeric()
    } else {
      550
    }
  })
  stp <- reactive({
    if (input$step != '') {
      input$step %>% as.numeric()
    } else {
      17
    }
  })
  fstp <- reactive({
    if (input$finer_step != '') {
      input$finer_step %>% as.numeric()
    } else {
      1
    }
  })
  pd_thrsd <- reactive({
    if (input$pd_threshold != '') {
      input$pd_threshold %>% as.numeric()
    } else {
      0.9
    }
  })
  ### Obtain threshold setting ----
  thrsd <- reactive({
    if (input$threshold != '') {
      input$threshold %>% as.numeric()
    } else {
      50
    }
  })
  ### Obtain saving setting ----
  sop <- reactive({
      input$save_opt
  })
  sop2 <- reactive({
      input$save_opt2
  })
  ### Chromosome visualisation setting ----
  show_ch <- reactive({
    req(input$show_chr)
  })


  ## 2.2 Main ----
  ### 2.2.1 Function to analyse recombination locations ----
  rec_analyse <- function(data, sc_order, chromosome, offspring, loc, alg, thrsd, radius, step, finer_step, finer_threshold) {
    dt_path <- system.file("dictionary", package = "RecView")
    dictionary <- read_tsv(file = list.files(path = dt_path, pattern = "dict_complete.tsv", full.names = T), col_types = cols()) %>% 
      filter(!(Paternal == "N" & Maternal == "N"))
    if ((chromosome != "Z") && (chromosome != "X")) {
      dictionary <- dictionary %>% filter(Chromosome_type != "Z", Chromosome_type != "X")
    } else if (chromosome == "Z") {
      dictionary <- dictionary %>% filter(Chromosome_type == "Z")
    } else if (chromosome == "X") {
      dictionary <- dictionary %>% filter(Chromosome_type == "X")
    }

    AB_line <- c("Yes", "Yes", "Yes")
    CD_line <- c("Yes", "Yes", "Yes")
    case <- 0
    if (!("A" %in% colnames(data))) AB_line[1] <- "No"
    if (!("B" %in% colnames(data))) AB_line[2] <- "No"
    if (!("AB" %in% colnames(data))) AB_line[3] <- "No"
    if (!("C" %in% colnames(data))) CD_line[1] <- "No"
    if (!("D" %in% colnames(data))) CD_line[2] <- "No"
    if (!("CD" %in% colnames(data))) CD_line[3] <- "No"

    if ("No" %in% AB_line) case <- 1
    if ("No" %in% CD_line) case <- 2
    
    goo_inference <- function(tb) {
      goo <- dictionary %>% 
        filter(GT_string == tb$GT_string[1]) %>%
        select(Paternal, Maternal)
      return(goo)
    }
    
    dd_ind <- data %>%
      select(id, CHROM, POS,
             A, B, C, D,
             AB, CD,
             !!offspring) %>% 
      nest(pos_info = c("id", "CHROM", "POS")) %>% 
      mutate(tmp_id = seq(1, nrow(.)), GT_string = str_c(.$A, .$B, .$C, .$D, .$AB, .$CD, get(offspring, .), sep = "_")) %>% 
      nest(input = !tmp_id) %>% 
      mutate(output = map(input, goo_inference)) %>% 
      select(-tmp_id) %>% 
      unnest(cols = c("input", "output")) %>%
      filter(!is.na(Paternal), !is.na(Maternal)) %>% 
      unnest(cols = "pos_info")
    
    # dd_ind <- data %>%
    #   select(id, CHROM, POS,
    #          A, B, C, D,
    #          AB, CD,
    #          !!offspring) %>%
    #   unite(col = "GT_string", A:!!offspring, sep = '_') %>%
    #   left_join(dictionary, by = "GT_string") %>%
    #   filter(!is.na(Paternal), !is.na(Maternal))

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

    dd_ind_chr <- semi_join(dd_ind, chr, by = c("CHROM" = "scaffold")) %>%
      left_join(chr %>% select(scaffold, size, orientation, accumulated_size), by = c("CHROM" = "scaffold")) %>%
      mutate(scaffold_orientation = paste0(CHROM, ' ', orientation))
    dd_ind_chr <- dd_ind_chr %>% arrange(factor(CHROM, levels = chr$scaffold))

    dd_ind_chr_minus <- dd_ind_chr %>%
      filter(orientation == "-") %>%
      mutate(POS_chr = accumulated_size + size - POS) %>%
      select(POS_chr, everything())

    dd_ind_chr_plus <- dd_ind_chr %>%
      filter(orientation == "+") %>%
      mutate(POS_chr = accumulated_size + POS) %>%
      select(POS_chr, everything())

    dd_ind_chr <- bind_rows(dd_ind_chr_minus, dd_ind_chr_plus)

    res_raw <- dd_ind_chr %>%
      select(POS_chr, CHROM, scaffold_orientation, Paternal, Maternal) %>%
      arrange(POS_chr)

    res_par <- res_raw %>%
      select(POS_chr, CHROM, scaffold_orientation, Paternal) %>%
      `colnames<-`(., c("POS_chr", "CHROM", "scaffold_orientation", "Grandparent_of_origin")) %>%
      filter(Grandparent_of_origin != "N") %>%
      mutate(Grandparent_of_origin_value = NaN, Side = "Paternal")

    res_mar <- res_raw %>%
      select(POS_chr, CHROM, scaffold_orientation, Maternal) %>%
      `colnames<-`(., c("POS_chr", "CHROM", "scaffold_orientation", "Grandparent_of_origin")) %>%
      filter(Grandparent_of_origin != "N") %>%
      mutate(Grandparent_of_origin_value = NaN, Side = "Maternal")

    res <- list()
    res[[1]] <- bind_rows(res_par, res_mar) %>%
      mutate(Grandparent_of_origin_value = if_else(Grandparent_of_origin == "D", 0, Grandparent_of_origin_value),
             Grandparent_of_origin_value = if_else(Grandparent_of_origin == "C", 1, Grandparent_of_origin_value),
             Grandparent_of_origin_value = if_else(Grandparent_of_origin == "B", 2, Grandparent_of_origin_value),
             Grandparent_of_origin_value = if_else(Grandparent_of_origin == "A", 3, Grandparent_of_origin_value),
             POS_chr_mb = POS_chr/1e6)

    res[[1]]$Side <- factor(res[[1]]$Side, levels = c("Paternal", "Maternal"))
    res[[1]]$scaffold_orientation <- factor(res[[1]]$scaffold_orientation, levels = unique(res[[1]]$scaffold_orientation))

    #### Automatic detection of recombination location ----
    if (loc == "Yes") {
      if (alg == "PD") {
        seqle_mod <- function(x) {
          n <- length(x)  
          y <- x[-1L] != x[-n] + 1
          i <- c(which(y|is.na(y)),n) 
          return(x[head(c(0L,i)+1L,-1L)])
        }
        
        find_local_maxima <- function(tb) {
          tb[which(tb$diff == max(tb$diff)),]
        }

        goo_proportion <- function(input, tb, side, symbol) {
          tb <- tb[input$seg_start[1]:input$seg_end[1],]
          n <- tb %>% filter(get(side, tb) == !!symbol) %>% nrow()
          frequency <- n / nrow(tb)
          return(frequency)
        }
        
        finer_running_difference <- function(input, goo_proportion_diff_tb, tb, side, symbol, radius, finer_step, finer_threshold) {
          i <- input$row_start[1]
          while (goo_proportion_diff_tb$diff[i] >= finer_threshold) {
            i <- i + 1
          }
          cut_point <- seq(goo_proportion_diff_tb$cut_point[ifelse((input$row_start[1] - 1) == 0, 1, input$row_start[1] - 1)] + finer_step, goo_proportion_diff_tb$cut_point[i] - finer_step, finer_step)
          seg_1_proportion <- tibble(cut_point = cut_point, seg_start = cut_point - radius, seg_end = cut_point - 1) %>% 
            nest(input = !cut_point) %>% 
            mutate(p_1 = map(input, goo_proportion, tb, side, symbol)) %>% 
            unnest(cols = c("input", "p_1"))
          seg_2_proportion <- tibble(cut_point = cut_point, seg_start = cut_point, seg_end = cut_point + radius - 1) %>% 
            nest(input = !cut_point) %>% 
            mutate(p_2 = map(input, goo_proportion, tb, side, symbol)) %>% 
            unnest(cols = c("input", "p_2"))
          out <- inner_join(seg_1_proportion, seg_2_proportion, by = "cut_point") %>% 
            mutate(diff = abs(p_1 - p_2))
          return(out)
        }
        
        running_difference <- function(tb, radius, step, finer_step = 1, finer_threshold = 0.95, side, symbol) {
          snp <- nrow(tb)
          cut_point <- seq(1, snp - radius * 2, step) + radius
          seg_1_proportion <- tibble(cut_point = cut_point, seg_start = cut_point - radius, seg_end = cut_point - 1) %>% 
            nest(input = !cut_point) %>% 
            mutate(p_1 = map(input, goo_proportion, tb, side, symbol)) %>% 
            unnest(cols = c("input", "p_1"))
          seg_2_proportion <- tibble(cut_point = cut_point, seg_start = cut_point, seg_end = cut_point + radius - 1) %>% 
            nest(input = !cut_point) %>% 
            mutate(p_2 = map(input, goo_proportion, tb, side, symbol)) %>% 
            unnest(cols = c("input", "p_2"))
          segment_goo_proportion_diff <- inner_join(seg_1_proportion, seg_2_proportion, by = "cut_point") %>% 
            mutate(diff = abs(p_1 - p_2))
          
          rows_above_thrshd <- seqle_mod(which(segment_goo_proportion_diff$diff >= finer_threshold)) %>% 
            discard(is.na)
          if (length(rows_above_thrshd) > 0) {
            segment_goo_proportion_finer <- tibble(row_start = rows_above_thrshd) %>% 
              mutate(id = seq(1, nrow(.))) %>% 
              nest(input = !id) %>% 
              mutate(output = map(input, finer_running_difference, segment_goo_proportion_diff, tb, side, symbol, radius, finer_step, finer_threshold)) %>% 
              mutate(local_maxima = map(output, find_local_maxima)) %>% 
              select(output, local_maxima)
          } else {
            segment_goo_proportion_finer <- tibble(output = tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric()),
                                                  local_maxima = tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric()))
          }
          
          out <- list(segment_goo_proportion_diff, segment_goo_proportion_finer)
          return(out)
        }

        dd_pat <- res_raw %>% filter(Paternal != "N")
        dd_mat <- res_raw %>% filter(Maternal != "N")

        if (nrow(dd_pat) > 0) {
          dd_pat_mod <- running_difference(tb = dd_pat, radius = radius, step = step, finer_step = finer_step, finer_threshold = finer_threshold, side = "Paternal", symbol = "A")
        } else {
          dd_pat_mod <- list()
          dd_pat_mod[[1]] <- tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric())
          dd_pat_mod[[2]] <- tibble(output = tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric()))
        }
        if (nrow(dd_mat) > 0) {
          dd_mat_mod <- running_difference(tb = dd_mat, radius = radius, step = step, finer_step = finer_step, finer_threshold = finer_threshold, side = "Maternal", symbol = "C")
        } else {
          dd_mat_mod <- list()
          dd_mat_mod[[1]] <- tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric())
          dd_mat_mod[[2]] <- tibble(output = tibble(cut_point = as.numeric(), seg_start.x = as.numeric(), seg_end.x= as.numeric(), p_1 = as.numeric(), seg_start.y = as.numeric(), seg_end.y = as.numeric(), p_2 = as.numeric(), diff = as.numeric()))
        }
        
        dd_pat_mod_ <- bind_rows(dd_pat_mod[[1]], dd_pat_mod[[2]] %>% select(output) %>% unnest(cols = "output")) %>% 
          mutate(Side = "Paternal", POS_chr = dd_pat$POS_chr[.$cut_point])
        dd_mat_mod_ <- bind_rows(dd_mat_mod[[1]], dd_mat_mod[[2]] %>% select(output) %>% unnest(cols = "output")) %>% 
          mutate(Side = "Maternal", POS_chr = dd_mat$POS_chr[.$cut_point])

        res[[2]] <- bind_rows(dd_pat_mod_, dd_mat_mod_) %>% 
          distinct(Side, POS_chr, diff) %>% 
          mutate(POS_chr_mb = POS_chr/1e6)
        res[[2]]$Side <- factor(res[[2]]$Side, levels = c("Paternal", "Maternal"))
        
        dd_pat_mod_local_maxima <- dd_pat_mod[[2]] %>% 
          select(local_maxima) %>% 
          unnest(cols = "local_maxima") %>% 
          mutate(Side = "Paternal", POS_chr = dd_pat$POS_chr[.$cut_point])
        dd_mat_mod_local_maxima <- dd_mat_mod[[2]] %>% 
          select(local_maxima) %>% 
          unnest(cols = "local_maxima") %>% 
          mutate(Side = "Maternal", POS_chr = dd_mat$POS_chr[.$cut_point])
        
        res[[3]] <- rbind(dd_pat_mod_local_maxima, dd_mat_mod_local_maxima) %>%
          distinct(Side, POS_chr, diff) %>% 
          mutate(POS_chr_mb = POS_chr/1e6)
        res[[3]]$Side <- factor(res[[3]]$Side, levels = c("Paternal", "Maternal"))
        
        if (nrow(res[[3]]) == 0) {
          res[[3]] <- tibble(Side = c("Paternal", "Maternal"),
                             POS_chr = c("-", "-"),
                             POS_chr_mb = c("-", "-"),
                             diff = c("-", "-"))
        }

        res1_mod <- res[[1]] %>% select(POS_chr, scaffold_orientation)
        res1_mod$POS_chr <- as.character(res1_mod$POS_chr)
        res1_mod$scaffold_orientation <- as.character(res1_mod$scaffold_orientation)
        res[[3]]$POS_chr <- as.character(res[[3]]$POS_chr)
        res[[3]] <- left_join(res[[3]], res1_mod, by = "POS_chr") %>%
          select(Side, scaffold_orientation, everything()) %>%
          replace_na(list(scaffold_orientation = "-")) %>%
          `colnames<-`(.,c("Origin", "Scaffold & Orientation", "Chromosomal position (bp)", "Chromosomal position (Mb)", "Abs. GoO prop. diff."))

      } else if (alg == "CCS") {
        check_continuity <- function(tb, side) {
          tb <- tb %>%
            filter(Side == !!side) %>%
            mutate(goo_minus = c(0, .$Grandparent_of_origin[1:(length(.$Grandparent_of_origin)-1)]))

          similarity <- ifelse(tb$Grandparent_of_origin == tb$goo_minus, 1, 0)

          zeros <- which(similarity == 0)
          sepa <- list()
          if (zeros[length(zeros)] != nrow(tb)) {
            sepa[[1]] <- c(zeros, nrow(tb))
            sepa[[2]] <- "incl"
          } else {
            sepa[[1]] <- zeros
            sepa[[2]] <- "excl"
          }

          if (sepa[[2]] == "incl") {
            tb_zero <- tb[sepa[[1]][1:(length(sepa[[1]])-1)],] %>% mutate(continuity = 0, abs_cumulative_continuity = 0)
          } else {
            tb_zero <- tb[sepa[[1]][1:length(sepa[[1]])],] %>% mutate(continuity = 0, abs_cumulative_continuity = 0)
          }

          tb_list <- list()
          i <- 1
          for (j in 2:length(sepa[[1]])) {
            if (sepa[[1]][j] - 1 != sepa[[1]][j-1])  {
              tb_list[[i]] <- tb[(sepa[[1]][j-1] + 1):(sepa[[1]][j]-1),]
              if (sepa[[1]][j] == nrow(tb)) {
                if (sepa[[2]] == "incl") {
                  tb_list[[i]] <- tb[(sepa[[1]][j-1] + 1):sepa[[1]][j],]
                } else {
                  tb_list[[i]] <- tb[(sepa[[1]][j-1] + 1):(sepa[[1]][j]-1),]
                }
              }
              if (unique(tb_list[[i]]$Grandparent_of_origin) %in% c("A","C")) {
                tb_list[[i]] <- tb_list[[i]] %>% mutate(continuity = seq(1,nrow(.)), abs_cumulative_continuity = rep(nrow(.), times = nrow(.)))
              } else if (unique(tb_list[[i]]$Grandparent_of_origin) %in% c("B","D")) {
                tb_list[[i]] <- tb_list[[i]] %>% mutate(continuity = seq(-1,-nrow(.),-1), abs_cumulative_continuity = rep(nrow(.), times = nrow(.)))
              }
              i <- i + 1
            }
          }

          tb <- reduce(tb_list, bind_rows) %>%
            bind_rows(tb_zero) %>%
            arrange(POS_chr) %>%
            select(-goo_minus)
          return(tb)
        }

        backtrack_position <- function(tb, threshold) {
          tb <- tb %>%
            mutate(index = seq(1:nrow(.))) %>%
            select(POS_chr, POS_chr_mb, continuity, abs_cumulative_continuity, index)

          tb_1 <- tb %>%
            filter(tb$continuity == threshold | tb$continuity == -threshold) %>%
            arrange(POS_chr)

          tb_2 <- list()
          if (nrow(tb_1) > 1) {
            k <- 1
            for (j in 2:nrow(tb_1)) {
              if (tb_1$continuity[j] == -tb_1$continuity[j - 1]) {
                tb_2[[k]] <- tb_1[(j-1):j,]
                k <- k + 1
              }
            }
          }

          rows_retrieve <- function(x, threshold, tb) {
            x <- x %>% arrange(POS_chr)
            t_index_1 <- x$abs_cumulative_continuity[1] - threshold + x$index[1]
            t_index_2 <- x$index[2] - threshold

            t_row_1 <- tb %>% filter(index == t_index_1)
            t_row_2 <- tb %>% filter(index == t_index_2)

            out <- tibble(POS_chr = round(mean(c(t_row_1$POS_chr, t_row_2$POS_chr))), Left_boundary = t_row_1$POS_chr, Right_boundary = t_row_2$POS_chr)
            return(out)
          }

          if (length(tb_2) > 0) {
            out <- lapply(tb_2, rows_retrieve, threshold = threshold, tb = tb) %>%
              reduce(bind_rows)
          } else {
            out <- tibble(POS_chr = as.numeric(), Left_boundary = as.numeric(), Right_boundary= as.numeric())
          }

          return(out)
        }

        tb_pat <- check_continuity(res[[1]], side = "Paternal")
        tb_mat <- check_continuity(res[[1]], side = "Maternal")
        res[[2]] <- bind_rows(tb_pat, tb_mat)

        tb_pat_mod <- backtrack_position(tb_pat, threshold = thrsd)
        if (nrow(tb_pat_mod) != 0) {
          tb_pat_mod <- tb_pat_mod %>%
            mutate(Origin = "Paternal", POS_chr_mb = POS_chr / 1e6) %>%
            select(Origin, POS_chr, POS_chr_mb, everything())
        } else {
          tb_pat_mod <- tibble(Origin = "Paternal",
                               POS_chr = "-",
                               POS_chr_mb = "-",
                               Left_boundary = "-",
                               Right_boundary = "-")
        }
        tb_pat_mod$POS_chr <- as.character(tb_pat_mod$POS_chr)
        tb_pat_mod$POS_chr_mb <- as.character(tb_pat_mod$POS_chr_mb)
        tb_pat_mod$Left_boundary <- as.character(tb_pat_mod$Left_boundary)
        tb_pat_mod$Right_boundary <- as.character(tb_pat_mod$Right_boundary)

        tb_mat_mod <- backtrack_position(tb_mat, threshold = thrsd)
        if (nrow(tb_mat_mod) != 0) {
          tb_mat_mod <- tb_mat_mod %>%
            mutate(Origin = "Maternal", POS_chr_mb = POS_chr / 1e6) %>%
            select(Origin, POS_chr, POS_chr_mb, everything())
        } else {
          tb_mat_mod <- tibble(Origin = "Maternal",
                               POS_chr = "-",
                               POS_chr_mb = "-",
                               Left_boundary = "-",
                               Right_boundary = "-")
        }
        tb_mat_mod$POS_chr <- as.character(tb_mat_mod$POS_chr)
        tb_mat_mod$POS_chr_mb <- as.character(tb_mat_mod$POS_chr_mb)
        tb_mat_mod$Left_boundary <- as.character(tb_mat_mod$Left_boundary)
        tb_mat_mod$Right_boundary <- as.character(tb_mat_mod$Right_boundary)

        res[[3]] <- bind_rows(tb_pat_mod, tb_mat_mod) %>%
          unite(col = "id", c("Origin", "Right_boundary"), remove = FALSE)
        res1_mod <- res[[1]] %>% select(POS_chr, Side, scaffold_orientation) %>%
          unite(col = "id", c("Side", "POS_chr"), remove = TRUE)
        res1_mod$scaffold_orientation <- as.character(res1_mod$scaffold_orientation)
        res[[3]] <- left_join(res[[3]], res1_mod, by = "id") %>%
          select(-id) %>%
          select(Origin, scaffold_orientation, everything()) %>%
          replace_na(list(scaffold_orientation = "-")) %>%
          `colnames<-`(.,c("Origin", "Scaffold & Orientation", "Chromosomal position (bp)", "Chromosomal position (Mb)", "Left boundary (bp)", "Right boundary (bp)"))

      }
      res[[3]]$Origin <- factor(res[[3]]$Origin, levels = c("Paternal", "Maternal"))
    }

    #### Informative SNP density ----
    interval <- 1e5
    digits <- (log10(interval) %>% floor()) + 1
    if (case != 1) {
      breaks_pat <- seq(0, ceiling(max(get("POS_chr", res[[1]] %>% filter(Side == "Paternal")))/10^digits) * 10^digits, interval)
      hist_pat <- hist(get("POS_chr", res[[1]] %>% filter(Side == "Paternal")), breaks = breaks_pat, plot = F)
      max_break <- max(hist_pat$breaks)

      snp_density_pat <- tibble(POS_chr = hist_pat$mids[1:length(hist_pat$mids)],
                                breaks = hist_pat$breaks[1:(length(hist_pat$breaks)-1)],
                                Counts = hist_pat$counts[1:length(hist_pat$counts)]) %>%
        mutate(Density = Counts / interval,
               Precision = interval / Counts,
               POS_chr_mb = POS_chr / 1e6,
               Side = "Paternal")
    } else {
      snp_density_pat <- tibble(POS_chr = NA, breaks = NA, Counts = NA,
                                Density = NA, Precision = NA, POS_chr_mb = NA,
                                Side = "Paternal")
    }

    if (case != 2) {
      breaks_mat <- seq(0, ceiling(max(get("POS_chr", res[[1]] %>% filter(Side == "Maternal")))/10^digits) * 10^digits, interval)
      hist_mat <- hist(get("POS_chr", res[[1]] %>% filter(Side == "Maternal")), breaks = breaks_mat, plot = F)
      max_break <- max(hist_mat$breaks)

      snp_density_mat <- tibble(POS_chr = hist_mat$mids[1:length(hist_mat$mids)],
                                breaks = hist_mat$breaks[1:(length(hist_mat$breaks)-1)],
                                Counts = hist_mat$counts[1:length(hist_mat$counts)]) %>%
        mutate(Density = Counts / interval,
               Precision = interval / Counts,
               POS_chr_mb = POS_chr / 1e6,
               Side = "Maternal")
    } else {
      snp_density_mat <- tibble(POS_chr = NA, breaks = NA, Counts = NA,
                                Density = NA, Precision = NA, POS_chr_mb = NA,
                                Side = "Maternal")
    }

    res[[4]] <- rbind(snp_density_pat, snp_density_mat) %>%
      mutate(Precision = na_if(Precision, Inf))
    res[[4]]$Side <- factor(res[[4]]$Side, levels = c("Paternal", "Maternal"))


    if (loc == "Yes") {
       left <- get("breaks", res[[4]]) %>% unique() %>% .[1:length(.)]
       right <- c(get("breaks", res[[4]]) %>% unique() %>% .[2:length(.)],
                  max_break)

       intervals_retrieve <- function(x, tb) {
         side <- x$Origin
         if (x$`Chromosomal position (bp)` != '-') {
          left_boolean <- as.numeric(x$`Chromosomal position (bp)`) >= left
          right_boolean <- as.numeric(x$`Chromosomal position (bp)`) < right
          boolean <- which((left_boolean & right_boolean) == TRUE)

          ee <- tb %>%
            filter(Side == side) %>%
            .[boolean,] %>%
            get("Precision",.) %>%
            round()

         } else {
           ee = "-"
         }
         out <- x %>%
           mutate(Esimated_error = ee)
         out$Esimated_error <- as.character(out$Esimated_error)
         return(out)
       }

       res[[3]] <- res[[3]] %>%
         mutate(index = seq(1:nrow(.))) %>%
         nest(input = !index) %>%
         mutate(output = map(input, intervals_retrieve, tb = res[[4]])) %>%
         select(-input) %>%
         unnest(cols = "output") %>%
         select(-index) %>%
         select(Origin, `Scaffold & Orientation`, `Chromosomal position (bp)`, Esimated_error, `Chromosomal position (Mb)`, everything()) %>%
         `colnames<-`(.,c("Origin", "Scaffold & Orientation", "Chromosomal position (bp)",
                          "Precision (bp)", "Chromosomal position (Mb)", colnames(res[[3]])[5:length(colnames(res[[3]]))]))

    }

    return(res)
  }

  ### 2.2.2 Function to plot results ----
  rec_plot_0 <- function(tb) {
    res <- tb$res
    chromosome <- tb$chromosome %>% unique()
    offspring <- tb$offspring %>% unique()
    resolution <- tb$resolution %>% unique()
    location <- tb$location %>% unique()
    algorithm <- tb$algorithm %>% unique()
    threshold <- tb$threshold %>% unique()

    p0 <- res[[1]] %>%
      mutate(Chromosome = chromosome,
             Offspring = offspring) %>%
      select(Offspring, Chromosome, Side, POS_chr, scaffold_orientation, Grandparent_of_origin) %>%
      `colnames<-`(.,c("Offspring", "Chromosome", "Origin", "Position", "Scaffold_Orientation", "Grandparent-of-origin"))

    return(p0)
  }


  rec_plot_1 <- function(tb) {
    res <- tb$res
    chromosome <- tb$chromosome %>% unique()
    offspring <- tb$offspring %>% unique()
    resolution <- tb$resolution %>% unique()
    location <- tb$location %>% unique()
    algorithm <- tb$algorithm %>% unique()
    threshold <- tb$threshold %>% unique()

    if (max(res[[1]]$POS_chr_mb) > 15) {
      break_step <- 10
    } else if (max(res[[1]]$POS_chr_mb) > 5) {
      break_step <- 2.5
    } else if (max(res[[1]]$POS_chr_mb) > 2) {
      break_step <- 1
    } else if (max(res[[1]]$POS_chr_mb) > 1) {
      break_step <- 0.25
    }

    if (nrow(res[[1]]) < 2000) {
      point_size <- 6
    } else {
      point_size <- exp(-8.407e-08 * nrow(res[[1]]) + 0.9498) %>% round(1)
    }
    
    macaron_mod <- c("#940214", "#3e1f16", "#123358", "#e8c001", "#015c7a", "#85ab87", "#d75624", "#dd395a")
    n_scaffold <- res[[1]]$scaffold_orientation %>% unique() %>% length()
    n_colour <- ceiling(n_scaffold / 2)
    colour_scheme <- vector()
    if (n_colour <= length(macaron_mod)) {
      macaron_mod_light <- colorspace::lighten(macaron_mod[1:n_colour], amount = 0.5) %>% rev()
    } else {
      macaron_mod <- colorRampPalette(macaron_mod)(n_colour)
      macaron_mod_light <- colorspace::lighten(macaron_mod, amount = 0.5) %>% rev()
    }
    for (i in 1:length(macaron_mod_light)) {
      colour_scheme <- c(colour_scheme, macaron_mod[i], macaron_mod_light[i])
    }

    ex_tb <- tibble(POS_chr_mb = 0, Grandparent_of_origin_value = -0.13, Side = "Maternal")
    ex_tb$Side <- factor(ex_tb$Side, levels = c("Paternal", "Maternal"))
    if (resolution == "Low") {
      p1 <- ggplot(data = res[[1]]) +
        scattermore::geom_scattermore(aes(x = POS_chr_mb, y = Grandparent_of_origin_value, color = scaffold_orientation),
                                      pointsize = point_size, position = position_jitter(height = 0.2), pixels = c(1500, 264)) +
        scattermore::geom_scattermore(data = ex_tb, aes(x = POS_chr_mb, y = Grandparent_of_origin_value), alpha = 0,
                                      pointsize = 0, pixels = c(1500, 264))
    } else {
      p1 <- ggplot(data = res[[1]]) +
        geom_point(aes(x = POS_chr_mb, y = Grandparent_of_origin_value, color = scaffold_orientation),
                   size = 0.1, position = position_jitter(height = 0.2)) +
        geom_point(data = ex_tb, aes(x = POS_chr_mb, y = Grandparent_of_origin_value), alpha = 0,
                                      size = 0)
    }
    p1 <- p1 + scale_x_continuous(breaks = seq(0, max(res[[1]]$POS_chr_mb), break_step),
                                  limits = c(-max(res[[1]]$POS_chr_mb)*0.01, max(res[[1]]$POS_chr_mb)*1.01), 
                                  expand = c(0,0)) +
      labs(x = "Chromosomal position (Mb)",
           y = "Grandparent-of-origin",
           title = paste0(offspring, '; Chr', chromosome),
           color = "Scaffold\nOrientation") +
      facet_wrap(~Side, scales = "free_y", nrow = 2, strip.position = "right") +
      scale_y_continuous(breaks = c(0,1,2,3), labels = c("MGM", "MGF", "PGM", "PGF")) +
      scale_color_manual(values = colour_scheme) +
      guides(color = guide_legend(override.aes = list(size=10))) +
      theme(plot.title = element_text(size = 20, face = "bold"),
            legend.position = "right",
            legend.title = element_text(size = 12, face = "bold", color = "black"),
            legend.text = element_text(size = 10, color = "black"),
            legend.margin = margin(0,0,0,0),
            plot.margin = margin(0.3,0,0.3,0),
            panel.background = element_blank(),
            panel.border = element_rect(color = "grey10", linewidth = 1, fill = NA),
            panel.grid.major.x = element_line(color = "grey60", size = 0.2),
            panel.grid.major.y = element_blank(),
            panel.spacing.y = unit(2.2, "lines"),
            axis.text.x = element_text(size = 16, color = "black"),
            axis.text.y = element_text(size = 20, color = "black"),
            axis.title = element_text(size = 20, face = "bold", color = "black"),
            strip.text = element_text(size = 16, color = "black"),
            strip.background = element_rect(fill = "grey90", color = "grey10", linewidth = 1))

    return(p1)
  }

  rec_plot_2 <- function(tb) {
    res <- tb$res
    chromosome <- tb$chromosome %>% unique()
    offspring <- tb$offspring %>% unique()
    resolution <- tb$resolution %>% unique()
    location <- tb$location %>% unique()
    algorithm <- tb$algorithm %>% unique()
    threshold <- tb$threshold %>% unique()

    if (max(res[[1]]$POS_chr_mb) > 15) {
      break_step <- 10
    } else if (max(res[[1]]$POS_chr_mb) > 5) {
      break_step <- 2.5
    } else if (max(res[[1]]$POS_chr_mb) > 2) {
      break_step <- 1
    } else if (max(res[[1]]$POS_chr_mb) > 1) {
      break_step <- 0.25
    }

    round_any <- function(x, accuracy, f=round) {f(x/ accuracy) * accuracy}
    
    if (algorithm == "PD") {
      p2 <- ggplot(data = res[[2]]) +
        geom_line(aes(x = POS_chr_mb, y = diff, group = Side), size = 0.2, color = "#C91959") +
        geom_point(aes(x = POS_chr_mb, y = diff), size = 0.5, color = "#C91959") +
        scale_x_continuous(breaks = seq(0, max(res[[1]]$POS_chr_mb), break_step)) +
        scale_y_continuous(limits = c(-0.05,1.05), breaks = seq(0,1,0.2)) +
        coord_cartesian(xlim = c(-max(res[[1]]$POS_chr_mb)*0.01, max(res[[1]]$POS_chr_mb)*1.01), expand = F) +
        facet_wrap(~Side, nrow = 2, strip.position = "right") +
        labs(x = "Chromosomal position (Mb)",
             y = "Abs. GoO prop. diff.",
             title = paste0(offspring, '; Chr', chromosome)) +
        theme(plot.title = element_text(size = 20, face = "bold"),
              plot.margin = margin(0.3,0,0.3,0, unit = "mm"),
              panel.background = element_blank(),
              panel.border = element_rect(color = "grey10", size = 1, fill = NA),
              panel.grid.major = element_line(color = "grey60", size = 0.2),
              panel.spacing.y = unit(2.2, "lines"),
              axis.text = element_text(size = 16, color = "black"),
              axis.title = element_text(size = 20, face = "bold", color = "black"),
              strip.text = element_text(size = 16, color = "black"),
              strip.background = element_rect(fill = "grey90", color = "grey10", size = 1))
      
      p2 <- list(p2)
    } else if (algorithm == "CCS") {
      tb <- res[[2]] %>% mutate(sign = ifelse(continuity >= 0, 1, -1))
      
      tb_list <- list()
      zeros <- which(tb$continuity == 0)
      k <- 1
      for (i in 1:(length(zeros) - 1)) {
        intermediate <- tb[c(zeros[i], zeros[i + 1] - 1),] %>% mutate(line_group = k)
        
        if(tb$abs_cumulative_continuity[zeros[i] + 1] >= threshold) {
          intermediate <- intermediate %>% mutate(sign = prod(sign))
        } else {
          intermediate <- intermediate %>% mutate(sign = 2)
        }
        
        tb_list[[k]] <- intermediate %>% select(Side, POS_chr_mb, sign, line_group)
        k <- k + 1
      }
      
      if(zeros[length(zeros)] != nrow(tb)) {
        intermediate <- tb[c(zeros[length(zeros)], nrow(tb)),] %>% mutate(line_group = length(tb_list) + 1)
        
        if(tb$abs_cumulative_continuity[zeros[length(zeros)] + 1] >= threshold) {
          intermediate <- intermediate %>% mutate(sign = prod(sign))
        } else {
          intermediate <- intermediate %>% mutate(sign = 2)
        }
        tb_list[[length(tb_list) + 1]] <- intermediate %>% select(Side, POS_chr_mb, sign, line_group)
      }
      
      tb_list <- reduce(tb_list, bind_rows) %>% mutate(CCS = "NA")
      tb_list$CCS[tb_list$sign == -1] <- paste0("<= -", threshold)
      tb_list$CCS[tb_list$sign == 1] <- paste0(">= ", threshold)
      tb_list$CCS[tb_list$sign == 2] <- paste0("< |", threshold, "|")
      
      tb_list$Side <- factor(tb_list$Side, levels = c("Paternal", "Maternal"))
      tb_list$CCS <- factor(tb_list$CCS, levels = c(paste0("<= -", threshold), paste0(">= ", threshold), paste0("< |", threshold, "|")))
      cols <- setNames(c("#81BECE", "#E45545", "black"), c(paste0("<= -", threshold), paste0(">= ", threshold), paste0("< |", threshold, "|")))
      
      p2 <- ggplot(data = tb_list) +
        geom_line(aes(x = POS_chr_mb, y = "CCS", group = line_group, color = CCS), linewidth = 13) +
        facet_wrap(~Side, scales = "free_y", nrow = 2, strip.position = "right") +
        scale_x_continuous(breaks = seq(0, max(res[[1]]$POS_chr_mb), break_step)) +
        coord_cartesian(xlim = c(-max(res[[1]]$POS_chr_mb)*0.01, max(res[[1]]$POS_chr_mb)*1.01), expand = F) +
        scale_color_manual(values = cols) +
        labs(x = "Chromosomal position (Mb)", 
             title = paste0(offspring, '; Chr', chromosome)) +
        guides(color = guide_legend(override.aes = list(size=8))) +
        theme(plot.title = element_text(size = 20, face = "bold"),
              legend.position = "right",
              legend.title = element_text(size = 12, face = "bold", color = "black"),
              legend.text = element_text(size = 10, color = "black"),
              legend.margin = margin(0,0,0,0),
              plot.margin = margin(0.3,0,0.3,0),
              panel.background = element_blank(),
              panel.border = element_rect(color = "grey10", size = 1, fill = NA),
              panel.grid.major.x = element_line(color = "grey60", size = 0.2),
              panel.grid.major.y = element_blank(),
              panel.spacing.y = unit(2.2, "lines"),
              axis.text.x = element_text(size = 16, color = "black"),
              axis.text.y = element_blank(),
              axis.title.x = element_text(size = 20, face = "bold", color = "black"),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.text = element_text(size = 16, color = "black"),
              strip.background = element_rect(fill = "grey90", color = "grey10", size = 1))

      p2 <- list(p2)
    }
    return(p2)
  }

  rec_plot_3 <- function(tb) {
    res <- tb$res
    chromosome <- tb$chromosome %>% unique()
    offspring <- tb$offspring %>% unique()
    resolution <- tb$resolution %>% unique()
    location <- tb$location %>% unique()
    algorithm <- tb$algorithm %>% unique()
    threshold <- tb$threshold %>% unique()

    p3 <- res[[3]] %>%
      select(-`Chromosomal position (Mb)`) %>%
      mutate(Offspring = offspring, Chromosome = chromosome) %>%
      select(Offspring, Chromosome, Origin, `Scaffold & Orientation`, `Chromosomal position (bp)`, `Precision (bp)`)
    return(p3)
  }

  rec_plot_4 <- function(tb) {
    res <- tb$res
    chromosome <- tb$chromosome %>% unique()
    offspring <- tb$offspring %>% unique()
    resolution <- tb$resolution %>% unique()
    location <- tb$location %>% unique()
    algorithm <- tb$algorithm %>% unique()
    threshold <- tb$threshold %>% unique()

    if (max(res[[1]]$POS_chr_mb) > 15) {
      break_step <- 10
    } else if (max(res[[1]]$POS_chr_mb) > 5) {
      break_step <- 2.5
    } else if (max(res[[1]]$POS_chr_mb) > 2) {
      break_step <- 1
    } else if (max(res[[1]]$POS_chr_mb) > 1) {
      break_step <- 0.25
    }

    pat_max <- get("Density", res[[4]] %>% filter(Side == "Paternal")) %>% max()
    mat_max <- get("Density", res[[4]] %>% filter(Side == "Maternal")) %>% max()
    range_upper <- max(pat_max*1e5, mat_max*1e5)

    p4 <- ggplot(data = res[[4]]) +
      geom_line(aes(x = POS_chr_mb, y = Density*1e5, group = Side), size = 0.2) +
      scale_x_continuous(breaks = seq(0, max(res[[1]]$POS_chr_mb), break_step)) +
      coord_cartesian(ylim = c(-range_upper*0.1, range_upper*1.1), xlim = c(-max(res[[1]]$POS_chr_mb)*0.01, max(res[[1]]$POS_chr_mb)*1.01), expand = F) +
      labs(x = "Chromosomal position (Mb)",
           y = expression(bold(atop("Density of", "informative SNPs" (''%*%"10"^"-5")))),
           title = paste0(offspring, '; Chr', chromosome)) +
      facet_wrap(~Side, nrow = 2, strip.position = "right") +
      theme(plot.title = element_text(size = 20, face = "bold"),
            plot.margin = margin(5,10,0,0),
            panel.background = element_blank(),
            panel.border = element_rect(color = "grey10", size = 1, fill = NA),
            panel.grid.major = element_line(color = "grey60", size = 0.2),
            panel.spacing.y = unit(2.2, "lines"),
            axis.text = element_text(size = 16, color = "black"),
            axis.title = element_text(size = 20, face = "bold", color = "black"),
            strip.text = element_text(size = 16, color = "black"),
            strip.background = element_rect(fill = "grey90", color = "grey10", size = 1))

    # if (location == "Yes") {
    #   res3_plot <- res[[3]] %>% filter(`Chromosomal position (Mb)` != "-") %>%
    #     `colnames<-`(.,c("Side", colnames(.)[-1]))
    #   res3_plot$`Chromosomal position (Mb)` <- as.numeric(res3_plot$`Chromosomal position (Mb)`)
    #   if (nrow(res3_plot) > 0) {
    #     p4 <- p4 +
    #       geom_point(data =res3_plot,
    #                  aes(x = `Chromosomal position (Mb)`, y = 0), pch = 23, fill = "firebrick2", size = 2.3)
    #   }
    # }

    return(p4)
  }

  ### 2.2.3 Conduct analysis ----
  rec_plot_event <- eventReactive(input$click, {
    if (!isFALSE(dd()) && !isFALSE(sc()) && length(ch()) > 0) {

      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Running analysis...")

      offspring <- of()
      for (i in 1:length(offspring)) {
        offspring[i] <- str_trim(offspring[i], side = "both")
      }

      dd_in <- dd()
      sc_in <- sc()
      loc_in <- loc()
      alg_in <- alg()
      thrsd_in <- thrsd()
      rad_in <- rad()
      stp_in <- stp()
      fstp_in <- fstp()
      fthrsd_in <- pd_thrsd()
      if (is.null(rsn())) {
        rsn_in <- "High"
      } else {
        rsn_in <- rsn()
      }

      ch_in <- ch()
      for (i in 1:length(ch_in)) {
        ch_in[i] <- str_trim(ch_in[i], side = "both")
      }

      avail_cores <- parallel::detectCores()
      need_cores <- length(offspring) * length(ch_in)
      if (need_cores <= avail_cores - 1) {
        use_cores <- need_cores
      } else {
        num <- seq(1:avail_cores - 1)
        multiplication <- ceiling(need_cores / num) %>% min()
        multi_num <- multiplication * num
        abs_diff <- abs(need_cores - multi_num)
        index <- which(abs_diff == min(abs_diff))
        if (length(index) == 1) {
          use_cores <- num[index]
        } else {
          use_cores <- num[max(index)]
        }
      }
      msg_out <- paste0("\n ------------------------------\n Note: There are ", avail_cores, " cores available.\n The analysis is run on ", use_cores, " cores.\n ------------------------------\n\n")
      cat(msg_out)
      doParallel::registerDoParallel(use_cores)
      
      p_list <- foreach(i = 1:length(ch_in), .combine = bind_rows) %:%
        foreach(j = 1:length(offspring), .combine = bind_rows) %dopar% {
          res <- rec_analyse(data = dd_in, sc_order = sc_in, chromosome = ch_in[i], offspring = offspring[j],
                             loc = loc_in, alg = alg_in, thrsd = thrsd_in,
                             radius = rad_in, step = stp_in, finer_step = fstp_in, finer_threshold = fthrsd_in)
          
          out <- tibble(index = paste(ch_in[i], offspring[j], sep = "_"),
                        chromosome = ch_in[i],
                        offspring = offspring[j],
                        res = res,
                        resolution = rsn_in,
                        location = loc_in,
                        algorithm = alg_in,
                        threshold = thrsd_in) %>%
            nest(input = !index) %>%
            mutate(fig_1 = map(input, rec_plot_1),
                   fig_4 = map(input, rec_plot_4),
                   fig_0 = map(input, rec_plot_0))
          
          if (loc_in == "Yes") {
            out <- out %>%
              mutate(fig_2 = map(input, rec_plot_2),
                     fig_3 = map(input, rec_plot_3))
          }
          return(out)
        }
      p_list <- p_list %>% separate(col = "index", into = c("chromosome", "offspring"), sep = "_")
      
      p_list_list <- list()
      for (i in 1:length(ch_in)) {
        p_list_mod <- p_list %>% filter(chromosome == ch_in[i])
        if (loc_in == "Yes") {
          p_list_list[[i]] <- list(p_list_mod$fig_1, p_list_mod$fig_2, p_list_mod$fig_3, p_list_mod$fig_4, p_list_mod$fig_0)
        } else {
          p_list_list[[i]] <- list(p_list_mod$fig_1, p_list_mod$fig_4, p_list_mod$fig_0)
        }
      }
      names(p_list_list) <- ch_in
      doParallel::stopImplicitCluster()
      return(p_list_list)
    }
  })

  
  #<<< Output to User Interface >>>#
  observeEvent(c(show_ch(), input$click), {
    num_off <- length(of())
    hei <- num_off * 450

    rpe <- rec_plot_event()

    if (length(rpe) == 1) {
      stamp <- format(Sys.time(), "%m%d%H%M")
      inference_filename <- paste0('Chromosome_', ch(), '_GoO_TimeStamp', stamp, '.csv')
      if (loc() == "Yes") {
        plot_filename <- paste0('Chromosome_', ch(), '_', alg(), '_TimeStamp', stamp, '.pdf')
      } else {
        plot_filename <- paste0('Chromosome_', ch(), '_TimeStamp', stamp, '.pdf')
      }
      table_filename <- paste0('Chromosome_', ch(), '_', alg(), '_RecPos_TimeStamp', stamp, '.csv')

      if (loc() == "Yes") {
        hei_mod <- 15
        wid_mod <- 15

        if ("GoO Inferences" %in% sop2()) {
          write_csv(reduce(rpe[[1]][[5]], rbind), file = inference_filename)
        }

        out_temp <- list()
        if ("Plots" %in% sop2()) {
          for (i in 1:length(rpe[[1]][[1]])) {
            goo_legend_1 <- ggpubr::get_legend(rpe[[1]][[1]][[i]])
            goo_legend_2 <- ggpubr::get_legend(rpe[[1]][[2]][[i]])
            out_temp[[i]] <- ggpubr::ggarrange(plotlist =  list(rpe[[1]][[1]][[i]] + theme(plot.title = element_text(size = 36, face = "bold"), legend.position = "none"), 
                                                                NULL,
                                                                rpe[[1]][[2]][[i]][[1]] + theme(plot.title = element_blank(), legend.position = "none"),
                                                                NULL,
                                                                rpe[[1]][[4]][[i]] + theme(plot.title = element_blank())), ncol = 1, heights = c(1.1,0.1,1.1,0.1,1), align = "v")
            
            out_temp[[i]] <- ggpubr::ggarrange(out_temp[[i]], 
                                               ggpubr::ggarrange(goo_legend_1, goo_legend_2, NULL, ncol = 1), 
                                               ncol = 2, widths = c(10,1))
          }
          
          ggpubr::ggexport(plotlist = out_temp, filename = plot_filename, verbose = FALSE, width = wid_mod, height = hei_mod)
        }

        if ("Locations" %in% sop2()) {
          write_csv(reduce(rpe[[1]][[3]], rbind), file = table_filename)
        }

      } else {
        hei_mod <- 9
        wid_mod <- 15

        if ("GoO Inferences" %in% sop()) {
          write_csv(reduce(rpe[[1]][[3]], rbind), file = inference_filename)
        }

        out_temp <- list()
        if ("Plots" %in% sop()) {
          for (i in 1:length(rpe[[1]][[1]])) {
            out_temp[[i]] <- ggpubr::ggarrange(plotlist =  list(rpe[[1]][[1]][[i]], 
                                                                NULL,
                                                                rpe[[1]][[2]][[i]]+ theme(plot.title = element_blank())), ncol = 1, heights = c(1,0.08,1), common.legend = T, align = "v", legend = "right")
          }
          ggpubr::ggexport(plotlist = out_temp, filename = plot_filename, verbose = FALSE, width = wid_mod, height = hei_mod)
        }
      }
    } else {
      for (i in 1:length(rpe)) {
        stamp <- format(Sys.time(), "%m%d%H%M")
        inference_filename <- paste0('Chromosome_', names(rpe)[i], '_GoO_TimeStamp', stamp, '.csv')
        if (loc() == "Yes") {
          plot_filename <- paste0('Chromosome_', names(rpe)[i], '_', alg(), '_TimeStamp', stamp, '.pdf')
        } else {
          plot_filename <- paste0('Chromosome_', names(rpe)[i], '_TimeStamp', stamp, '.pdf')
        }
        table_filename <- paste0('Chromosome_', names(rpe)[i], '_', alg(), '_RecPos_TimeStamp', stamp, '.csv')

        if (loc() == "Yes") {
          hei_mod <- 15
          wid_mod <- 15

          if ("GoO Inferences" %in% sop2()) {
            write_csv(reduce(rpe[[i]][[5]], rbind), file = inference_filename)
          }

          out_temp <- list()
          if ("Plots" %in% sop2()) {
            for (j in 1:length(rpe[[i]][[1]])) {
              goo_legend_1 <- ggpubr::get_legend(rpe[[i]][[1]][[j]])
              goo_legend_2 <- ggpubr::get_legend(rpe[[i]][[2]][[j]])
              out_temp[[j]] <- ggpubr::ggarrange(plotlist =  list(rpe[[i]][[1]][[j]] + theme(plot.title = element_text(size = 36, face = "bold"), legend.position = "none"), 
                                                                  NULL,
                                                                  rpe[[i]][[2]][[j]][[1]] + theme(plot.title = element_blank(), legend.position = "none"), 
                                                                  NULL,
                                                                  rpe[[i]][[4]][[j]] + theme(plot.title = element_blank())), ncol = 1, heights = c(1.1,0.1,1.1,0.1,1), align = "v")
              
              out_temp[[j]] <- ggpubr::ggarrange(out_temp[[j]], 
                                                 ggpubr::ggarrange(goo_legend_1, goo_legend_2, NULL, ncol = 1), 
                                                 ncol = 2, widths = c(10,1))
            }
              
            ggpubr::ggexport(plotlist = out_temp, filename = plot_filename, verbose = FALSE, width = wid_mod, height = hei_mod)
          }

          if ("Locations" %in% sop2()) {
            write_csv(reduce(rpe[[i]][[3]], rbind), file = table_filename)
          }

        } else {
          hei_mod <- 9
          wid_mod <- 15

          if ("GoO Inferences" %in% sop()) {
            write_csv(reduce(rpe[[i]][[3]], rbind), file = inference_filename)
          }

          out_temp <- list()
          if ("Plots" %in% sop()) {
            for (j in 1:length(rpe[[i]][[1]])) {
              out_temp[[j]] <- ggpubr::ggarrange(plotlist = list(rpe[[i]][[1]][[j]], 
                                                                 NULL,
                                                                 rpe[[i]][[2]][[j]] + theme(plot.title = element_blank())), ncol = 1, heights = c(1,0.08,1), common.legend = T, align = "v", legend = "right")
            }
            ggpubr::ggexport(plotlist = out_temp, filename = plot_filename, verbose = FALSE, width = wid_mod, height = hei_mod)
          }
        }
      }
    }

    recfig <- list()
    show_chr_in <- show_ch()
    if (length(rpe[[show_chr_in]]) == 3) {
      recfig[[1]] <- renderPlot({ggpubr::ggarrange(plotlist = rpe[[show_chr_in]][[1]], ncol = 1, align = "v")}, height = hei * 0.9)
      recfig[[2]] <- renderPlot({ggpubr::ggarrange(plotlist = rpe[[show_chr_in]][[2]], ncol = 1, align = "v")}, height = hei * 0.9)
    } else if (length(rpe[[show_chr_in]]) > 3) {
      if (length(rpe[[show_chr_in]][[2]][[1]]) == 1) {
        for (i in 1:length(rpe[[show_chr_in]][[2]])) {
          rpe[[show_chr_in]][[2]][[i]] <- ggpubr::ggarrange(plotlist =  rpe[[show_chr_in]][[2]][[i]], ncol = 1)
        }
      } else {
        for (i in 1:length(rpe[[show_chr_in]][[2]])) {
          plot_scale <- rpe[[show_chr_in]][[2]][[i]][[3]]
          part1 <- ggpubr::ggarrange(plotlist = rpe[[show_chr_in]][[2]][[i]][[1]], ncol = 1, align = "v", heights = c(plot_scale + 0.02, 1 - 2*plot_scale - 0.02, plot_scale))
          part2 <- ggpubr::ggarrange(plotlist = rpe[[show_chr_in]][[2]][[i]][[2]], ncol = 1, align = "v", heights = c(plot_scale + 0.02, 1 - 2*plot_scale - 0.02, plot_scale))
          rpe[[show_chr_in]][[2]][[i]] <- ggpubr::ggarrange(part1, part2, ncol = 1)
        }
      }
      recfig[[1]] <- renderPlot({ggpubr::ggarrange(plotlist = rpe[[show_chr_in]][[1]], ncol = 1, align = "v")}, height = hei * 0.9)
      recfig[[2]] <- renderPlot({ggpubr::ggarrange(plotlist = rpe[[show_chr_in]][[2]], ncol = 1, align = "v")}, height = hei)
      recfig[[3]] <- DT::renderDataTable(reduce(rpe[[show_chr_in]][[3]], rbind),
                                         options = list(
                                           pageLength = 25,
                                           searchHighlight = TRUE
                                         ),
                                         filter = 'top')
      recfig[[4]] <- renderPlot({ggpubr::ggarrange(plotlist = rpe[[show_chr_in]][[4]], ncol = 1, align = "v")}, height = hei * 0.9)
    }

    output$ui_out <- renderUI({
      nTabs <- length(recfig)
      if (nTabs != 0) {
        if (nTabs == 1) {
          output$plot1 <- recfig[[1]]
          tablist <- list(tabPanel(title = "Message", plotOutput("plot1", height = paste0(hei,'px'))))
        } else if (nTabs == 2) {
          output$plot1 <- recfig[[1]]
          output$plot2 <- recfig[[2]]
          tablist <- list(tabPanel(title = "GoO inference", plotOutput("plot1", height = paste0(hei,'px'))),
                          tabPanel(title = "Density", plotOutput("plot2", height = paste0(hei,'px'))))
        } else if (nTabs > 2) {
          output$plot1 <- recfig[[1]]
          output$plot2 <- recfig[[2]]
          output$plot3 <- recfig[[3]]
          output$plot4 <- recfig[[4]]
          tablist <- list(tabPanel(title = "GoO inference", plotOutput("plot1", height = paste0(hei,'px'))),
                          tabPanel(title = "Location inference", plotOutput("plot2", height = paste0(hei,'px'))),
                          tabPanel(title = "Location table", DT::DTOutput("plot3", height = '450px')),
                          tabPanel(title = "Density", plotOutput("plot4", height = paste0(hei,'px')))
                          )
        }
        do.call(tabsetPanel, tablist)
      }
    })
  })
}
