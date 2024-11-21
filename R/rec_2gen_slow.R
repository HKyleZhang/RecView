#' Locate recombination with two-generation pedigree
#' @description rec_2gen() old codes; very slow...
#' 
#' @noRd
rec_2gen_slow <- function(data, sc_order, chromosome, offspring) {
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
  
  result_collection <- tibble(off_1 = as.character(), off_2 = as.character(), off_3 = as.character())
  for(i in 1:(length(offspring)-2)) {
    for (j in (i+1):(length(offspring)-1) ) {
      for (k in (j+1):length(offspring)){
        result_collection <- add_row(result_collection, off_1 = offspring[i], off_2 = offspring[j], off_3 = offspring[k])
      }
    }
  }
  result_collection_pat <- result_collection %>% 
    mutate(id = seq(1,nrow(.))) %>% 
    nest(input = !id) %>% 
    mutate(output = map(input, trioff_analyse, data_chr, "Paternal")) %>% 
    unnest(cols = "output") %>% 
    select(-id)
  
  result_collection_mat <- result_collection %>% 
    mutate(id = seq(1,nrow(.))) %>% 
    nest(input = !id) %>% 
    mutate(output = map(input, trioff_analyse, data_chr, "Maternal")) %>% 
    unnest(cols = "output") %>% 
    select(-id)
  
  result_collection <- bind_rows(result_collection_pat, result_collection_mat)
  
  comparison_result <- reduce(result_collection$comparison_result, bind_rows) %>% 
    distinct(Side, POS_chr, x_group, y_group, Chromosome, .keep_all = TRUE) %>% 
    arrange(Side, POS_chr)
  comparison_result$CHROM <- factor(comparison_result$CHROM, levels = unique(comparison_result$CHROM))
  comparison_result$Side <- factor(comparison_result$Side, levels = c("Paternal", "Maternal"))
  comparison_result$x_group <- factor(comparison_result$x_group, levels = offspring)
  comparison_result$y_group <- factor(comparison_result$y_group, levels = offspring)
  
  scaffold_indication <- comparison_result %>% 
    unite(col = "group", c("x_group", "y_group"), sep = "_<>_") %>% 
    summarise(min = min(POS_chr), max = max(POS_chr), .by = c("CHROM", "group")) %>% 
    gather(key = "minmax", value = "POS_chr", c("min", "max")) %>% 
    separate(col = "group", into = c("x_group", "y_group"), sep = "_<>_")
  scaffold_indication$CHROM <- factor(scaffold_indication$CHROM, levels = unique(comparison_result$CHROM))
  
  # Plotting setting-related >>>>
  macaron_mod <- c("#940214", "#3e1f16", "#123358", "#e8c001", "#015c7a", "#85ab87", "#d75624", "#dd395a")
  n_scaffold <- length(unique(comparison_result$CHROM))
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
  colour_scheme <- colour_scheme[1:n_scaffold]
  
  if (nrow(comparison_result) == 0) {
    x_axis_max <- 0
  } else {
    x_axis_max <- max(comparison_result$POS_chr)/1e6
  }
  if (x_axis_max > 15) {
    break_step <- 10
  } else if (x_axis_max > 5) {
    break_step <- 2.5
  } else if (x_axis_max > 2) {
    break_step <- 1
  } else if (x_axis_max > 1) {
    break_step <- 0.25
  }
  
  if (nrow(comparison_result) < 2000) {
    point_size <- 6
  } else {
    point_size <- exp(-8.407e-08 * nrow(comparison_result) + 0.9498) %>% round(1)
  }
  
  theme_custom <- theme(panel.grid.major.y = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.background = element_rect(fill = "grey95"),
                        strip.text = element_text(face = "bold"),
                        axis.title.y = element_blank())
  # <<<<
  
  p <- ggplot() +
    scattermore::geom_scattermore(data = comparison_result, aes(x = POS_chr/1e6, y = Chromosome, color = Side), pointsize = point_size, position = position_jitter(height = 0.2), pixels = c(1500, 264)) +
    geom_line(data = scaffold_indication, aes(x = POS_chr/1e6, y = -0.5, color = CHROM, group = 1), linewidth = 4) +
    scale_x_continuous(breaks = seq(0, x_axis_max, break_step),
                       limits = c(-x_axis_max*0.01, x_axis_max*1.01), 
                       expand = c(0,0)) +
    scale_y_continuous(breaks = c(-0.5, 0,1), labels = c("Scaffold", "Concord", "Discord")) +
    scale_color_manual(values = c("grey15", "grey45", colour_scheme)) +
    guides(color = guide_legend(override.aes = list(size=8))) +
    labs(x = "Position (Mb)", color = "Legend") +
    facet_grid(rows = vars(y_group), cols = vars(x_group)) +
    theme_bw() +
    theme_custom
  
  # p_pat <- ggplot() +
  #   scattermore::geom_scattermore(data = comparison_result %>% filter(Side == "Paternal"), aes(x = POS_chr/1e6, y = Chromosome, color = CHROM), pointsize = point_size, position = position_jitter(height = 0.2), pixels = c(1500, 264)) +
  #   scale_x_continuous(breaks = seq(0, x_axis_max, break_step),
  #                      limits = c(-x_axis_max*0.01, x_axis_max*1.01), 
  #                      expand = c(0,0)) +
  #   scale_y_continuous(breaks = c(0,1), labels = c("Concord", "Discord")) +
  #   scale_color_manual(values = colour_scheme) +
  #   guides(color = guide_legend(override.aes = list(size=8))) +
  #   labs(x = "Position (Mb)", color = "Scaffold") +
  #   facet_grid(rows = vars(y_group), cols = vars(x_group)) +
  #   theme_bw() +
  #   theme_custom
  # 
  # p_mat <- ggplot() +
  #   scattermore::geom_scattermore(data = comparison_result %>% filter(Side == "Maternal"), aes(x = POS_chr/1e6, y = Chromosome, color = CHROM), pointsize = point_size, position = position_jitter(height = 0.2), pixels = c(1500, 264)) +
  #   scale_x_continuous(breaks = seq(0, x_axis_max, break_step),
  #                      limits = c(-x_axis_max*0.01, x_axis_max*1.01), 
  #                      expand = c(0,0)) +
  #   scale_y_continuous(breaks = c(0,1), labels = c("Concord", "Discord")) +
  #   scale_color_manual(values = colour_scheme) +
  #   guides(color = guide_legend(override.aes = list(size=8))) +
  #   labs(x = "Position (Mb)", color = "Scaffold") +
  #   facet_grid(rows = vars(y_group), cols = vars(x_group)) +
  #   theme_bw() +
  #   theme_custom
  # 
  # p <- ggpubr::ggarrange(ggpubr::text_grob(label = "Paternal\nchromosome", rot = 90, size = 15,face = "bold"), p_pat,
  #                        ggpubr::text_grob(label = "Maternal\nchromosome", rot = 90, size = 15,face = "bold"), p_mat,
  #                        ncol = 2, nrow = 2,
  #                        widths = c(0.1,1))
  
  position_result <- reduce(result_collection$position_result, bind_rows) %>%
    arrange(Offspring, Recombination_Position) %>% 
    summarise(Mean_bp = mean(Recombination_Position), SD_bp = sd(Recombination_Position), .by = c("Offspring", "id", "Side")) 
  
  position_result_dup <- position_result %>% 
    group_by(Mean_bp, Side) %>% 
    tally() %>% 
    filter(n > 1)
  
  position_result_pat <- position_result %>% 
    filter(Side == "Paternal") %>% 
    filter(!(Mean_bp %in% position_result_dup$Mean_bp[position_result_dup$Side == "Paternal"]))
  
  position_result_mat <- position_result %>% 
    filter(Side == "Maternal") %>% 
    filter(!(Mean_bp %in% position_result_dup$Mean_bp[position_result_dup$Side == "Maternal"]))
  
  position_result <- bind_rows(position_result_pat, position_result_mat) %>% 
    select(Side, Offspring, Mean_bp, SD_bp)
  
  return(list(p, position_result))
}
#' 
#' @noRd
trioff_analyse <- function(input, gt_data_in, side) {
  offspring_in <- as.data.frame(t(input)) %>% pull(V1)
  comparison_labels <- c(paste(offspring_in[1], offspring_in[2], sep = "_<>_"),
                         paste(offspring_in[1], offspring_in[3], sep = "_<>_"),
                         paste(offspring_in[2], offspring_in[3], sep = "_<>_"))
  
  if (side == "Paternal") {
    informative_gts <- c("1_0_0", "1_0_1", "1_2_1", "1_2_2")
  } else {
    informative_gts <- c("0_1_0", "0_1_1", "2_1_1", "2_1_2")
  }
  
  data_mod <- gt_data_in %>% 
    select(id, CHROM, orientation, POS, POS_chr, AB, CD, all_of(offspring_in)) %>% 
    arrange(POS_chr) %>% 
    unite(col = "GTS_1", c("AB", "CD", offspring_in[1]), remove = FALSE) %>% 
    unite(col = "GTS_2", c("AB", "CD", offspring_in[2]), remove = FALSE) %>%
    unite(col = "GTS_3", c("AB", "CD", offspring_in[3]), remove = FALSE) %>% 
    filter(GTS_1 %in% informative_gts, 
           GTS_2 %in% informative_gts,
           GTS_3 %in% informative_gts) %>% 
    mutate("{comparison_labels[1]}" := abs(get(offspring_in[1],.) - get(offspring_in[2],.)),
           "{comparison_labels[2]}" := abs(get(offspring_in[1],.) - get(offspring_in[3],.)),
           "{comparison_labels[3]}" := abs(get(offspring_in[2],.) - get(offspring_in[3],.)))
  
  cp12 <- changepoint::cpt.mean(pull(data_mod, comparison_labels[1]))@cpts
  cp13 <- changepoint::cpt.mean(pull(data_mod, comparison_labels[2]))@cpts
  cp23 <- changepoint::cpt.mean(pull(data_mod, comparison_labels[3]))@cpts
  
  cp12 <- cp12[which(cp12 != nrow(data_mod))]
  cp13 <- cp13[which(cp13 != nrow(data_mod))]
  cp23 <- cp23[which(cp23 != nrow(data_mod))]
  
  cp1 <- tibble(Offspring = offspring_in[1], Recombination_Position = data_mod$POS_chr[cp12[cp12 %in% cp13]]) %>% 
    arrange(Recombination_Position)
  if (nrow(cp1) > 0) {
    cp1 <- cp1 %>% mutate(id = seq(1,nrow(.)))
  } else {
    cp1 <- mutate(cp1, id = as.numeric())
  }
  
  cp2 <- tibble(Offspring = offspring_in[2], Recombination_Position = data_mod$POS_chr[cp12[cp12 %in% cp23]]) %>% 
    arrange(Recombination_Position)
  if (nrow(cp2) > 0) {
    cp2 <- cp2 %>% mutate(id = seq(1,nrow(.)))
  } else {
    cp2 <- mutate(cp2, id = as.numeric())
  }
  
  cp3 <- tibble(Offspring = offspring_in[3], Recombination_Position = data_mod$POS_chr[cp13[cp13 %in% cp23]]) %>% 
    arrange(Recombination_Position)
  if (nrow(cp3) > 0) {
    cp3 <- cp3 %>% mutate(id = seq(1,nrow(.)))
  } else {
    cp3 <- mutate(cp3, id = as.numeric())
  }
  
  rec_pos <- bind_rows(cp1, cp2, cp3) %>% 
    mutate(Side = side)
  
  if (side == "Paternal") {
    data_mod <- data_mod %>% 
      mutate(Side = side) %>% 
      gather(key = "group", value = "Chromosome", all_of(comparison_labels)) %>% 
      separate(col = "group", into = c("x_group", "y_group"), sep = "_<>_")
  } else {
    data_mod <- data_mod %>% 
      mutate(Side = side) %>% 
      gather(key = "group", value = "Chromosome", all_of(comparison_labels)) %>% 
      separate(col = "group", into = c("y_group", "x_group"), sep = "_<>_")
  }
  
  out <- tibble(comparison_result = list(data_mod), position_result = list(rec_pos))
  return(out)
}