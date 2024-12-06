#' Plotting function for S3 class "RecView"
#' @description Plotting function for S3 class "RecView".
#' @usage plot(obj)
#' @param obj S3 Class "RecView".
#'
#' @export
plot.RecView <- function(obj) {
  cat('Plotting...')
  
  offspring <- pull(obj$offspring, "Offspring")
  comparison_result <- obj$pairwise_comparison %>% 
    filter(!is.na(Concord_Discord))
  
  group_tb <- tibble(group = unique(comparison_result$group)) %>% 
    separate(col = "group", into = c("a_group", "b_group"), sep = "_<>_", remove = FALSE) %>% 
    select(group, a_group, b_group)
  
  comparison_result_pat <- comparison_result %>% 
    filter(Chromosome_origin == "Paternal") %>% 
    left_join(group_tb %>% `colnames<-`(.,c("group", "x_group", "y_group")), by = "group")
  
  comparison_result_mat <- comparison_result %>% 
    filter(Chromosome_origin == "Maternal") %>% 
    left_join(group_tb %>% `colnames<-`(.,c("group", "y_group", "x_group")), by = "group")
  
  comparison_result <- bind_rows(comparison_result_pat, comparison_result_mat) %>% 
    distinct(Chromosome_origin, POS_chr, x_group, y_group, Concord_Discord, .keep_all = TRUE) %>% 
    arrange(Chromosome_origin, POS_chr) %>% 
    nest(input = !c("x_group", "y_group", "Chromosome_origin"))
  
  for (i in 1:nrow(comparison_result)) {
    if (comparison_result$Chromosome_origin[i] == "Paternal") {
      x_group_tmp <- comparison_result$x_group[i]
      y_group_tmp <- comparison_result$y_group[i]
      
      if (which(offspring == x_group_tmp) > which(offspring == y_group_tmp)) {
        comparison_result$x_group[i] <- y_group_tmp
        comparison_result$y_group[i] <- x_group_tmp
      }
    } else {
      x_group_tmp <- comparison_result$x_group[i]
      y_group_tmp <- comparison_result$y_group[i]
      
      if (which(offspring == x_group_tmp) < which(offspring == y_group_tmp)) {
        comparison_result$x_group[i] <- y_group_tmp
        comparison_result$y_group[i] <- x_group_tmp
      }
    }
  }
  comparison_result <- comparison_result %>% 
    unnest(cols = "input")
  comparison_result$CHROM <- factor(comparison_result$CHROM, levels = unique(comparison_result$CHROM))
  comparison_result$Chromosome_origin <- factor(comparison_result$Chromosome_origin, levels = c("Paternal", "Maternal"))
  comparison_result$x_group <- factor(comparison_result$x_group, levels = offspring)
  comparison_result$y_group <- factor(comparison_result$y_group, levels = offspring)
  
  scaffold_indication <- comparison_result %>% 
    unite(col = "group", c("x_group", "y_group"), sep = "_<>_") %>% 
    summarise(min = min(POS_chr), max = max(POS_chr), .by = c("CHROM", "group")) %>% 
    gather(key = "minmax", value = "POS_chr", c("min", "max")) %>% 
    separate(col = "group", into = c("x_group", "y_group"), sep = "_<>_")
  scaffold_indication$CHROM <- factor(scaffold_indication$CHROM, levels = unique(comparison_result$CHROM))
  scaffold_indication$x_group <- factor(scaffold_indication$x_group, levels = offspring)
  scaffold_indication$y_group <- factor(scaffold_indication$y_group, levels = offspring)
  
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
  
  extra_tb <- tibble(x_group = offspring, y_group = offspring, 
                     x = list(tibble(x = c(0, x_axis_max))), y = list(tibble(y = c(1.3, -0.5)))) %>% 
    unnest(cols = c("x", "y"))
  extra_tb$x_group <- factor(extra_tb$x_group, levels = offspring)
  extra_tb$y_group <- factor(extra_tb$y_group, levels = offspring)
  
  theme_custom <- theme(panel.grid.major.y = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.background = element_rect(fill = "grey95"),
                        strip.text = element_text(face = "bold"),
                        axis.title.y = element_blank())
  # <<<<
  
  p <- ggplot() +
    scattermore::geom_scattermore(data = comparison_result, aes(x = POS_chr/1e6, y = Concord_Discord, color = Chromosome_origin), pointsize = point_size, position = position_jitter(height = 0.2), pixels = c(1500, 264)) +
    geom_line(data = scaffold_indication, aes(x = POS_chr/1e6, y = -0.5, color = CHROM, group = 1), linewidth = 4) +
    geom_line(data = extra_tb, aes(x = x, y = y, group = 1)) +
    scale_x_continuous(breaks = seq(0, x_axis_max, break_step),
                       limits = c(-x_axis_max*0.01, x_axis_max*1.01), 
                       expand = c(0,0)) +
    scale_y_continuous(breaks = c(-0.5, 0,1), labels = c("Scaffold", "Concord", "Discord")) +
    scale_color_manual(values = c("black", "#AEB6E5", colour_scheme)) +
    guides(color = guide_legend(override.aes = list(size=8))) +
    labs(x = "Position (Mb)", color = "Legend") +
    facet_grid(rows = vars(y_group), cols = vars(x_group)) +
    theme_bw() +
    theme_custom
  
  return(p)
}
