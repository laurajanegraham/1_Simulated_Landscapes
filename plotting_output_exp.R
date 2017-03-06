library(ggplot2)
library(readr)
library(cowplot)
library(dplyr)
library(tidyr)

plot_output <- function(modname, es_output, a) {
  es_mod <- filter(es_output, a_val == a)
  es_mod$window_size <- factor(es_mod$window_size, levels=c(3, 5, 9, 15))

  # get the means of each of the landscapes
  es_mod_summary <- group_by(es_mod, h_val, p_val, window_size) %>% 
    summarise(mean_val = mean(es_mean),
              min_val = min(es_mean),
              max_val = max(es_mean),
              se_val = sd(es_mean)/sqrt(n()))
  
  # plot heatmap of the es_means
  hm <- ggplot(data = es_mod_summary, aes(x=p_val, y=h_val, fill=mean_val)) + 
    geom_tile() + scale_fill_gradient2("Mean ES value", low = "white", high = "blue") + facet_wrap(~window_size)
  save_plot(paste0("plots/", modname, "_heatmap_exp_aval_", a, ".png"), hm, base_width = 10, base_height = 10)
  
  # plot means and standard errors of all the runs with a panel per proportion,
  # h_val on the x axis and grouped by window size (is it reduced uncertainty with
  # reduced window size)
  se_plot <- ggplot(data = es_mod_summary, aes(x = h_val, y = mean_val, color = window_size)) + 
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(ymin=mean_val - se_val, ymax=mean_val + se_val), width=.05) + 
    facet_wrap(~p_val) + ylab("Mean ES value (SE)")
  save_plot(paste0("plots/", modname, "_seplot_exp_aval_", a, ".png"), se_plot, base_width = 10, base_height = 10)
}

modname <- "es_mod_rec_exp"
es_mod <- read_csv(paste0("results/", modname, ".csv"))
es_mod$ex_keys <- gsub("\\{", "", es_mod$ex_keys)
es_mod$ex_keys <- gsub("\\}", "", es_mod$ex_keys)
es_mod <- separate(es_mod, ex_keys, into = c("lc", "a"), ",") %>%
  separate(a, into=c("a", "a_val"), ":") %>%
  select(-a)

es_mod$a_val <- factor(es_mod$a_val, labels=c("-1", "-5", "-50", "1", "~0", "5", "50"))

for(exp in unique(es_mod$a_val)) {
  plot_output(modname, es_mod, exp)
}
