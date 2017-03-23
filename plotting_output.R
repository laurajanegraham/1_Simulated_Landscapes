library(ggplot2)
library(readr)
library(cowplot)
library(dplyr)

plot_output <- function(modname, varname) {
  # read in python output
  es_mod <- read_csv(paste0("results/", modname, ".csv"))
  
  # plot the results
  # window_size needs to be factor
  es_mod$window_size = factor(es_mod$window_size, levels=c("1", "3", "5", "9", "15"))
  
  # set column to the required variable
  es_mod[,"varname"] <- es_mod[,varname]
  
  # get the means of each of the landscapes
  es_mod_summary <- group_by(es_mod, h_val, p_val, window_size) %>% 
    summarise(mean_val = mean(varname),
              min_val = min(varname),
              max_val = max(varname),
              se_val = sd(varname)/sqrt(n()))
  
  # plot heatmap of the es_means
  hm <- ggplot(data = es_mod_summary, aes(x=p_val, y=h_val, fill=mean_val)) + 
    geom_tile() + scale_fill_gradient2(varname, low = "white", high = "blue") + facet_wrap(~window_size)
  save_plot(paste0("plots/", modname, "_heatmap_", varname, ".png"), hm, base_width = 10, base_height = 10)
  
  # plot means and standard errors of all the runs with a panel per proportion,
  # h_val on the x axis and grouped by window size (is it reduced uncertainty with
  # reduced window size)
  se_plot <- ggplot(data = es_mod_summary, aes(x = h_val, y = mean_val, color = window_size)) + 
    geom_point() +
    geom_line() + 
    geom_errorbar(aes(ymin=mean_val - se_val, ymax=mean_val + se_val), width=.05) + 
    facet_wrap(~p_val) + ylab(paste0(varname, "(SE)"))
  save_plot(paste0("plots/", modname, "_seplot_", varname, ".png"), se_plot, base_width = 10, base_height = 10)
}

plot_output("es_mod_rec", "es_mean")
plot_output("es_mod_rec", "es_total")
plot_output("es_mod_rec", "es_var")
plot_output("es_mod_poll_agri", "es1_mean")
plot_output("es_mod_poll_agri", "es1_total")
plot_output("es_mod_poll_agri", "es1_var")
plot_output("es_mod_poll_agri", "es2_mean")
plot_output("es_mod_poll_agri", "es2_total")
plot_output("es_mod_poll_agri", "es2_var")

# example landscape plotting - takes ages
es_mod_example <- read_csv("es_mod_example.csv")
es_mod_example$LC_type <- ifelse(es_mod_example$value ==0, "Agriculture", "Forest")
# plot the landscapes and save to results folder
p <- ggplot(data = es_mod_example, aes(x = x, y = y, fill = LC_type)) + geom_tile() + facet_grid(h ~ p) + scale_fill_manual(values = c("orange", "darkgreen"))
save_plot("results/es_mod_example.png", p, base_height = 10)
rm(es_mod_example)