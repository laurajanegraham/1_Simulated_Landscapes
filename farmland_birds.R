# plot and condense the output from the simulations

library(plyr)
library(tidyr)
library(dplyr)

if(!dir.exists("farmland_birds/plots")) dir.create("farmland_birds/plots")

results <- ldply(list.files("farmland_birds/results/", pattern="output", full.names=TRUE), read.csv) %>%
  group_by(h_val, npp, w1, w2) %>%
  summarise(es_mean_se = sd(es_mean)/sqrt(n()),
            es_mean = mean(es_mean),
            es_var_se = sd(es_var)/sqrt(n()),
            es_var = var(es_var)) %>%
  
fs <- list.files("farmland_birds/results/", full.names = TRUE)
file.remove(fs)
save(results, file="farmland_birds/results/results.rda")

