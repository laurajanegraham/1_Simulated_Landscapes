# plot and condense the output from the simulations

library(plyr)
library(tidyr)
library(dplyr)

results <- ldply(list.files("farmland_birds/results/", pattern="output", full.names=TRUE), read.csv) %>%
  group_by(h_val, npp, w1, w2) %>%
  summarise(es_mean_se = sd(es_mean)/sqrt(n()),
            es_mean = mean(es_mean),
            es_var_se = sd(es_var)/sqrt(n()),
            es_var = var(es_var))
  
fs <- list.files("farmland_birds/results/", full.names = TRUE)
file.remove(fs)
save(results, file="farmland_birds/results/results.rda")


# plots
load("farmland_birds/results/results.rda")
library(lattice)

results <- mutate(results, window_size = paste0("window = ", w1))
results$w1 <- factor(results$w1, levels=c("3", "9", "15"))
pdf(file=paste0("farmland_birds/plots/npp_vs_config.pdf"),  
    width=15, 
    height=5)
print(wireframe(es_mean~h_val*npp|w1, zlab = list("ES", rot=90), drape = TRUE, data=results, col.regions = colorRampPalette(c("blue", "pink"))(100)))
dev.off()

