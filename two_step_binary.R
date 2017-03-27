# plot and condense the output from the simulations

library(plyr)
library(tidyr)
library(dplyr)
library(lattice)

if(!dir.exists("two_step_binary/plots")) dir.create("two_step_binary/plots")

results <- ldply(list.files("two_step_binary/results/", pattern="output", full.names=TRUE), read.csv) %>%
  group_by(h_val, p_val, w1, w2, f1, f2) %>%
  summarise(es_mean_se = sd(es_mean)/sqrt(n()),
            es_mean = mean(es_mean),
            es_var_se = sd(es_var)/sqrt(n()),
            es_var = var(es_var)) %>%
  mutate(window_size = paste0("w1 = ", w1, ", w2 = ", w2),
         func = paste0("f1 = ", f1, ", f2 = ", f2))

for(w in unique(results$window_size)) {
  wfile = gsub(" = ", "_", w)
  wfile = gsub(", ", "_", wfile)
  res <- filter(results, window_size == w)
  pdf(file=paste0("two_step_binary/plots/", wfile, ".pdf"),  
      width=10, 
      height=10)
  print(wireframe(es_mean~h_val*p_val|func, zlab = list("ES", rot=90), drape = TRUE, data=res, col.regions = colorRampPalette(c("blue", "pink"))(100)))
  dev.off()
}

fs <- list.files("two_step_binary/results/", full.names = TRUE)
file.remove(fs)
save(results, file="two_step_binary/results/results.rda")

