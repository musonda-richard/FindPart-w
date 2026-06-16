library(ggplot2)
library(tidyr)
library(dplyr)

sha_best <- as.character(read.csv("kmeans-rgs-AIC.csv")[1,3])
df_best_estimates <- read.csv(paste0(sha_best, "_estimates.csv")) %>% arrange(desc(k))

lineage_order <- rev(df_best_estimates$variant)

k <- as.vector(read.csv("kmeans-rgs-AIC.csv")[1,1])
df <- read.csv(paste0("consensus_k", k, ".csv"))
df_long <- pivot_longer(df, cols = -Lineage, names_to = "Lineage2", values_to = "Value") %>% 
  rename(Lineage1=Lineage)

# apply ordering
df_long$Lineage1 <- factor(df_long$Lineage1, levels = lineage_order)
df_long$Lineage2 <- factor(df_long$Lineage2, levels = lineage_order)

p <- ggplot(df_long, aes(x = Lineage2, y = Lineage1, fill = Value)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
  coord_fixed(ratio = 1) +
  labs(x = "Lineage", y = "Lineage", fill = "Consensus") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"), panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA))

ggsave("consensus_matrix.tiff", plot = p, device = "tiff", width = 12, height = 8, dpi = 600)

