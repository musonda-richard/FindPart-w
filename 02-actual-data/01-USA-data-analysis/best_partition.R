library(dplyr)

df <- read.csv("best_estimates.csv")
df_selected <- select(df,variant, k, k_lb, k_ub)
df_sig <- mutate(df_selected, across(where(is.numeric), ~ floor(.x * 1000) / 1000))
df_connected <- mutate(df_sig, k_ci_95 = paste0(k_lb, "—", k_ub)) %>% select(-c(k_lb,k_ub))
df_grouped <- group_by(df_connected, k, k_ci_95)
df_variant_combined <- summarise(df_grouped, variant = paste(variant, collapse = ","),.groups = "drop")
df_variant_combined_sort <- arrange(df_variant_combined, desc(k)) %>% 
  select(variant, k, k_ci_95)

write.csv(df_variant_combined_sort, "best_partition.csv", row.names = F)
