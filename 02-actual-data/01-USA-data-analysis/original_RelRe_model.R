library(dplyr)

df <- read.csv("estimates.csv")
df_selected <- select(df,variant, k, k_lb, k_ub) 
df_sorted <- arrange(df_selected, desc(k))
df_sig <- mutate(df_sorted, across(where(is.numeric), ~ floor(.x * 1000) / 1000))
df_connect_c.i <-  mutate(df_sig,k_ci_95 = paste0(k_lb, "—", k_ub))
df_ci_selected <- select(df_connect_c.i,variant, k, k_ci_95)

write.csv(df_ci_selected, "original_RelRe_model.csv", row.names = F)
