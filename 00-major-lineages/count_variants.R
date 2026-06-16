require(tidyverse)
library(lubridate)
start <- ymd("2024-10-01")
end <- ymd("2024-12-31")

df_major_variants <- read.csv("major_variants.csv") %>%
  select(variants)
df_all_var_spread <- read.csv("full_variant_counts.csv")
df_all_var_spread$date <- ymd(df_all_var_spread$date)
df_all_var_spread <- filter(df_all_var_spread,date >= start & date <= end)

major_variants <- df_major_variants$variants

df_all_var_long <- pivot_longer(df_all_var_spread ,cols = -date, names_to = "variants", values_to = "counts")

# Mark non-major variants as "Others"
df_all_var_long <- mutate(df_all_var_long, variants = if_else(variants %in% major_variants, variants, "Others"))

df_others_filtered <- filter(df_all_var_long, variants != "Others")

df_wider <- pivot_wider(df_others_filtered, names_from = variants, values_from = counts)

df_date_renamed <- rename(df_wider, date_from = date)

df_date_till_added <- mutate(df_date_renamed, date_till = date_from) %>%
  select(date_from, date_till, everything())

write.csv(df_date_till_added, "count_variants.csv", row.names = F) 
