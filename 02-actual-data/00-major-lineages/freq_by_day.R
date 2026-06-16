require(tidyverse)
require(lubridate)
start <- "2024-10-01"
end <- "2024-12-31" 

df_date_spread <- read.csv("full_variant_counts.csv")
df_date_spread$date <- as.Date(df_date_spread$date)
df_date_spread <- filter(df_date_spread, date >= as.Date(start) & date <= as.Date(end))
df_day_from <- rename(df_date_spread, date_from = date)
df_day_to <- mutate(df_day_from, date_till = date_from) %>%
  select(date_from, date_till, everything())

# normalization 
vec_day_totals <- rowSums(df_day_to[,3:ncol(df_day_to)])
normalised_df <- cbind(df_day_to[,1:2],
                       df_day_to[,3:ncol(df_day_to)]/ vec_day_totals)

write.csv(normalised_df, "freq_by_day.csv", row.names = F)
