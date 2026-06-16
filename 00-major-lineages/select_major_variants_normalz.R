df <- read.csv("freq_by_day.csv")
ave_freq <- colSums(df[, 3:ncol(df)], na.rm = TRUE) / nrow(df)
ave_sort <- sort(ave_freq, decreasing = TRUE)
cum_freq <- cumsum(ave_sort)

threshold <- 0.90
keep_90 <- cum_freq <= threshold

cum_freq_90 <- cum_freq[keep_90]
ave_sort_90 <- ave_sort[keep_90]

df_90 <- data.frame(variants = names(cum_freq_90), mean_freq = ave_sort_90, cum_mean_freq = cum_freq_90)

write.csv(df_90, "major_variants.csv", row.names = FALSE)

df_num <- data.frame(variants = names(ave_sort), mean_freq = ave_sort, cum_mean_freq = cum_freq, num = seq_along(ave_sort))

major_num <- nrow(df_90)
major_y <- df_num$cum_mean_freq[major_num]

x_ticks <- sort(unique(c(pretty(df_num$num), major_num)))

# plot selected major variants
tiff("S2_Fig.tif", width = 6, height = 4.5, units = "in", res = 300, compression = "lzw", bg = "white")
plot(df_num$num, df_num$cum_mean_freq, type = "p" ,col = "black", xlab = "Number of lineages (ranked by sequence counts)",
     ylab = "Cumulative percentage of sequqnce counts (%)", ylim = c(0, 1), xaxt = "n", yaxt = "n")
axis(side = 1, at = x_ticks, labels = x_ticks)
axis(side = 2, at = seq(0, 1, by = 0.1), labels = paste0(seq(0, 100, by = 10)))
abline(h = threshold, lty = 2)
abline(v = major_num, lty = 2)

dev.off()
