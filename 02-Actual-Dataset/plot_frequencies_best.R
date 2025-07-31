require(dplyr)
baseline = "KP.3.1.1"
df_frequency <- read.csv("best_frequencies.csv") %>% select(-c(average_c,average_k))
df_frequency$date <- as.Date(df_frequency$date)
df_frequency <- select(df_frequency, date, baseline, everything())
dates <- df_frequency$date
n <- ncol(df_frequency) - 1
color <- rep(rainbow(10),10)

df_count <- read.csv("count_variants.csv") %>% select(-date_till)
df_count$date_from <- as.Date(df_count$date_from)
variant_columns <- df_count[ ,-1]
n <- ncol(variant_columns)
f_counts <- variant_columns / rowSums(variant_columns, na.rm = T)
f_counts <- cbind(date=df_count[,1], f_counts)
f_counts <- select(f_counts, date, baseline, everything())
f_counts$date <- as.Date(f_counts$date)
weeks <- f_counts$date

pdf("trajectory_best.pdf",height = 4, width = 6)
plot(dates, rep(0,length(dates)), xlab="", xaxt="n", type="n", ylab="Relative frequency",
     ylim=c(0,1), las = 2, main = "")
axis(side=1, at=dates, labels=format(dates,"%h %d"), tick=T, las=2)

pch_values <- rep(0:18,10)

for (i in 1:n) {
  points(weeks, f_counts[,i+1], lty=1, lwd=1, col=color[i], cex=0.6, pch=pch_values[i])
  lines(dates, df_frequency[,i+1], lty=1, lwd=1, col=color[i], cex=0.8, pch=pch_values[i])
}
legend("top", legend=colnames(df_frequency)[2:(n+1)], col=color[1:n], lty=1,
       lwd=1.0, bty="n", cex=0.8, ncol=4, pch=pch_values[1:n])
dev.off()
