AIC_RelRe <- read.csv("../01-USA-data-analysis/loglikelihood.csv")$AIC[1]; num_mle_RelRe <- 1

AIC_findpart <- read.csv("../01-USA-data-analysis/rgs-AIC.csv")[1,5]; num_mle_findpart <- nrow(read.csv("../01-USA-data-analysis/rgs-AIC.csv"))

AIC_kmeans19 <- read.csv("../02-kmeans/01-kmean-19/kmeans-rgs-AIC.csv")[1,6]; num_mle_kmeans19 <- nrow(read.csv("../02-kmeans/01-kmean-19/kmeans-rgs-AIC.csv"))

AIC_kmeans100 <- read.csv("../02-kmeans/02-kmean-100/kmeans-rgs-AIC.csv")[1,6]; num_mle_kmeans100 <- nrow(read.csv("../02-kmeans/02-kmean-100/kmeans-rgs-AIC.csv"))

AIC_hierarchical <- read.csv("../03-hierarchical/hierarchical-rgs-AIC.csv")[1,7]; num_mle_hierarchical <- nrow(read.csv("../03-hierarchical/hierarchical-rgs-AIC.csv"))

df <- data.frame(Method = c("FindPart-1", "Centroid-linkage hierarchical", "k-means-19", "k-means-100", "Original RelRe"), Calculations = c(num_mle_findpart, num_mle_hierarchical, num_mle_kmeans19, num_mle_kmeans100, num_mle_RelRe), AIC = c(AIC_findpart, AIC_hierarchical, AIC_kmeans19, AIC_kmeans100, AIC_RelRe))

tiff("method_comparison.tiff", width = 7, height = 5, units = "in", res = 300)

plot(x = df$AIC, y = df$Calculations, xlab = "AIC", ylab = "Number of MLEs", pch = 19, cex = 1.5, las = 1, xlim = c(min(df$AIC) - 1, max(df$AIC) + 2), ylim = c(0, max(df$Calculations) * 1.1), xaxt = "n")

axis(side = 1, at = sort(unique(c(pretty(df$AIC), AIC_findpart, AIC_RelRe))), labels = round(sort(unique(c(pretty(df$AIC), AIC_findpart, AIC_RelRe))), 2), cex.axis = 0.8)

text(x = df$AIC, y = df$Calculations, labels = df$Method, pos = 3, cex = 0.8)

box(lwd = 2)

dev.off()