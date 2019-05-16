# examine H statistics for double copy vs. single copy loci,
# i.e. for miscanthus tags aligned to sorghum reference vs. miscanthus
# reference.
# Can we use this statistic to distinguish good loci from paralogs?
load("190515counts_matrices.RData")
tagtab <- read.csv("190515paralog_tags.csv", stringsAsFactors = FALSE)
Rcpp::sourceCpp('simpson.cpp')
library(ggplot2)
library(viridis)

# set up alleles2loc vectors
alleles2locM <- match(tagtab$Miscanthus, unique(tagtab$Miscanthus))
dim(diploid_mat)
length(alleles2locM)
alleles2locM[1:20]
tagtab[1:20,2:3]

alleles2locS <- match(tagtab$Sorghum, unique(tagtab$Sorghum))
alleles2locS[1:20]

# compare miscanthus and sorghum for diploids
Depth_dip_M <- tapply(colSums(diploid_mat), alleles2locM, sum)

H_dip_Misc <- HReadDepth(diploid_mat, alleles2locM, max(alleles2locM))
Hind_dip_Misc <- colMeans(H_dip_Misc[[1]], na.rm = TRUE)

ggplot(mapping = aes(x = H_dip_Misc[[2]], y = Hind_dip_Misc,
                     col = log(Depth_dip_M))) +
  geom_point() +
  scale_color_viridis(option = "magma") +
  geom_abline(slope = 1/2, intercept = 0) +
  labs(x = "He", y = "Hind", title = "Diploids with Miscanthus reference")

hist(log2(Hind_dip_Misc/H_dip_Misc[[2]]))
mean(Hind_dip_Misc/H_dip_Misc[[2]] > 1/2, na.rm = TRUE) # 23% exceed expectations (30% when corrected for size)

Depth_dip_S <- tapply(colSums(diploid_mat), alleles2locS, sum)

H_dip_Sorg <- HReadDepth(diploid_mat, alleles2locS, max(alleles2locS))
Hind_dip_Sorg <- colMeans(H_dip_Sorg[[1]], na.rm = TRUE)

ggplot(mapping = aes(x = H_dip_Sorg[[2]], y = Hind_dip_Sorg,
                     col = log(Depth_dip_S))) +
  geom_point() +
  scale_color_viridis(option = "magma") +
  geom_abline(slope = 1/2, intercept = 0) +
  labs(x = "He", y = "Hind", title = "Diploids with Sorghum reference")

hist(log2(Hind_dip_Sorg/H_dip_Sorg[[2]]))
mean(Hind_dip_Sorg/H_dip_Sorg[[2]] > 1/2, na.rm = TRUE) # 32% exceed expectations (38% previously)

# put into data frame to compare distributions
dip_df <- data.frame(Reference = c(rep("Miscanthus", length(Hind_dip_Misc)),
                                   rep("Sorghum", length(Hind_dip_Sorg))),
                     He = c(H_dip_Misc[[2]], H_dip_Sorg[[2]]),
                     Hind = c(Hind_dip_Misc, Hind_dip_Sorg),
                     Depth = c(Depth_dip_M, Depth_dip_S))
dip_df <- dip_df[which(dip_df$Depth > 0 & dip_df$He > 0),]

ggplot(dip_df, aes(x = Depth, y = Hind/He, color = Reference)) +
  geom_point() +
#  geom_density2d(color = "black") +
  geom_smooth(color = "black") +
  coord_trans(x = "log") +
  facet_wrap(~ Reference)

ggplot(dip_df, aes(x = log(Hind/He/log(Depth)), by = Reference,
                   color = Reference)) +
  geom_density()

ggplot(dip_df, aes(x = log2(Hind/He), by = Reference,
                   color = Reference)) +
  geom_density() +
  geom_vline(xintercept = -1)

# how to properly scale by depth?
# He is prob of sampling two different reads in population

# Hind_mat is taxa x locus
# Depth_mat is taxa x locus
# He is a vector by locus
scale_Hind <- function(Hind_mat, He, Depth_mat){
  divisor <- sweep(Depth_mat, 2, He, function(x, y) (x - 1)/x * y)
  return(Hind_mat / divisor)
}

scale_Hind(matrix(c(0, 0.5, 0.2, 0.3), nrow = 2, ncol = 2),
           c(0.25, 0.4),
           matrix(c(5, 3, 10, 13), nrow = 2, ncol = 2))

locDepth_dip_Misc <- t(rowsum(t(diploid_mat), alleles2locM))

scaled_Hind_dip_Misc <- scale_Hind(H_dip_Misc[[1]],
                                   H_dip_Misc[[2]],
                                   locDepth_dip_Misc)
hist(log2(colMeans(scaled_Hind_dip_Misc, na.rm = TRUE)), breaks = 50)

locDepth_dip_Sorg <- t(rowsum(t(diploid_mat), alleles2locS))

scaled_Hind_dip_Sorg <- scale_Hind(H_dip_Sorg[[1]],
                                   H_dip_Sorg[[2]],
                                   locDepth_dip_Sorg)

# If I have one read He is reduced to zero.
# If I have two reads, and He is 0.5, I get 0.5 * 0 + 0.5 * 0.5 which is 0.25
# If I have three reads, I get 1/4 * 0 + 3/4 * 0.4444 = 1/3

dip_df <- data.frame(Reference = c(rep("Miscanthus", length(Hind_dip_Misc)),
                                   rep("Sorghum", length(Hind_dip_Sorg))),
                     He = c(H_dip_Misc[[2]], H_dip_Sorg[[2]]),
                     Hind = c(Hind_dip_Misc, Hind_dip_Sorg),
                     Depth = c(Depth_dip_M, Depth_dip_S))

dip_df$Scaled_ratio <- c(colMeans(scaled_Hind_dip_Misc, na.rm = TRUE),
                         colMeans(scaled_Hind_dip_Sorg, na.rm = TRUE))

dip_df <- dip_df[which(dip_df$Depth > 0 & dip_df$He > 0),]

head(dip_df)

ggplot(dip_df, aes(x = log2(Scaled_ratio), by = Reference,
                   color = Reference)) +
  geom_density() +
  geom_vline(xintercept = -1)

ggplot(dip_df, aes(x = log(Depth), y = log2(Scaled_ratio))) +
  geom_point() +
  geom_density2d(col = "blue") +
  geom_hline(yintercept = -1, color = "red") + # expected max ratio
  geom_vline(xintercept = log(5 * nrow(diploid_mat)), col = "green") + # reasonable minimum depth
  facet_wrap(~ Reference)

ggplot(dip_df[dip_df$Depth >= 5 * nrow(diploid_mat),], 
       aes(x = log2(Scaled_ratio), by = Reference,
           color = Reference)) +
  geom_density() +
  geom_vline(xintercept = -1)
