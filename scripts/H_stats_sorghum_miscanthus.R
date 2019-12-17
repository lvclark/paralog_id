# examine H statistics for double copy vs. single copy loci,
# i.e. for miscanthus tags aligned to sorghum reference vs. miscanthus
# reference.
# Can we use this statistic to distinguish good loci from paralogs?
load("workspaces/190515counts_matrices.RData")
tagtab <- read.csv("marker_CSV/190515paralog_tags.csv", stringsAsFactors = FALSE)
Rcpp::sourceCpp('src/simpson.cpp')
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

# sanity check, and filter markers that don't actually have two alignment locations
keep <- logical(nrow(tagtab))
for(i in 1:max(alleles2locS)){
  thiscol <- which(alleles2locS == i)
  if(length(unique(alleles2locM[thiscol])) != 2){
    print(i)
    print(unique(tagtab$Sorghum[thiscol]))
    print(unique(tagtab$Miscanthus[thiscol]))
  } else {
    keep[thiscol] <- TRUE
  }
}

tagtab <- tagtab[keep,]
diploid_mat <- diploid_mat[,keep]
tetraploid_mat <- tetraploid_mat[,keep]

alleles2locM <- match(tagtab$Miscanthus, unique(tagtab$Miscanthus))
alleles2locS <- match(tagtab$Sorghum, unique(tagtab$Sorghum))

# function to get He based on depth ratio instead of counts
He_depth_ratio <- function(alleleDepth, alleles2loc){
  nloc <- max(alleles2loc)
  out <- numeric(nloc)
  for(i in 1:nloc){
    thesecol <- which(alleles2loc == i)
    depthRatios <- sweep(alleleDepth[,thesecol, drop = FALSE],
                         1, rowSums(alleleDepth[,thesecol, drop = FALSE]), "/")
    freq <- colMeans(depthRatios, na.rm = TRUE)
    out[i] <- 1 - sum(freq ^ 2)
  }
  return(out)
}

# get depth ratios
DepthRatio <- function(alleleDepth, alleles2loc){
  locDepth <- t(apply(alleleDepth, 1, function(x) tapply(x, alleles2loc, sum)))
  expandedLocDepth <- locDepth[,as.character(alleles2loc), drop = FALSE]
  depthRatio <- alleleDepth/expandedLocDepth
  return(depthRatio)
}

depthrat_dip_M <- DepthRatio(diploid_mat, alleles2locM)
depthrat_dip_S <- DepthRatio(diploid_mat, alleles2locS)

# compare miscanthus and sorghum for diploids
Depth_dip_M <- tapply(colSums(diploid_mat), alleles2locM, sum)

HindHe_dip_Misc <- HindHeByLoc(diploid_mat, depthrat_dip_M, alleles2locM, max(alleles2locM))
He_dip_Misc <- He_depth_ratio(diploid_mat, alleles2locM)
Hind_dip_Misc <- HindHe_dip_Misc * He_dip_Misc

ggplot(mapping = aes(x = He_dip_Misc, y = Hind_dip_Misc,
                     col = log(Depth_dip_M))) +
  geom_point() +
  scale_color_viridis(option = "magma") +
  geom_abline(slope = 1/2, intercept = 0) +
  labs(x = "He", y = "Hind", title = "Diploids with Miscanthus reference")

hist(log2(HindHe_dip_Misc))
hist(HindHe_dip_Misc[HindHe_dip_Misc < 2], breaks = 50)
mean(HindHe_dip_Misc > 1/2, na.rm = TRUE) # 24% exceed expectations (30% when corrected for size)

Depth_dip_S <- tapply(colSums(diploid_mat), alleles2locS, sum)

HindHe_dip_Sorg <- HindHeByLoc(diploid_mat, depthrat_dip_S, alleles2locS, max(alleles2locS))
He_dip_Sorg <- He_depth_ratio(diploid_mat, alleles2locS)
Hind_dip_Sorg <- HindHe_dip_Sorg * He_dip_Sorg

ggplot(mapping = aes(x = He_dip_Sorg, y = Hind_dip_Sorg,
                     col = log(Depth_dip_S))) +
  geom_point() +
  scale_color_viridis(option = "magma") +
  geom_abline(slope = 1/2, intercept = 0) +
  labs(x = "He", y = "Hind", title = "Diploids with Sorghum reference")

hist(log2(HindHe_dip_Sorg))
hist(HindHe_dip_Sorg[HindHe_dip_Sorg < 2], breaks = 50)
mean(HindHe_dip_Sorg > 1/2, na.rm = TRUE) # 32% exceed expectations (38% corrected for size)

# put into data frame to compare distributions
dip_df <- data.frame(Reference = c(rep("Miscanthus", length(Hind_dip_Misc)),
                                   rep("Sorghum", length(Hind_dip_Sorg))),
                     He = c(He_dip_Misc, He_dip_Sorg),
                     Hind = c(Hind_dip_Misc, Hind_dip_Sorg),
                     HindHe = c(HindHe_dip_Misc, HindHe_dip_Sorg),
                     Depth = c(Depth_dip_M, Depth_dip_S))
dip_df <- dip_df[which(dip_df$Depth > 0 & dip_df$He > 0),]

ggplot(dip_df, aes(x = Depth, y = HindHe, color = Reference)) +
  geom_point() +
#  geom_density2d(color = "black") +
  geom_smooth(color = "black") +
  coord_trans(x = "log") +
  facet_wrap(~ Reference)

ggplot(dip_df, aes(x = log(HindHe/log(Depth)), by = Reference,
                   color = Reference)) +
  geom_density()

ggplot(dip_df, aes(x = log2(HindHe), by = Reference,
                   color = Reference)) +
  geom_density() +
  geom_vline(xintercept = -1)

ggplot(dip_df, aes(x = HindHe, by = Reference,
                   color = Reference)) +
  geom_density() +
  geom_vline(xintercept = 0.5)

# how to properly scale by depth?
# He is prob of sampling two different reads in population

# (Deleted old code for scaling since that is now built into Rcpp function)

# If I have one read He is reduced to zero.
# If I have two reads, and He is 0.5, I get 0.5 * 0 + 0.5 * 0.5 which is 0.25
# If I have three reads, I get 1/4 * 0 + 3/4 * 0.4444 = 1/3

head(dip_df)

ggplot(dip_df, aes(x = log2(HindHe), by = Reference,
                   color = Reference)) +
  geom_density() +
  geom_vline(xintercept = -1)

ggplot(dip_df, aes(x = log(Depth), y = log2(HindHe))) +
  geom_point() +
  geom_density2d(col = "blue") +
  geom_hline(yintercept = -1, color = "red") + # expected max ratio
  geom_vline(xintercept = log(5 * nrow(diploid_mat)), col = "green") + # reasonable minimum depth
  facet_wrap(~ Reference)

ggplot(dip_df[dip_df$Depth >= 5 * nrow(diploid_mat) & dip_df$HindHe > 0,], 
       aes(x = HindHe, by = Reference,
           color = Reference)) +
  geom_density(lwd = 1) +
  geom_vline(xintercept = 0.5, lty = 2) +
  scale_x_continuous(expression(H[ind] / H[E]),
                   c(0.25, 0.5, 0.75, 1), limits = c(0, 1.5))
ggsave("hindhegraph_190904.tiff", compression = "lzw")
