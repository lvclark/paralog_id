# examine H statistics for double copy vs. single copy loci,
# i.e. for miscanthus tags aligned to sorghum reference vs. miscanthus
# reference.
# Can we use this statistic to distinguish good loci from paralogs?
load("workspaces/190515counts_matrices.RData")
tagtab <- read.csv("marker_CSV/190515paralog_tags.csv", stringsAsFactors = FALSE)
library(ggplot2)
library(viridis)
library(polyRAD)

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

# RADdata objects and Hind/He
radM2 <- RADdata(diploid_mat, alleles2locM,
                 data.frame(row.names = unique(tagtab$Miscanthus)),
                 list(2), 0.001, tagtab$TagSeq)
hhM2 <- HindHe(radM2)
rm(radM2)

radS2 <- RADdata(diploid_mat, alleles2locS,
                 data.frame(row.names = unique(tagtab$Sorghum)),
                 list(2), 0.001, tagtab$TagSeq)
hhS2 <- HindHe(radS2)
rm(radS2)

radM4 <- RADdata(tetraploid_mat, alleles2locM,
                 data.frame(row.names = unique(tagtab$Miscanthus)),
                 list(4), 0.001, tagtab$TagSeq)
hhM4 <- HindHe(radM4)
rm(radM4)

radS4 <- RADdata(tetraploid_mat, alleles2locS,
                 data.frame(row.names = unique(tagtab$Sorghum)),
                 list(4), 0.001, tagtab$TagSeq)
hhS4 <- HindHe(radS4)
rm(radS4)

# compare miscanthus and sorghum for diploids
Depth_dip_M <- tapply(colSums(diploid_mat), alleles2locM, sum)

HindHe_dip_Misc <- colMeans(hhM2, na.rm = TRUE)

hist(HindHe_dip_Misc[HindHe_dip_Misc < 2], breaks = 50)

mean(HindHe_dip_Misc > 1/2, na.rm = TRUE) # 30% exceed expectations

Depth_dip_S <- tapply(colSums(diploid_mat), alleles2locS, sum)

HindHe_dip_Sorg <- colMeans(hhS2, na.rm = TRUE)

hist(HindHe_dip_Sorg[HindHe_dip_Sorg < 2], breaks = 50)
mean(HindHe_dip_Sorg > 1/2, na.rm = TRUE) # 39% exceed expectations

# put into data frame to compare distributions
dip_df <- data.frame(Reference = c(rep("Miscanthus", length(HindHe_dip_Misc)),
                                   rep("Sorghum", length(HindHe_dip_Sorg))),
                     HindHe = c(HindHe_dip_Misc, HindHe_dip_Sorg),
                     Depth = c(Depth_dip_M, Depth_dip_S))
dip_df <- dip_df[which(dip_df$Depth > 0 & !is.na(dip_df$HindHe)),]

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
