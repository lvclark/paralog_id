# examine H statistics for double copy vs. single copy loci,
# i.e. for miscanthus tags aligned to sorghum reference vs. miscanthus
# reference.
# Can we use this statistic to distinguish good loci from paralogs?
load("workspaces/190515counts_matrices.RData")
tagtab <- read.csv("marker_CSV/190515paralog_tags.csv", stringsAsFactors = FALSE)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(viridis)
library(dplyr)
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

# compare for tetraploids
Depth_tet_M <- tapply(colSums(tetraploid_mat), alleles2locM, sum)

HindHe_tet_Misc <- colMeans(hhM4, na.rm = TRUE)

Depth_tet_S <- tapply(colSums(tetraploid_mat), alleles2locS, sum)

HindHe_tet_Sorg <- colMeans(hhS4, na.rm = TRUE)

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

p2 <- ggplot(dip_df[dip_df$Depth >= 5 * nrow(diploid_mat) & dip_df$HindHe > 0,], 
       aes(x = HindHe, by = Reference,
           color = Reference)) +
  geom_density(lwd = 1) +
  geom_vline(xintercept = 0.5, lty = 2) +
  scale_x_continuous(expression(H[ind] / H[E]),
                   c(0.25, 0.5, 0.75, 1), limits = c(0, 1.5))
#ggsave("hindhegraph_190904.tiff", compression = "lzw")

# add tetraploid data frame

tet_df <- data.frame(Reference = c(rep("Miscanthus", length(HindHe_tet_Misc)),
                                   rep("Sorghum", length(HindHe_tet_Sorg))),
                     HindHe = c(HindHe_tet_Misc, HindHe_tet_Sorg),
                     Depth = c(Depth_tet_M, Depth_tet_S))
tet_df <- tet_df[which(tet_df$Depth > 0 & !is.na(tet_df$HindHe)),]

p4 <- ggplot(tet_df[tet_df$Depth >= 5 * nrow(tetraploid_mat) & tet_df$HindHe > 0,], 
       aes(x = HindHe, by = Reference,
           color = Reference)) +
  geom_density(lwd = 1) +
  geom_vline(xintercept = 0.75, lty = 2) +
  scale_x_continuous(expression(H[ind] / H[E]),
                     c(0.25, 0.5, 0.75, 1), limits = c(0, 1.5))

# plot together ####
all_df <- data.frame(Reference = c(rep("Miscanthus", length(HindHe_dip_Misc)),
                                   rep("Sorghum", length(HindHe_dip_Sorg))),
                     HindHe_diploids = c(HindHe_dip_Misc, HindHe_dip_Sorg),
                     HindHe_tetraploids = c(HindHe_tet_Misc, HindHe_tet_Sorg),
                     Depth = c(Depth_dip_M + Depth_tet_M, Depth_dip_S + Depth_tet_S))
all_df <- all_df[which(all_df$Depth > 0 & !is.na(all_df$HindHe_diploids) &
                         !is.na(all_df$HindHe_tetraploids)),]

sum(all_df$Reference == "Miscanthus") # 21590
sum(all_df$Reference == "Sorghum")    # 17860

#save(all_df, file = "workspaces/all_df.RData")

all_df_filt <- filter(all_df, #HindHe_diploids < 2, HindHe_tetraploids < 2,
                      Depth >= 5 * (nrow(diploid_mat) + nrow(tetraploid_mat)))

sum(all_df_filt$Reference == "Miscanthus") # 11516
sum(all_df_filt$Reference == "Sorghum")    #  8820

p24 <- ggplot(all_df_filt, aes(x = HindHe_diploids, y = HindHe_tetraploids,
           color = Depth / (nrow(diploid_mat) + nrow(tetraploid_mat)))) +
  geom_point() +
  geom_density_2d(color = "black") +
  scale_color_viridis(trans = "log", breaks = c(7, 20, 50, 150, 400)) +
  geom_vline(xintercept = 0.5, lty = 2) +
  geom_hline(yintercept = 0.75, lty = 2) +
  facet_wrap(~ Reference) +
  labs(x = expression(H[ind] / H[E] ~ ", diploids"),
       y = expression(H[ind] / H[E] ~ ", tetraploids"),
       color = "Mean depth")

#tiff("191218miscanthus_v_sorghum.tiff", width = 6.5 * 300, height = 5 * 300,
#     res = 300, compression = "lzw")
grid.arrange(arrangeGrob(p2 + theme(legend.position="none") + ggtitle("Diploids"), 
                         p4 + ggtitle("Tetraploids"),
                         p24, layout_matrix = matrix(c(1,3,2,3), nrow = 2, ncol = 2),
                         widths = c(0.75, 1)))
#dev.off()
