# determine inbreeding after running process_sam_multi

library(polyRAD)
library(ggplot2)
library(viridis)

Msa1 <- readProcessSamMulti("large_datasets/Msa_split_1_align.csv")

diploids <- readLines("marker_csv/diploids.txt")
tetraploids <- readLines("marker_csv/tetraploids.txt")

Msa1_2x <- SubsetByTaxon(Msa1, diploids)
Msa1_4x <- SubsetByTaxon(Msa1, tetraploids)

hh2 <- HindHe(Msa1_2x)
hh4 <- HindHe(Msa1_4x)

hist(colMeans(hh2, na.rm = TRUE), breaks = 30)
hist(colMeans(hh4, na.rm = TRUE), breaks = 30)

InbreedingFromHindHe(0.4, 2) # 0.2
InbreedingFromHindHe(0.65, 4) # 0.13

# use 0.15 as rough estimate for both

# distinguishing ploidy ####
set.seed(715)
Msa1 <- readProcessSamMulti("large_datasets/Msa_split_1_align.csv",
                            expectedLoci = 1e4)
hhAll <- HindHe(Msa1)
hhInd <- rowMeans(hhAll, na.rm = TRUE)
hist(hhInd)
hist(hhInd[diploids])
hist(hhInd[tetraploids])

nTaxa(Msa1)
length(diploids) + length(tetraploids)

# get a table with both ploidy and proportion ancestry M. sinensis for plots
accessions <- read.csv("marker_csv/all_accession_names.csv", stringsAsFactors = FALSE)
all(names(hhInd) %in% accessions$Accession)

accessions <- accessions[match(names(hhInd), accessions$Accession),]
identical(names(hhInd), accessions$Accession)

accessions$HindHe <- hhInd
accessions <- accessions[accessions$Ploidy %in% c("2x", "3x", "4x"),]
table(accessions$Ploidy)
summary(accessions$Prop_Msi)
hist(accessions$Prop_Msi)
accessions[which(accessions$Prop_Msi > 0.5),]

accessions <- accessions[which(accessions$Prop_Msi < 0.5),]

accessions$Depth <- rowMeans(Msa1$locDepth[accessions$Accession,])

# plot Hind/He vs. ploidy and proportion Msi ancestry
# tiff("191219hindhe_by_ploidy.tiff", width = 6.5 * 300, height = 4 * 300, res = 300,
#      compression = "lzw")
#pdf("Fig4_HindHe_by_ind.pdf", width = 6.7, height = 4)
ggplot(accessions, aes(x = Depth, y = HindHe, color = Prop_Msi)) +
  geom_point() +
  facet_wrap(~ Ploidy) +
  scale_color_viridis() + 
  labs(x = "Mean read depth", y = expression(H[ind] / H[E]),
       color = "Hybrid ancestry") +
  #coord_trans(x = "log")
  scale_x_continuous(trans = "log2", breaks = c(5, 10, 20, 40, 80, 160, 320)) +
  geom_hline(data = data.frame(Ploidy = c("2x", "3x", "4x"), ExpectVal = c(0.5, 2/3, 0.75)),
             mapping = aes(yintercept = ExpectVal), lty = 2)
#dev.off()
