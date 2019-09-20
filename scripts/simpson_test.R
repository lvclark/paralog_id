source('src/hindhe.R') # have to do manually to get file structure right
load("~/Excellence in Breeding/Msa_4x_Chr05_GWAS/Msa_tetraploids_Chr05_RADdata.RData") # tetraploid Msa data
library(polyRAD)
library(ggplot2)
library(viridis)

# alternative dataset; use diploid Msi data
myfile <- system.file("extdata", "Msi01genes.vcf", package = "polyRAD")
myRAD <- VCF2RADdata(myfile, possiblePloidies = list(2))

# test Rcpp estimation of HindHe, which should resemble Python estimation
myHindHe <- colMeans(HindHe(myRAD), na.rm = TRUE)
hist(myHindHe, breaks = 50)
abline(v = 0.75, col = "red")
abline(v = 0.5, col = "blue")
# for tetraploids, peak at 0.63, indicating F of 0.16

# old stuff for old version of functions ####
## compare Simpson index (~=He) of reads from pop vs. reads from individuals

test <- HReadDepth(myRAD$alleleDepth, myRAD$alleles2loc, nLoci(myRAD))
Hpop <- test[[2]]
names(Hpop) <- GetLoci(myRAD)
Hindmat <- test[[1]]
dimnames(Hindmat) <- list(GetTaxa(myRAD), GetLoci(myRAD))
Hind <- colMeans(Hindmat, na.rm = TRUE)

plot(Hpop, Hind)
abline(a = 0, b = 1, col = "red")
abline(a = 0, b = 0.5, col = "blue")
abline(a = 0, b = 2/3, col = "purple")
abline(a = 0, b = 3/4, col = "green")

TotDepth <- colSums(myRAD$locDepth)

ggplot(mapping = aes(x = Hpop, y = Hind, col = log(TotDepth))) +
  geom_point() +
  scale_color_viridis() +
  geom_abline(slope = 0.75, intercept = 0) +
  geom_abline(slope = 0.5, intercept = 0)

ggplot(mapping = aes(y = Hind/Hpop, x = log(TotDepth))) +
  geom_point() +
  geom_density2d() +
  geom_abline(slope = 0, intercept = 0.75, col = "purple")

hist(Hind/Hpop/log(TotDepth), breaks = 50)

hist(Hind/Hpop, breaks = 100, col = "lightgrey")
abline(v = 3/4, col = "purple", lwd = 2)

# What is the expected maximum for Hind/Hpop?
# 0.5 in a biallelic diploid
# max Hind/Hpop = (ploidy - 1)/ploidy 
# --> for autopolyploids, and based on some exploration in my notes

# How does it relate the the number of individuals with more than 4 alleles?
Nmorethan4 <- integer(nLoci(myRAD))
for(L in 1:nLoci(myRAD)){
  theseal <- which(myRAD$alleles2loc == L)
  Nmorethan4[L] <- sum(rowSums(myRAD$alleleDepth[,theseal] > 0) > 4)
}

ggplot(mapping = aes(y = Hind/Hpop, x = log(TotDepth), col = log(Nmorethan4))) +
  geom_point() +
  scale_color_viridis() +
  geom_abline(slope = 0, intercept = 0.75)

hist(Nmorethan4)

plot(Hind/Hpop, log1p(Nmorethan4))
abline(v = 3/4, col = "purple")
## --> Number of individuals w/ more than 4 alleles just not very useful

# investigate some that look paralogish
para <- which(Hind/Hpop > 0.75)
hist(Hpop[para])
Hpop[para[1]]
Hpop[para[2]]

myRAD$alleleDepth[1:50, myRAD$alleles2loc == para[1]]
Hind[para[1]]/Hpop[para[1]]
hist(myRAD$depthRatio[,myRAD$alleles2loc == para[1]])
# --> Just over threshold, doesn't look too bad from depth

myRAD$alleleDepth[1:50, myRAD$alleles2loc == para[2]]
Hind[para[2]]/Hpop[para[2]]
hist(myRAD$depthRatio[,myRAD$alleles2loc == para[2]])
# --> Two alleles that are very common; looks pretty suspicious

myRAD$alleleDepth[1:50, myRAD$alleles2loc == para[3]]
Hind[para[3]]/Hpop[para[3]]
hist(myRAD$depthRatio[,myRAD$alleles2loc == para[3]])
# --> Close to threshold.  A lot of common alleles, depth ratios that don't make sense

# Distribution across taxa
TotDepthT <- rowSums(myRAD$locDepth)
HindT <- rowMeans(Hindmat, na.rm = TRUE)
plot(log(TotDepthT), HindT)

# Look at mix of diploids and tetraploids to see if distinguishable ####
library(VariantAnnotation)
myRAD <- VCF2RADdata("~/DOE Msa study/Seq/GBSv2_180110/180208Msa_filtered.vcf.bgz",
                     svparam = ScanVcfParam(fixed = "ALT", info = NA, geno = "AD",
                                            which = GRanges("05", IRanges(1, 112e6))),
                     yieldSize = NA)

myHindHe <- HindHe(myRAD)

TotDepthT <- rowSums(myRAD$locDepth)

myHindHeByInd <- rowMeans(myHindHe, na.rm = TRUE)

accessions <- read.csv("~/DOE Msa study/Seq/GBSv2_180110/all_accession_names.csv",
                       stringsAsFactors = FALSE)
ploidies <- accessions$Ploidy[match(names(myHindHeByInd), accessions$Accession)]

ggplot(mapping = aes(x = TotDepthT, y = myHindHeByInd, color = ploidies)) +
  geom_point() +
  scale_x_log10()

hist(myHindHeByInd, breaks = 50)
# bimodal distribution, not perfect but could be very helpful

ggplot(data.frame(Depth = TotDepthT, HindHe = myHindHeByInd, Ploidy = ploidies)[ploidies %in% c("2x", "3x", "4x"),],
       mapping = aes(x = Depth, y = HindHe, color = Ploidy)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~ Ploidy)
