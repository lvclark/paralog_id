# Simulate Mendelian and paralogous loci and test methods for distinguishing them

library(polyRAD)
library(ggplot2)

# test distributions for allele frequency
hist(rgamma(1000, shape = 0.3, scale = 1) / 10, breaks = 50)

mean(rgamma(1000, shape = 0.3, scale = 1) / 10 > 0.01) # 48%

# simulate allele frequencies ####
nloci <- 1e3L

shape <- 0.3
scale <- 1

minmaf <- 0.01

set.seed(1229)

# number of alleles per locus
allelesPerLoc <- sample(2:8, nloci, replace = TRUE)

# list for allele frequency output in loop
alFreqList <- vector(mode = "list", length = nloci)

for(i in seq_len(nloci)){
  af <- rgamma(allelesPerLoc[i] - 1L, shape = shape, scale = scale)
  while(sum(af) >= 1 | all(af < minmaf)){
    af <- rgamma(allelesPerLoc[i] - 1L, shape = shape, scale = scale)
  }
  alFreqList[[i]] <- c(af, 1 - sum(af))
}

# generate a second set of allele freqs with no MAF restrictions,
# to represent the paralogous locus.
nloci2 <- nloci %/% 2L
allelesPerLoc2 <- sample(1:8, nloci2, replace = TRUE)
alFreqList2 <- vector(mode = "list", length = nloci2)

for(i in seq_len(nloci2)){
  af <- rgamma(allelesPerLoc2[i] - 1L, shape = shape, scale = scale)
  while(sum(af) >= 1){
    af <- rgamma(allelesPerLoc2[i] - 1L, shape = shape, scale = scale)
  }
  alFreqList2[[i]] <- c(af, 1 - sum(af))
}

# simulate locus depth ####

# parameters and from SimulateRADReads.R; https://doi.org/10.13012/B2IDB-9729830_V2
mndepth_shape = 3.2
mndepth_scale = 8
inddepth_scale = 10

nsam <- 200L

locDepth <- matrix(NA, nrow = nsam, ncol = nloci + nloci2)

for(i in seq_len(nloci + nloci2)){
  # randomly select a mean read depth for this locus
  meanDepth <- ceiling(rgamma(1, shape = mndepth_shape, scale = mndepth_scale))
  inddepth_shape <- meanDepth/inddepth_scale
  # randomly generate total locus depth for each taxon
  indDepth <- as.integer(round(rgamma(nsam, shape = inddepth_shape, scale = inddepth_scale)))
  
  locDepth[,i] <- indDepth
}

# Simulate genotypes ####
inbr <- 0.3 # inbreeding level

alleles2loc1 <- rep(seq_len(nloci), times = allelesPerLoc)
alleles2loc2 <- rep(seq_len(nloci2), times = allelesPerLoc2)

genoDip1 <- SimGenotypes(unlist(alFreqList), alleles2loc1, nsam,
                         inbreeding = inbr, ploidy = 2L)
genoDip2 <- SimGenotypes(unlist(alFreqList2), alleles2loc2, nsam,
                         inbreeding = inbr, ploidy = 2L)

genoTet1 <- SimGenotypes(unlist(alFreqList), alleles2loc1, nsam,
                         inbreeding = inbr, ploidy = 4L)
genoTet2 <- SimGenotypes(unlist(alFreqList2), alleles2loc2, nsam,
                         inbreeding = inbr, ploidy = 4L)

# Simulate allelic read depth ####
locDepth1 <- locDepth[,1:nloci]
locDepth2 <- locDepth[,(1:nloci2) + nloci]
colnames(locDepth1) <- as.character(seq_len(nloci))
colnames(locDepth2) <- as.character(seq_len(nloci2))

ADdip1 <- SimAlleleDepth(locDepth1, genoDip1, alleles2loc1,
                         overdispersion = 20, contamRate = 0.001)
ADdip2 <- SimAlleleDepth(locDepth2, genoDip2, alleles2loc2,
                         overdispersion = 20, contamRate = 0.001)

ADtet1 <- SimAlleleDepth(locDepth1, genoTet1, alleles2loc1,
                         overdispersion = 20, contamRate = 0.001)
ADtet2 <- SimAlleleDepth(locDepth2, genoTet2, alleles2loc2,
                         overdispersion = 20, contamRate = 0.001)

rownames(ADdip1) <- rownames(ADdip2) <- rownames(ADtet1) <- rownames(ADtet2) <-
  paste0("sam", seq_len(nsam))

colnames(ADdip1) <- colnames(ADtet1) <-
  paste0("loc", alleles2loc1, "_", unlist(lapply(allelesPerLoc, seq_len)))

colnames(ADdip2) <- colnames(ADtet2) <-
  paste0("loc", alleles2loc2, "_", unlist(lapply(allelesPerLoc2, seq_len)) +
           allelesPerLoc[alleles2loc2])

# Build RADdata objects and call genotypes ####
# random allele nucleotides
alnuc <- do.call(paste0,
                 lapply(1:10, function(x) sample(c("A", "C", "G", "T"), ncol(ADdip1) + ncol(ADdip2), replace = TRUE)))

RADdip <- RADdata(cbind(ADdip1, ADdip2), c(alleles2loc1, alleles2loc2),
                  data.frame(row.names = paste("loc", seq_len(nloci))),
                  list(2L), 0.001, alnuc)

RADtet <- RADdata(cbind(ADtet1, ADtet2), c(alleles2loc1, alleles2loc2),
                  data.frame(row.names = paste("loc", seq_len(nloci))),
                  list(4L), 0.001, alnuc)

RADdip <- IterateHWE(RADdip)
RADtet <- IterateHWE(RADtet)

# Run statistics on loci ####
source("scripts/HoHe.R")

hhDip <- colMeans(HindHe(RADdip), na.rm = TRUE)
hhTet <- colMeans(HindHe(RADtet), na.rm = TRUE)

oeDip <- colMeans(HoHe(GetProbableGenotypes(RADdip, omit1allelePerLocus = FALSE)[[1]],
                       RADdip$alleles2loc, 2L))

hapDip <- colMeans(HapPerGen(RADdip$alleleDepth, RADdip$alleles2loc) > 2L)
hapTet <- colMeans(HapPerGen(RADtet$alleleDepth, RADtet$alleles2loc) > 4L)

ggplot(mapping = aes(x = hhDip,
                     fill = rep(c("Paralog", "Mendelian"),
                                 times = c(nloci2, nloci - nloci2)))) +
  geom_density(alpha = 0.5) +
  labs(x = "Hind/He", fill = "Locus type")

ggplot(mapping = aes(x = oeDip,
                     fill = rep(c("Paralog", "Mendelian"),
                                times = c(nloci2, nloci - nloci2)))) +
  geom_density(alpha = 0.5) +
  labs(x = "Ho/He", fill = "Locus type")

ggplot(mapping = aes(x = hapDip,
                     fill = rep(c("Paralog", "Mendelian"),
                                times = c(nloci2, nloci - nloci2)))) +
  geom_density(alpha = 0.5) +
  labs(x = "Proportion samples with more haplotypes than expected", fill = "Locus type") +
  scale_x_continuous(trans = "log1p")

ggplot(mapping = aes(x = hapTet,
                     fill = rep(c("Paralog", "Mendelian"),
                                times = c(nloci2, nloci - nloci2)))) +
  geom_density(alpha = 0.5) +
  labs(x = "Proportion samples with more haplotypes than expected", fill = "Locus type") +
  scale_x_continuous(trans = "log1p")
