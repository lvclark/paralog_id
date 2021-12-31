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
# to represent paralogous loci.
nloci2 <- nloci * 2L
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
inbr <- 0.1 # inbreeding level

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

# combine pairs of adjacent loci into single loci
colnames(ADdip2) <- colnames(ADtet2) <-
  paste0("loc", (alleles2loc2 + 1L) %/% 2L + nloci, "_",
         unlist(lapply(allelesPerLoc2[seq(1, nloci2 - 1, by = 2)] +
                         allelesPerLoc2[seq(2, nloci2, by = 2)], seq_len)))

# Build RADdata objects and call genotypes ####
# random allele nucleotides
alnuc <- do.call(paste0,
                 lapply(1:10, function(x) sample(c("A", "C", "G", "T"), ncol(ADdip1) + ncol(ADdip2), replace = TRUE)))

alleles2loc_combined <- c(alleles2loc1, (alleles2loc2 + 1L) %/% 2L + nloci)

RADdip <- RADdata(cbind(ADdip1, ADdip2), alleles2loc_combined,
                  data.frame(row.names = paste0("loc", seq_len(nloci * 2))),
                  list(2L), 0.001, alnuc)

RADtet <- RADdata(cbind(ADtet1, ADtet2), alleles2loc_combined,
                  data.frame(row.names = paste0("loc", seq_len(nloci * 2))),
                  list(4L), 0.001, alnuc)

RADdip <- IterateHWE(RADdip)
RADtet <- IterateHWE(RADtet)

# consider also simulating with null alleles to make the data messier

# Run statistics on loci ####
source("scripts/HoHe.R")

hhDip <- colMeans(HindHe(RADdip), na.rm = TRUE)
hhTet <- colMeans(HindHe(RADtet), na.rm = TRUE)

genoDip <- GetProbableGenotypes(RADdip, omit1allelePerLocus = FALSE, multiallelic = "na")[[1]]
genoTet <- GetProbableGenotypes(RADtet, omit1allelePerLocus = FALSE, multiallelic = "na")[[1]]

mean(is.na(genoDip[,RADdip$alleles2loc > nloci]))  # 54% of paralogs don't have allele copy num adding up
mean(is.na(genoDip[,RADdip$alleles2loc <= nloci])) # only 3% of non-paralogs have that issue
mean(is.na(genoTet[,RADtet$alleles2loc > nloci]))  # 55%
mean(is.na(genoTet[,RADtet$alleles2loc <= nloci])) # 26%

oeDip <- colMeans(HoHe(genoDip, RADdip$alleles2loc, 2L), na.rm = TRUE)
oeTet <- colMeans(HoHe(genoTet, RADtet$alleles2loc, 4L), na.rm = TRUE)

hapDip <- colMeans(HapPerGen(RADdip$alleleDepth, RADdip$alleles2loc) > 2L)
hapTet <- colMeans(HapPerGen(RADtet$alleleDepth, RADtet$alleles2loc) > 4L)

depthDip <- colMeans(GetLocDepth(RADdip))
depthTet <- colMeans(GetLocDepth(RADtet))

Zdip <- abs(Zscore(RADdip$alleleDepth, genoDip, RADdip$alleles2loc, 2L))
Ztet <- abs(Zscore(RADtet$alleleDepth, genoTet, RADtet$alleles2loc, 4L))

ggplot(mapping = aes(x = hhDip,
                     fill = rep(c("Mendelian", "Paralog"),
                                 each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Hind/He", fill = "Locus type")

ggplot(mapping = aes(x = oeDip,
                     fill = rep(c("Mendelian", "Paralog"),
                                 each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Ho/He", fill = "Locus type")

ggplot(mapping = aes(x = hapDip,
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Proportion samples with more haplotypes than expected", fill = "Locus type") +
  scale_x_continuous(trans = "log1p")

ggplot(mapping = aes(x = hapTet,
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Proportion samples with more haplotypes than expected", fill = "Locus type") +
  scale_x_continuous(trans = "log1p")

ggplot(mapping = aes(x = Zdip,
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Z score", fill = "Locus type")

ggplot(mapping = aes(x = Ztet,
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Z score", fill = "Locus type")

# Summarize approach efficiency ####
# Get 95th percentile for Mendelian markers. How many paralogs would be filtered at that cutoff?
quantile(hhDip[1:nloci], 0.95)

summtab <- data.frame(Approach = rep(c("Hind/He", "Ho/He", "Haplotypes", "Depth", "Z-score"), each = 2),
                      Ploidy = rep(c("Diploid", "Tetraploid"), times = 5),
                      Mendelian95 = NA_real_,
                      ParalogsFiltered = NA_real_,
                      SE = NA_real_)

for(i in seq_len(nrow(summtab))){
  x <- list(hhDip, hhTet, oeDip, oeTet, hapDip, hapTet, depthDip, depthTet, Zdip, Ztet)[[i]]
  q <- quantile(x[1:nloci], 0.95,  na.rm = TRUE)
  p <- mean(x[(1:nloci) + nloci] > q, na.rm = TRUE)
  summtab$Mendelian95[i] <- q
  summtab$ParalogsFiltered[i] <- p
  summtab$SE[i] <- sqrt(p * (1 - p) / sum(!is.na(x[(1:nloci) + nloci])))
}

summtab
