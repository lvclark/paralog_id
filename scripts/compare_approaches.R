# Simulate Mendelian and paralogous loci and test methods for distinguishing them

library(polyRAD)

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
  while(sum(af >= 1) | all(af < minmaf)){
    af <- rgamma(allelesPerLoc[i] - 1L, shape = shape, scale = scale)
  }
  alFreqList[[i]] <- c(af, 1 - sum(af))
}

# generate a second set of allele freqs with no MAF restrictions,
# to represent the paralogous locus.
nloci2 <- nloci %/% 2L
allelesPerLoc2 <- sample(1:8, nloci2)
alFreqList2 <- vector(mode = "list", length = nloci2)

for(i in seq_len(nloci2)){
  af <- rgamma(allelesPerLoc2[i] - 1L, shape = shape, scale = scale)
  while(sum(af >= 1)){
    af <- rgamma(allelesPerLoc2[i] - 1L, shape = shape, scale = scale)
  }
  alFreqList2[[i]] <- c(af, 1 - sum(af))
}
