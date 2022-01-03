# Simulate Mendelian and paralogous loci and test methods for distinguishing them

library(polyRAD)
library(ggplot2)
library(tictoc)

tic.clear()
tic.clearlog() # start fresh recording timings

# function to extract timings from log
extractTimings <- function(log){
  out <- sapply(log, function(x) x$toc - x$tic)
  names(out) <- sapply(log, function(x) x$msg)
  return(out)
}

# object to store timings
tictocres <-
  array(list(), dim = c(7,2),
        dimnames = list(c("IterateHWE", "HindHe", "GetProbableGenotypes",
                          "HoHe", "Haplotypes", "Depth", "Zscore"), c("Dip", "Tet")))

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
  af <- rgamma(allelesPerLoc[i] - 1L, shape = shape, scale = scale) / 10 + minmaf
  while(sum(af) >= 1){
    af <- rgamma(allelesPerLoc[i] - 1L, shape = shape, scale = scale) / 10 + minmaf
  }
  alFreqList[[i]] <- c(af, 1 - sum(af))
}

# generate a second set of allele freqs with no MAF restrictions,
# to represent paralogous loci.
nloci2 <- nloci * 2L
allelesPerLoc2 <- sample(1:8, nloci2, replace = TRUE)
alFreqList2 <- vector(mode = "list", length = nloci2)

for(i in seq_len(nloci2)){
  af <- rgamma(allelesPerLoc2[i] - 1L, shape = shape, scale = scale) / 10 + minmaf
  while(sum(af) >= 1){
    af <- rgamma(allelesPerLoc2[i] - 1L, shape = shape, scale = scale) / 10 + minmaf
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
inbr <- seq(0, 1, by = 0.1) # inbreeding level

alleles2loc1 <- rep(seq_len(nloci), times = allelesPerLoc)
alleles2loc2 <- rep(seq_len(nloci2), times = allelesPerLoc2)

genoDip1 <- lapply(inbr,
                   function(x){
                     SimGenotypes(unlist(alFreqList), alleles2loc1, nsam,
                         inbreeding = x, ploidy = 2L)
                   })
genoDip2 <- lapply(inbr,
                   function(x){
                     SimGenotypes(unlist(alFreqList2), alleles2loc2, nsam,
                         inbreeding = x, ploidy = 2L)
                   })

genoTet1 <- lapply(inbr,
                   function(x){
                     SimGenotypes(unlist(alFreqList), alleles2loc1, nsam,
                         inbreeding = x, ploidy = 4L)
                   })
genoTet2 <- lapply(inbr,
                   function(x){
                     SimGenotypes(unlist(alFreqList2), alleles2loc2, nsam,
                         inbreeding = x, ploidy = 4L)
                   })

# Simulate allelic read depth ####
locDepth1 <- locDepth[,1:nloci]
locDepth2 <- locDepth[,(1:nloci2) + nloci]
colnames(locDepth1) <- as.character(seq_len(nloci))
colnames(locDepth2) <- as.character(seq_len(nloci2))

ADdip1 <- lapply(genoDip1,
                 function(x){
                   out <- SimAlleleDepth(locDepth1, x, alleles2loc1,
                         overdispersion = 20, contamRate = 0.001)
                   rownames(out) <- paste0("sam", seq_len(nsam))
                   colnames(out) <- paste0("loc", alleles2loc1, "_", unlist(lapply(allelesPerLoc, seq_len)))
                   out
                 })
ADdip2 <- lapply(genoDip2,
                 function(x){
                   out <- SimAlleleDepth(locDepth2, x, alleles2loc2,
                         overdispersion = 20, contamRate = 0.001)
                   rownames(out) <- paste0("sam", seq_len(nsam))
                   # combine pairs of adjacent loci into single loci
                   colnames(out) <- paste0("loc", (alleles2loc2 + 1L) %/% 2L + nloci, "_",
                                           unlist(lapply(allelesPerLoc2[seq(1, nloci2 - 1, by = 2)] +
                                                           allelesPerLoc2[seq(2, nloci2, by = 2)], seq_len)))
                   out
                 })

ADtet1 <- lapply(genoTet1,
                 function(x){
                   out <- SimAlleleDepth(locDepth1, x, alleles2loc1,
                         overdispersion = 20, contamRate = 0.001)
                   rownames(out) <- paste0("sam", seq_len(nsam))
                   colnames(out) <- paste0("loc", alleles2loc1, "_", unlist(lapply(allelesPerLoc, seq_len)))
                   out
                 })
ADtet2 <- lapply(genoTet2,
                 function(x){
                   out <- SimAlleleDepth(locDepth2, x, alleles2loc2,
                         overdispersion = 20, contamRate = 0.001)
                   rownames(out) <- paste0("sam", seq_len(nsam))
                   # combine pairs of adjacent loci into single loci
                   colnames(out) <- paste0("loc", (alleles2loc2 + 1L) %/% 2L + nloci, "_",
                                           unlist(lapply(allelesPerLoc2[seq(1, nloci2 - 1, by = 2)] +
                                                           allelesPerLoc2[seq(2, nloci2, by = 2)], seq_len)))
                   out
                 })

# Build RADdata objects and call genotypes ####
# random allele nucleotides
alnuc <- do.call(paste0,
                 lapply(1:10, function(x) sample(c("A", "C", "G", "T"), ncol(ADdip1[[1]]) + ncol(ADdip2[[1]]), replace = TRUE)))

alleles2loc_combined <- c(alleles2loc1, (alleles2loc2 + 1L) %/% 2L + nloci)

tic.clearlog()

RADdip <- lapply(1:11,
                 function(x){
                   out <- RADdata(cbind(ADdip1[[x]], ADdip2[[x]]), alleles2loc_combined,
                  data.frame(row.names = paste0("loc", seq_len(nloci * 2))),
                  list(2L), 0.001, alnuc)
                   tic(paste("IterateHWE in diploids at inbreeding", inbr[x]))
                   out <- IterateHWE(out)
                   toc(log = TRUE)
                   return(out)
                 })

tictocres[["IterateHWE", "Dip"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()

RADtet <- lapply(1:11,
                 function(x){
                   out <- RADdata(cbind(ADtet1[[x]], ADtet2[[x]]), alleles2loc_combined,
                  data.frame(row.names = paste0("loc", seq_len(nloci * 2))),
                  list(4L), 0.001, alnuc)
                   tic(paste("IterateHWE in tetraploids at inbreeding", inbr[x]))
                   out <- IterateHWE(out)
                   toc(log = TRUE)
                   return(out)
                 })

tictocres[["IterateHWE", "Tet"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()

# consider also simulating with null alleles to make the data messier

# Run statistics on loci ####
source("scripts/HoHe.R")

hhDip <- lapply(1:11, function(x){
  tic(paste("HindHe in diploids at inbreeding", inbr[x]))
  out <- colMeans(HindHe(RADdip[[x]]), na.rm = TRUE)
  toc(log = TRUE)
  return(out)
} )
tictocres[["HindHe", "Dip"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()
hhTet <- lapply(1:11, function(x){
  tic(paste("HindHe in tetraploids at inbreeding", inbr[x]))
  out <- colMeans(HindHe(RADtet[[x]]), na.rm = TRUE)
  toc(log = TRUE)
  return(out)
} )
tictocres[["HindHe", "Tet"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()

genoDip <- lapply(1:11,
                  function(x){
                    tic(paste("GetProbableGenotypes in diploids at inbreeding", inbr[x]))
                    out <- GetProbableGenotypes(RADdip[[x]], omit1allelePerLocus = FALSE, multiallelic = "na")[[1]]
                    toc(log = TRUE)
                    return(out)
                  })
tictocres[["GetProbableGenotypes", "Dip"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()
genoTet <- lapply(1:11,
                  function(x){
                    tic(paste("GetProbableGenotypes in tetraploids at inbreeding", inbr[x]))
                    out <- GetProbableGenotypes(RADtet[[x]], omit1allelePerLocus = FALSE, multiallelic = "na")[[1]]
                    toc(log = TRUE)
                    return(out)
                  })
tictocres[["GetProbableGenotypes", "Tet"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()

mean(is.na(genoDip[[2]][,RADdip$alleles2loc > nloci]))  # 36% of paralogs don't have allele copy num adding up
mean(is.na(genoDip[[2]][,RADdip$alleles2loc <= nloci])) # only 1% of non-paralogs have that issue
mean(is.na(genoTet[[2]][,RADtet$alleles2loc > nloci]))  # 40%
mean(is.na(genoTet[[2]][,RADtet$alleles2loc <= nloci])) # 10%

oeDip <- lapply(1:11,
                function(x){
                  tic(paste("HoHe in diploids at inbreeding", inbr[x]))
                  out <- colMeans(HoHe(genoDip[[x]], RADdip[[x]]$alleles2loc, 2L), na.rm = TRUE)
                  toc(log = TRUE)
                  return(out)
                })
tictocres[["HoHe", "Dip"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()
oeTet <- lapply(1:11,
                function(x){
                  tic(paste("HoHe in tetraploids at inbreeding", inbr[x]))
                  out <- colMeans(HoHe(genoTet[[x]], RADtet[[x]]$alleles2loc, 4L), na.rm = TRUE)
                  toc(log = TRUE)
                  return(out)
                })
tictocres[["HoHe", "Tet"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()

hapDip <- lapply(1:11,
                 function(x){
                   tic(paste("Haplotype count in diploids at inbreeding", inbr[x]))
                   out <- colMeans(HapPerGen(RADdip[[x]]$alleleDepth, RADdip[[x]]$alleles2loc) > 2L)
                   toc(log = TRUE)
                   return(out)
                 })
tictocres[["Haplotypes", "Dip"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()
hapTet <- lapply(1:11,
                 function(x){
                   tic(paste("Haplotype count in tetraploids at inbreeding", inbr[x]))
                   out <- colMeans(HapPerGen(RADtet[[x]]$alleleDepth, RADtet[[x]]$alleles2loc) > 4L)
                   toc(log = TRUE)
                   return(out)
                 })
tictocres[["Haplotypes", "Tet"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()

tic("Mean locus depth")
depthDip <- colMeans(GetLocDepth(RADdip[[1]]))
toc(log = TRUE)
tictocres[["Depth", "Dip"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()
#depthTet <- colMeans(GetLocDepth(RADtet))

Zdip <- lapply(1:11,
               function(x){
                 tic(paste("Z-score in diploids at inbreeding", inbr[x]))
                 out <- abs(Zscore(RADdip[[x]]$alleleDepth, genoDip[[x]], RADdip[[x]]$alleles2loc, 2L))
                 toc(log = TRUE)
                 return(out)
               })
tictocres[["Zscore", "Dip"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()
Ztet <- lapply(1:11,
               function(x){
                 tic(paste("Z-score in tetraploids at inbreeding", inbr[x]))
                 out <- abs(Zscore(RADtet[[x]]$alleleDepth, genoTet[[x]], RADtet[[x]]$alleles2loc, 4L))
                 toc(log = TRUE)
                 return(out)
               })
tictocres[["Zscore", "Tet"]] <- extractTimings(tic.log(format = FALSE))
tic.clearlog()

ggplot(mapping = aes(x = hhDip[[2]],
                     fill = rep(c("Mendelian", "Paralog"),
                                 each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Hind/He", fill = "Locus type")

ggplot(mapping = aes(x = oeDip[[2]],
                     fill = rep(c("Mendelian", "Paralog"),
                                 each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Ho/He", fill = "Locus type")

ggplot(mapping = aes(x = hapDip[[2]],
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Proportion samples with more haplotypes than expected", fill = "Locus type") +
  scale_x_continuous(trans = "log1p")

ggplot(mapping = aes(x = hapTet[[2]],
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Proportion samples with more haplotypes than expected", fill = "Locus type") +
  scale_x_continuous(trans = "log1p")

ggplot(mapping = aes(x = Zdip[[2]],
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Z score", fill = "Locus type")

ggplot(mapping = aes(x = Ztet[[2]],
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Z score", fill = "Locus type")

# Summarize approach efficiency ####
# Get 95th percentile for Mendelian markers. How many paralogs would be filtered at that cutoff?
quantile(hhDip[1:nloci], 0.95)

summtab <- data.frame(Approach = rep(c("Hind/He", "Ho/He", "Haplotypes", "Z-score", "Depth"), each = 22),
                      Ploidy = rep(rep(c("Diploid", "Tetraploid"), each = 11), times = 5),
                      Inbreeding = rep(inbr, times = 10),
                      Mendelian95 = NA_real_,
                      ParalogsFiltered = NA_real_,
                      SE = NA_real_,
                      ComputationTime = c(tictocres[["HindHe", "Dip"]],
                                          tictocres[["HindHe", "Tet"]],
                                          tictocres[["HoHe", "Dip"]] + tictocres[["IterateHWE", "Dip"]] + tictocres[["GetProbableGenotypes", "Dip"]],
                                          tictocres[["HoHe", "Tet"]] + tictocres[["IterateHWE", "Tet"]] + tictocres[["GetProbableGenotypes", "Tet"]],
                                          tictocres[["Haplotypes", "Dip"]],
                                          tictocres[["Haplotypes", "Tet"]],
                                          tictocres[["Zscore", "Dip"]] + tictocres[["IterateHWE", "Dip"]] + tictocres[["GetProbableGenotypes", "Dip"]],
                                          tictocres[["Zscore", "Tet"]] + tictocres[["IterateHWE", "Tet"]] + tictocres[["GetProbableGenotypes", "Tet"]],
                                          rep(tictocres[["Depth", "Dip"]], 22)))

summtab$Approach <- factor(summtab$Approach, levels = unique(summtab$Approach))

for(i in seq_len(nrow(summtab))){
  testindex <- (i + 10L) %/% 11L
  if(testindex < 9){
    inbreedingindex <- (i + 10L) %% 11L + 1L
    x <- list(hhDip, hhTet, oeDip, oeTet, hapDip, hapTet, Zdip, Ztet)[[testindex]][[inbreedingindex]]
  } else {
    x <- depthDip
  }
  
  q <- quantile(x[1:nloci], 0.95,  na.rm = TRUE)
  p <- mean(x[(1:nloci) + nloci] > q, na.rm = TRUE)
  summtab$Mendelian95[i] <- q
  summtab$ParalogsFiltered[i] <- p
  summtab$SE[i] <- sqrt(p * (1 - p) / sum(!is.na(x[(1:nloci) + nloci])))
}

summtab

median(unlist(alFreqList))

median(locDepth)

# Summary figure comparing methods
ggplot(summtab, aes(x = Inbreeding, y = ParalogsFiltered, color = Approach, linetype = Ploidy)) +
  geom_ribbon(aes(ymin = ParalogsFiltered - SE,
                  ymax = ParalogsFiltered + SE,
                  group = paste(Approach, Ploidy)),
              fill = "lightgrey") +
  geom_point() +
  geom_line() +
  scale_x_continuous(n.breaks = 11, minor_breaks = NULL) +
  scale_color_brewer(palette = "Dark2")

# Why does Z-score dip down and back up for tetraploids?
sapply(genoTet, function(x) sum(x == 3L, na.rm = TRUE))
sapply(genoTet, function(x) sum(x[,1:ncol(genoTet1[[1]])] == 3L, na.rm = TRUE))
sapply(Ztet, function(x) sum(is.na(x)))

ggplot(mapping = aes(x = Ztet[[1]],
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Z score", fill = "Locus type")

ggplot(mapping = aes(x = Ztet[[6]],
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Z score", fill = "Locus type")

ggplot(mapping = aes(x = Ztet[[11]],
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Z score", fill = "Locus type")

# Increasing proportion of heterozygous genotypes that are due to error
# and thus have high Z-scores.  Then paralogs catch up.

# Why does Ho/He decrease in efficiency at high inbreeding?
# 95th percentile for Mendelian stays the same
sapply(oeTet, function(x) median(x[1:nloci], na.rm = TRUE))
# median goes down

ggplot(mapping = aes(x = oeTet[[11]],
                     fill = rep(c("Mendelian", "Paralog"),
                                each = nloci))) +
  geom_density(alpha = 0.5) +
  labs(x = "Ho/He", fill = "Locus type")
