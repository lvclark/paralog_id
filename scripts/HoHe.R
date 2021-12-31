# Functions for other paralog detection approaches

# Observed over Expected heterozygosity
HoHe <- function(genmat, alleles2loc, ploidy){
  alleleFreq <- colMeans(genmat, na.rm = TRUE) / ploidy
  nloc <- max(alleles2loc)
  nsam <- nrow(genmat)
  genfreq <- genmat / ploidy
  sc <- ploidy / (ploidy - 1)
  
  freqsums <- tapply(alleleFreq, alleles2loc, sum)
  if(!isTRUE(all.equal(unname(as.vector(freqsums)), rep(1, nloc)))){
    stop("Allele frequencies don't add up. Was omit1allelePerLocus not set to FALSE? Is ploidy correct?")
  }
  
  He <- tapply(alleleFreq, alleles2loc,
               function(x) 1 - sum(x ^ 2))
  Ho <- (1 - t(rowsum(t(genfreq ^ 2), alleles2loc))) * sc
  out <- sweep(Ho, 2, He, "/")
  
  return(out)
}

# Number of haplotypes per genotype
# Based on Willis et al. (2017; doi: 10.1111/1755-0998.12647), count the number
# of haplotypes with more than 2 reads at each sample and locus.
# Outputs a sample x locus matrix; do colMeans(out > ploidy) as a metric for
# locus filtering.
HapPerGen <- function(countsmat, alleles2loc, minreads = 3L){
  tfmat <- t(countsmat >= minreads)
  mode(tfmat) <- "double"
  
  out <- t(rowsum(tfmat, alleles2loc))
  
  return(out)
}

# Z-score for read depth ratios
# Based on McKinney et al. (2017; doi: 10.1111/1755-0998.12613)
# Find heterozygous genotypes, sum up read counts, and see how
# much it deviates from expectations.
# To expand to multiallelic genotypes and ploidies higher than two, look just
# at genotypes where the most common allele has ploidy - 1 copies.
Zscore <- function(countsmat, genmat, alleles2loc, ploidy){
  nloc <- max(alleles2loc)
  out <- numeric(nloc)
  ht <- ploidy - 1L
  p <- ht / ploidy
  q <- 1 / ploidy
  pq <- p * q
  
  for(i in seq_len(nloc)){
    thesecol <- which(alleles2loc == i)
    thesegen <- genmat[,thesecol]
    thesecounts <- countsmat[,thesecol]
    commonAl <- which.max(colSums(thesegen == ht, na.rm = TRUE))
    hetsam <- which(thesegen[,commonAl] == ht)
    readsAl <- sum(thesecounts[hetsam,commonAl])
    readsLoc <- sum(thesecounts[hetsam,])
    sd <- sqrt(readsLoc * pq)
    z <- (p * readsLoc - readsAl) / sd
    out[i] <- z
  }
  
  return(out)
}
