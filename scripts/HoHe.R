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
