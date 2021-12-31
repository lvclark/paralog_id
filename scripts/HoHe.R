# Functions for other paralog detection approaches

# Observed over Expected heterozygosity
HoHe <- function(genmat, alleles2loc, ploidy){
  alleleFreq <- colMeans(genmat) / ploidy
  nloc <- max(alleles2loc)
  nsam <- nrow(genmat)
  genfreq <- genmat / ploidy
  sc <- ploidy / (ploidy - 1)
  
  freqsums <- tapply(alleleFreq, alleles2loc, sum)
  if(!isTRUE(all.equal(unname(as.vector(freqsums)), rep(1, nloc)))){
    stop("Allele frequencies don't add up. Was omit1allelePerLocus not set to FALSE? Is ploidy correct?")
  }
  
  out <- matrix(NA_real_, nrow = nsam, ncol = nloc)
  
  # Consider implementation in Rcpp if using for more than just this test
  for(L in seq_len(nloc)){
    thesecol <- which(alleles2loc == L)
    He <- 1 - sum(alleleFreq[thesecol] ^ 2)
    for(s in seq_len(nsam)){
      # probability of drawing two different alleles, without replacement
      Ho <- (1 - sum(genfreq[s,thesecol] ^ 2)) * sc
      out[s,L] <- Ho / He
    }
  }
  
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
