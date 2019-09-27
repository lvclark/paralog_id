# functions relating to getting Hind/He in mapping populations

# grab some internal functions, to make it easy to incorporate into polyRAD later
.makeGametes <- polyRAD:::.makeGametes
.gameteProb <- polyRAD:::.gameteProb
.progenyProb <- polyRAD:::.progenyProb

# Function to get probabilities that, if two alleles are sampled from a progeny,
# they will both be from parent 1, both from parent 2, or from different parents.
# If there is backcrossing, parent 1 is the recurrent parent.
.progAlProbs <- function(ploidy, gen_backcrossing, gen_selfing){
  if(ploidy %% 2 != 0){
    stop("Even ploidy needed.")
  }
  p2 <- ploidy/2
  # Frequencies of all possible gametes from all possible genotypes
  allgamprobs <- .gameteProb(.makeGametes(0:ploidy, ploidy), ploidy)
  
  # first get frequencies of progeny with each possible allele copy number from parent 2.
  # F1: all individuals are 1:1 parent 1 and parent 2.
  progfreqs <- matrix(c(rep(0, p2), 1, rep(0, p2)), nrow = ploidy + 1, ncol = 1)
  # Backcrossing
  if(gen_backcrossing > 0){
    # parent 1 only produces gametes with 0 copies of parent 2 alleles
    p1gamprobs <- matrix(c(1, rep(0, p2)), nrow = p2 + 1, ncol = ploidy + 1)
    # progeny frequency after crossing any genotype with parent 1
    allBCprogprobs <- .progenyProb(p1gamprobs, allgamprobs)
    for(g in 1:gen_backcrossing){
      progfreqs <- allBCprogprobs %*% progfreqs
    }
  }
  # Selfing
  if(gen_selfing > 0){
    # progeny frequency after self-fertilizing any genotype
    allSelfprogprobs <- .progenyProb(allgamprobs, allgamprobs)
    for(g in 1:gen_selfing){
      progfreqs <- allSelfprogprobs %*% progfreqs
    }
  }
  
  # For each possible genotype, probability of sampling, with replacement,
  # two alleles from p1, two alleles from p2, one allele from each (in that order).
  allprobmat <- matrix(c((ploidy:0) / ploidy * c((ploidy - 1):0, 0) / (ploidy - 1),
                         (0:ploidy) / ploidy * c(0, 0:(ploidy - 1)) / (ploidy - 1),
                         (ploidy:0) / ploidy * c(0:(ploidy - 1), ploidy - 1) / (ploidy - 1) +
                           (0:ploidy)/ ploidy * c(ploidy - 1, (ploidy - 1):0) / (ploidy - 1)),
                       byrow = TRUE, nrow = 3, ncol = ploidy + 1)
  return(allprobmat %*% progfreqs)
}
