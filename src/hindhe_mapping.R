# functions relating to getting Hind/He in mapping populations
library(polyRAD) # 1.2 with some functions updated and made separate

# Function for making multiallelic gametes.  Autopolyploid segretation only.
# geno is a vector of length ploidy with values indicating allele identity.
# rnd is how many alleles will be in the gamete.
.makeGametes2 <- function(geno, rnd = length(geno)/2){
  nal <- length(geno) # number of remaining alleles in genotype
  if(rnd == 1){
    return(matrix(geno, nrow = nal, ncol = 1))
  } else {
    outmat <- matrix(0L, nrow = choose(nal, rnd), ncol = rnd)
    currrow <- 1
    for(i in 1:(nal - rnd + 1)){
      submat <- .makeGametes2(geno[-(1:i)], rnd - 1)
      theserows <- currrow + 1:nrow(submat) - 1
      currrow <- currrow + nrow(submat)
      outmat[theserows,1] <- geno[i]
      outmat[theserows,2:ncol(outmat)] <- submat
    }
    return(outmat)
  }
}
# take output of .makeGametes2 and get progeny probabilities.
# Return a list, where the first item is a matrix with progeny genotypes in
# rows, and the second item is a vector of corresponding probabilities.
.progenyProb2 <- function(gam1, gam2){
  allprog <- cbind(gam1[rep(1:nrow(gam1), each = nrow(gam2)),],
                   gam2[rep(1:nrow(gam2), times = nrow(gam1)),])
  allprog <- t(apply(allprog, 1, sort))
  # set up output
  outprog <- matrix(allprog[1,], nrow = 1, ncol = ncol(allprog))
  outprob <- 1/nrow(allprog)
  # tally up unique genotypes
  for(i in 2:nrow(allprog)){
    thisgen <- allprog[i,]
    existsYet <- FALSE
    for(j in 1:nrow(outprog)){
      if(identical(thisgen, outprog[j,])){
        outprob[j] <- outprob[j] + 1/nrow(allprog)
        existsYet <- TRUE
        break
      }
    }
    if(!existsYet){
      outprog <- rbind(outprog, thisgen, deparse.level = 0)
      outprob <- c(outprob, 1/nrow(allprog))
    }
  }
  return(list(outprog, outprob))
}
# Consolidate a list of outputs from .progenyProb2 into a single output.
# It is assumed the all probabilities in pplist have been multiplied by
# some factor so they sum to one.
.consolProgProb <- function(pplist){
  outprog <- pplist[[1]][[1]]
  outprob <- pplist[[1]][[2]]
  if(length(pplist) > 1){
    for(i in 2:length(pplist)){
      for(pr1 in 1:nrow(pplist[[i]][[1]])){
        thisgen <- pplist[[i]][[1]][pr1,]
        existsYet <- FALSE
        for(pr2 in 1:nrow(outprog)){
          if(identical(thisgen, outprog[pr2,])){
            outprob[pr2] <- outprob[pr2] + pplist[[i]][[2]][pr1]
            existsYet <- TRUE
            break
          }
        }
        if(!existsYet){
          outprog <- rbind(outprog, thisgen, deparse.level = 0)
          outprob <- c(outprob, pplist[[i]][[2]][pr1])
        }
      }
    }
  }
  
  return(list(outprog, outprob))
}
# Probability that, depending on the generation, two alleles sampled in a progeny
# are different locus copies from the same parent, or from different parents.
# (No double reduction.)
# Can take a few seconds to run if there are many generations, but it is only
# intended to be run once for the whole dataset.
.progAlProbs <- function(ploidy, gen_backcrossing, gen_selfing){
  # set up parents; number indexes locus copy
  p1 <- 1:ploidy
  p2 <- p1 + ploidy
  # create F1 progeny probabilities
  progprob <- .progenyProb2(.makeGametes2(p1), .makeGametes2(p2))
  # backcross
  if(gen_backcrossing > 0){
    gam1 <- .makeGametes2(p1)
    for(g in 1:gen_backcrossing){
      allprogprobs <- lapply(1:nrow(progprob[[1]]),
                             function(x){
                               pp <- .progenyProb2(gam1, .makeGametes2(progprob[[1]][x,]))
                               pp[[2]] <- pp[[2]] * progprob[[2]][x]
                               return(pp)
                             })
      progprob <- .consolProgProb(allprogprobs)
    }
  }
  # self-fertilize
  if(gen_selfing > 0){
    for(g in 1:gen_selfing){
      allprogprobs <- lapply(1:nrow(progprob[[1]]),
                             function(x){
                               gam <- .makeGametes2(progprob[[1]][x,])
                               pp <- .progenyProb2(gam, gam)
                               pp[[2]] <- pp[[2]] * progprob[[2]][x]
                               return(pp)
                             })
      progprob <- .consolProgProb(allprogprobs)
    }
  }
  
  # total probability that (without replacement, from individual progeny):
  diffp1 <- 0 # two different locus copies, both from parent 1, are sampled
  diffp2 <- 0 # two different locus copies, both from parent 2, are sampled
  diff12 <- 0 # locus copies from two different parents are sampled
  
  ncombo <- choose(ploidy, 2) # number of ways to choose 2 alleles from a genotype
  # examine each progeny genotype
  for(p in 1:nrow(progprob[[1]])){
    thisgen <- progprob[[1]][p,]
    thisprob <- progprob[[2]][p]
    for(m in 1:(ploidy - 1)){
      al1 <- thisgen[m]
      for(n in (m + 1):ploidy){
        al2 <- thisgen[n]
        if(al1 == al2) next
        if((al1 %in% p1) && (al2 %in% p1)){
          diffp1 <- diffp1 + (thisprob / ncombo)
          next
        }
        if((al1 %in% p2) && (al2 %in% p2)){
          diffp2 <- diffp2 + (thisprob / ncombo)
          next
        }
        if(((al1 %in% p1) && (al2 %in% p2)) ||
           ((al2 %in% p1) && (al1 %in% p2))){
          diff12 <- diff12 + (thisprob / ncombo)
          next
        }
        stop("Allele indexing messed up.")
      }
    }
  }
  
  return(c(diffp1, diffp2, diff12))
}

# Function to get Hind/He matrix in a mapping population
HindHeMapping <- function(object, ...){
  UseMethod("HindHeMapping", object)
}
HindHeMapping.RADdata <- function(object, n.gen.backcrossing = 0,
                                  n.gen.intermating = 0, n.gen.selfing = 0,
                                  ploidy = object$possiblePloidies[[1]],
                                  minLikelihoodRatio = 10){
  if(length(ploidy) != 1){
    stop("Current implementation for autopolyploids only")
  }
  if(n.gen.intermating > 0){
    stop("If the most recent generation was produced by random mating among progeny, use HindHe instead.")
  }
  donorParent <- GetDonorParent(object)
  recurrentParent <- GetRecurrentParent(object)
  object <- EstimateParentalGenotypes(object, donorParent = donorParent,
                                      recurrentParent = recurrentParent,
                                      n.gen.backcrossing = n.gen.backcrossing,
                                      n.gen.selfing = n.gen.selfing,
                                      n.gen.intermating = n.gen.intermating,
                                      minLikelihoodRatio = minLikelihoodRatio,
                                      donorParentPloidies = ploidy,
                                      recurrentParentPloidies = ploidy)
  # Get probabilities of pairs of alleles from a random progeny coming from
  # different locus copies one parent or the other, or from different parents.
  progAlProbs <- .progAlProbs(ploidy, n.gen.backcrossing, n.gen.selfing)
  
  # Identify loci where multiallelic genotypes can be determined
  goodLocDon <- tapply(object$likelyGeno_donor[as.character(ploidy),],
                       object$alleles2loc,
                       function(x) !any(is.na(x)) && sum(x) == ploidy)
  goodLocRec <- tapply(object$likelyGeno_recurrent[as.character(ploidy),],
                       object$alleles2loc,
                       function(x) !any(is.na(x)) && sum(x) == ploidy)
  keeploc <- which(goodLocDon & goodLocRec)
  
  # Get within- and across- parent probabilties of sampling two different alleles.
  
  
}
