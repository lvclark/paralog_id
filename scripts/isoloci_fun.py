from statistics import mean
from numpy.random import choice

# functions for sorting out isoloci

def GiniSimpson(counts, N = None):
  "Estimate the Gini-Simpson index, not corrected for size, on a list of integers."
  if N == None:
    N = sum(counts)
  if N == 0:
    return None
  freq = [c/N for c in counts]
  return 1.0 - sum([f**2 for f in freq])

def HindHe(countsmat):
  '''For one marker, estimate Hind/He.
  countsmat is a list of lists, with the first dimension representing alleles
  and the second dimension representing taxa. Each cell is read depth.
  Hind/He will be corrected for size on a per-taxon basis by multiplying by
  N/(N-1), then averaged across taxa.'''
  countsmatT = list(zip(*countsmat)) # transpose the matrix
  depthByInd = [sum(x) for x in countsmatT]
  nind = len(depthByInd)
  if all([d == 0 for d in depthByInd]):
    return None
  
  depthRatios = [[c / depthByInd[i] for c in countsmatT[i]] for i in range(nind) if depthByInd[i] > 0]
  meanDepthRatios = [mean(x) for x in zip(*depthRatios)] # mean by allele
  He = GiniSimpson(meanDepthRatios, N = 1)
  
  HindHeByInd = [GiniSimpson(countsmatT[i], N = depthByInd[i]) * \
                 depthByInd[i] / (depthByInd[i] - 1)/ He for i in range(nind) if depthByInd[i] > 0]
  return mean(HindHeByInd)
  
def SwapHap(NMmat, hapAssign, seqlen, base = 0.5):
  '''Look at the number of mutations between haplotypes and the various
  potential paralogs, and determine probabilities for moving a haplotype from
  one isolocus to another during random sampling.
  There should be a greater probability of getting moved to an isolocus where
  the reference is more similar to the haplotype.
  Return a new version of hapAssign where one haplotype has been moved from one
  isolocus to another, at random using the above probabilities.
  
  NMmat is a list of lists, the first dimension corresponding to alignment
  locations and the second dimension corresponding to haplotypes.  Values
  indicate the number of mutations between that haplotype and the reference.
  
  hapAssign is a list of lists, with the first dimension corresponding to
  alignment.  Each list contains indices of haplotypes that have been assigned
  to that location.
  
  seqlen is the sequence length, and base should be about twice the maximum
  possibley dissimilarity between haplotype and reference.'''
  
  nloc = len(NMmat) # number of isoloci
  # Build a list of potential movements of haplotypes from one isolocus to
  # another, and their weights.
  swaplist = []
  
  for loc1 in range(nloc): # current locus
    otherloc = [l for l in range(nloc) if l != loc1] # loci other than this one
    for hap in hapAssign[loc1]: # haplotype index
      for loc2 in otherloc: # destination locus
        swapval = base - (NMmat[loc2][hap] - NMmat[loc1][hap]) / seqlen
        swaplist.append([hap, loc1, loc2, swapval])

  # normalize weights to probabilities
  weights = [s[3] for s in swaplist]
  tot = sum(weights)
  probs = [w/tot for w in weights]
  
  # pick a swap to make
  thisswap = swaplist[choice(len(swaplist), size = 1, p = probs)[0]]
  # move the alleles
  hapAssign[thisswap[1]].remove(thisswap[0])
  hapAssign[thisswap[2]].append(thisswap[0])
  
  return hapAssign
