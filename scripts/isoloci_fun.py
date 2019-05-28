from statistics import mean
from numpy.random import choice
import math

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
                 depthByInd[i] / (depthByInd[i] - 1)/ He for i in range(nind) if depthByInd[i] > 1]
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

def InitHapAssign(NMmat):
  '''Generate initial assignments of haplotypes to loci based on number of
  mutations between haplotype and the locus.'''
  nloc = len(NMmat)
  nhap = len(NMmat[0])
  hapAssign = [[] for i in range(nloc)] # to store output
  
  for h in range(nhap):
    theseNM = [nm[h] for nm in NMmat]
    minNM = min(theseNM)
    bestLoc = [l for l in range(nloc) if theseNM[l] == minNM]
    if len(bestLoc) > 1:
      bestLoc = choice(bestLoc, 1)
    bestLoc = bestLoc[0] # convert list to number
    hapAssign[bestLoc].append(h) # add haplotype to locus
  
  return hapAssign
  
def HindHeByIsolocus(countsmat, hapAssign):
  '''For a given set of assignments of haplotypes to isoloci, estimate Hind/He
  for each isolocus.  countsmat and hapAssign are as defined above.'''
  splitcounts = [[countsmat[h] for h in isolocus] for isolocus in hapAssign]
  return [HindHe(c) for c in splitcounts]

def AnnealLocus(countsmat, NMmat, seqlen, base = 0.5, maxreps = 100, T0 = 0.5,
                rho = 0.95, logcon = None):
  '''Perform simulated annealing on one group of haplotypes to split it into isoloci.
  logcon is a file connection that is already open, for logging progress.'''
  hapAssign = InitHapAssign(NMmat) # initial assignment of haplotypes to isoloci
  hindhe = HindHeByIsolocus(countsmat, hapAssign)
  if logcon != None:
    logcon.write("Initial Hind/He: {}\n".format(" ".join([str(h) for h in hindhe])))
  hindhe_mean = mean([h for h in hindhe if h != None])
  
  # number of swaps to attempt per temperature
  # roughly allow each allele to get moved to each isolocus
  swapspertemp = math.factorial(len(NMmat)) * len(NMmat[0])
  Ti = T0
  
  for k in range(maxreps):
    didswap = False
    for r in range(swapspertemp):
      hapAssign_new = SwapHap(NMmat, hapAssign, seqlen, base = base)
      hindhe_new = HindHeByIsolocus(countsmat, hapAssign_new)
      hindhe_mean_new = mean([h for h in hindhe_new if h != None])
      doswap = hindhe_mean_new <= hindhe_mean # the new set is better
      if not doswap:
        # new set is worse, decide whether to keep
        swapprob = math.exp((hindhe_mean - hindhe_mean_new)/Ti)
        doswap = choice(range(2), size = 1, p = [swapprob, 1 - swapprob])[0] == 0
      if doswap: # everything that happens if we do a swap
        didswap = True
        hapAssign = hapAssign_new
        hindhe = hindhe_new
        hindhe_mean = hindhe_mean_new
        if logcon != None:
          logcon.write("Temperature: {}, Current Hind/He: {}\n".format(Ti, " ".join(hindhe)))
    if not didswap:
      break # no swaps happened, algorithm is done
    Ti = Ti * rho
  
  return hapAssign
