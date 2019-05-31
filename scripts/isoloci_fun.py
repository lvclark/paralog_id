from statistics import mean, StatisticsError
from numpy.random import choice
import math

# functions for sorting out isoloci

def GiniSimpson(counts, N = None):
  "Estimate the Gini-Simpson index, not corrected for size, on a list of integers."
  if N == None:
    N = sum(counts)
  if N == 0:
    return None
  for c in counts: # check to see if it should be 0, before doing calc.
    if c == N:
      return 0.0
    elif c != 0:
      break
  freq = [c/N for c in counts]
  return 1.0 - sum([f**2 for f in freq])

def HindHe(countsmat):
  '''For one marker, estimate Hind/He.
  countsmat is a list of lists, with the first dimension representing alleles
  and the second dimension representing taxa. Each cell is read depth.
  Hind/He will be corrected for size on a per-taxon basis by multiplying by
  N/(N-1), then averaged across taxa.'''
  if len(countsmat) < 2:
    return None

  countsmatT = list(zip(*countsmat)) # transpose the matrix
  depthByInd = [sum(x) for x in countsmatT]
  nind = len(depthByInd)
  if all([d == 0 for d in depthByInd]):
    return None
  
  depthRatios = [[c / depthByInd[i] for c in countsmatT[i]] for i in range(nind) if depthByInd[i] > 0]
  meanDepthRatios = [sum(x)/nind for x in zip(*depthRatios)] # mean by allele
  He = GiniSimpson(meanDepthRatios, N = 1)
  assert He > 0
  
  HindHeByInd = [GiniSimpson(countsmatT[i], N = depthByInd[i]) * \
                 depthByInd[i] / (depthByInd[i] - 1)/ He for i in range(nind) if depthByInd[i] > 1]

  m = sum(HindHeByInd)/nind
  return m
  
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
  
def MeanNMperLoc(NMmat, hapAssign):
  '''Get mean number of mutations per locus across all haplotypes versus their
  assigned locus.'''
  nLoc = len(NMmat)
  nHap = len(NMmat[0])
  out = sum([NMmat[L][h] for L in range(nLoc) for h in hapAssign[L]])/nHap
  return out
  
def IndexHapAssign(hapAssign):
  '''Generate an index for a given set of assignments of haplotypes to isoloci,
  so that we can track solutions that have already been examined.'''
  nLoc = len(hapAssign)
  return sum([i * nLoc ** h for i in range(nLoc) for h in hapAssign[i]])

def AnnealLocus(countsmat, NMmat, seqlen, expHindHe, base = 0.5, maxreps = 100,
                T0 = 0.1, rho = 0.95, logcon = None):
  '''Perform simulated annealing on one group of haplotypes to split it into isoloci.
  logcon is a file connection that is already open, for logging progress.'''
  hapAssign = InitHapAssign(NMmat) # initial assignment of haplotypes to isoloci
  hindhe = HindHeByIsolocus(countsmat, hapAssign)
  NM_mean = MeanNMperLoc(NMmat, hapAssign)
  if logcon != None:
    logcon.write("Initial Hind/He: {}\n".format(" ".join([str(h) for h in hindhe])))
    logcon.write("Initial average NM: {}\n".format(NM_mean))
  # if already fixed at each isolocus, don't do simulated annealing
  if all([h == None or h < expHindHe for h in hindhe]):
    return hapAssign
  # get the mean amount by which each Hind/He is greater than expectations
  hindhe_mean = mean([max([0, h - expHindHe]) for h in hindhe if h != None])
  
  # number of swaps to attempt per temperature
  # roughly allow each allele to get moved to each isolocus
  swapspertemp = math.factorial(len(NMmat)) * len(NMmat[0])
  Ti = T0 # current temperature
  
  # set aside objects to hold best solution
  hapAssign_best = hapAssign
  hindhe_mean_best = hindhe_mean
  NM_mean_best = NM_mean
  
  # are we likely to explore much of the solution space?  If so, track solutions.
  nsol = len(NMmat) ** len(NMmat[0]) # number of possible solutions
  tracksol = swapspertemp * maxreps <= nsol # true if we could explore whole space
  if tracksol:
    explored = [False for i in range(nsol)]
    all_hindhe = [[] for i in range(nsol)]
    all_hindhe_mean = [0 for i in range(nsol)]
    haInd = IndexHapAssign(hapAssign)
    explored[haInd] = True
    all_hindhe[haInd] = hindhe
    all_hindhe_mean[haInd] = hindhe_mean
  
  # start algorithm
  for k in range(maxreps):
    didswap = False
    for r in range(swapspertemp):
      hapAssign_new = SwapHap(NMmat, hapAssign, seqlen, base = base)
      if tracksol: # solution tracking
        haInd = IndexHapAssign(hapAssign)
        if explored[haInd]: # solution already explored
          hindhe_new = all_hindhe[haInd]
          hindhe_mean_new = all_hindhe_mean[haInd]
      if (not tracksol) or (not explored[haInd]): # Hind/He per isolocus
        hindhe_new = HindHeByIsolocus(countsmat, hapAssign_new)
        if tracksol:
          all_hindhe[haInd] = hindhe_new
      if all([h == None for h in hindhe_new]): # if we get all isoloci fixed, stop algorithm
        hapAssign = hapAssign_new
        didswap = False
        if logcon != None:
          logcon.write("All isoloci fixed.\n")
        break
      if (not tracksol) or (not explored[haInd]): # Mean Hind/He minus max expected Hind/He
        hindhe_mean_new = mean([max([0, h - expHindHe]) for h in hindhe_new if h != None])
        if tracksol:
          all_hindhe_mean[haInd] = hindhe_mean_new
          explored[haInd] = True
      doswap = hindhe_mean_new <= hindhe_mean # the new set is better
      if not doswap:
        # new set is worse, decide whether to keep
        swapprob = math.exp((hindhe_mean - hindhe_mean_new)/Ti)
        doswap = choice(range(2), size = 1, p = [swapprob, 1 - swapprob])[0] == 0
      if hindhe_mean_new == 0 and hindhe_mean == 0:
        # both the old and new have gotten below the expected value.
        # decide whether to swap based on whether the number of mutations is lower.
        doswap = MeanNMperLoc(NMmat, hapAssign_new) < MeanNMperLoc(NMmat, hapAssign)
      if doswap: # everything that happens if we do a swap
        didswap = True
        hapAssign = hapAssign_new
        hindhe = hindhe_new
        hindhe_mean = hindhe_mean_new
        NM_mean = MeanNMperLoc(NMmat, hapAssign)
        # update best solution of necessary
        if hindhe_mean < hindhe_mean_best or (hindhemean == hindhe_mean_best and NM_mean < NM_mean_best):
          hapAssign_best = hapAssign
          hindhe_mean_best = hindhe_mean
          NM_mean_best = NM_mean
        #if logcon != None:
        #  logcon.write("Temperature: {}, Current Hind/He: {}\n".format(Ti, " ".join([str(h) for h in hindhe])))
    if not didswap:
      if hindhe_mean == hindhe_mean_best and NM_mean == NM_mean_best:
        break # no swaps happened, optimal solution found, algorithm is done
      else: # restart algorithm at best solution
        hapAssign = hapAssign_best
        hindhe_mean = hindhe_mean_best
        NM_mean = NM_mean_best
    Ti = Ti * rho
  
  if logcon != None:
    logcon.write("Final Hind/He: {}\n".format(" ".join([str(h) for h in hindhe])))
    logcon.write("Final average NM: {}\n".format(MeanNMperLoc(NMmat, hapAssign)))
    logcon.write("Final temperature: {}\n".format(Ti))
  
  return hapAssign
