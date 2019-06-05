from statistics import mean, StatisticsError
from numpy.random import choice
from scipy.stats import kendalltau
from copy import deepcopy
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

def SwapHap(NMmat, hapAssign, seqlen, corrgrps, corrP, base = 0.5):
  '''Look at the number of mutations between haplotypes and the various
  potential paralogs, and determine probabilities for moving a haplotype from
  one isolocus to another during random sampling.
  There should be a greater probability of getting moved to an isolocus where
  the reference is more similar to the haplotype.
  Additionally, if there are groups based on allele correlations, and if the
  allele to be swapped is in one of those groups, move the whole group, with
  a probability based on the p-value used to make the groups.
  Return a new version of hapAssign where one haplotype has been moved from one
  isolocus to another, at random using the above probabilities.

  NMmat is a list of lists, the first dimension corresponding to alignment
  locations and the second dimension corresponding to haplotypes.  Values
  indicate the number of mutations between that haplotype and the reference.

  hapAssign is a list of lists, with the first dimension corresponding to
  alignment.  Each list contains indices of haplotypes that have been assigned
  to that location.

  seqlen is the sequence length, and base should be about twice the maximum
  possible dissimilarity between haplotype and reference.'''

  nloc = len(NMmat) # number of isoloci
  hapAssign = deepcopy(hapAssign)
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
  # determine whether it was in a group, and whether to move the whole group
  if len(corrgrps) > 0:
    thisgrp = [grp for grp in corrgrps if thisswap[0] in grp]
    if len(thisgrp) == 1 and \
    choice([True, False], size = 1, p = [1 - corrP, corrP])[0]:
      thisgrp = thisgrp[0]
      thisgrp.difference_update({thisswap[0]}) # this allele will be moved later
      for h in thisgrp:
        [hapAssign[i].remove(h) for i in range(nloc) if h in hapAssign[i]]
        hapAssign[thisswap[2]].append(h)

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

def AlleleAssociations(countsmat):
  '''Generate a square matrix of p-values for alleles being negatively
  associated with each other.'''
  nHap = len(countsmat)
  nInd = len(countsmat[0])
  outP = [[1.0 for i in range(nHap)] for j in range(nHap)]

  for h1 in range(nHap - 1):
    # tot1 and tot2 are per individual total read depth, omitting h1 and h2, respectively
    tot1 = [sum([countsmat[h][i] for h in range(nHap) if h != h1]) for i in range(nInd)]
    for h2 in range(h1 + 1, nHap):
      tot2 = [sum([countsmat[h][i] for h in range(nHap) if h != h2]) for i in range(nInd)]
      # depth ratios that ignore the other allele being considered
      rat1 = [countsmat[h1][i]/tot2[i] if tot2[i] > 0 else 0.0 for i in range(nInd)]
      rat2 = [countsmat[h2][i]/tot1[i] if tot1[i] > 0 else 0.0 for i in range(nInd)]
      # omit ones that are both zero (missing data)
      rat1, rat2 = zip(*[(r1, r2) for r1, r2, t1, t2 in zip(rat1, rat2, tot1, tot2) if t1 > 0 or t2 > 0])
      # skip if there is too much missing data
      if len(rat1) < 10:
        continue
      # perform test for association
      kout = kendalltau(rat1, rat2, nan_policy = 'raise', method = 'asymptotic')
      # convert p-value to one-tailed
      if kout[0] <= 0:
        pout = kout[1]/2
      else:
        pout = 1 - kout[1]/2
      # add to matrix
      outP[h1][h2] = pout
      outP[h2][h1] = pout
  return(outP)

def GroupByAlAssociations(countsmat, expHindHe, startP = 0.1):
  '''Find groups of alleles that are significantly negatively associated, and
  don't exceed the expected value of Hind/He.'''
  nHap = len(countsmat)
  pvals = AlleleAssociations(countsmat)
  currP = startP

  while True:
    grps = []
    for h1 in range(nHap - 1):
      for h2 in range(h1 + 1, nHap):
        if pvals[h1][h2] > currP:
          continue # doesn't meet threshold, don't add to group
        if len(grps) == 0: # first group
          grps.append({h1, h2})
          continue
        # determine if it is already in any groups
        h1grp = [i for i in range(len(grps)) if h1 in grps[i]]
        h2grp = [i for i in range(len(grps)) if h2 in grps[i]]
        assert len(h1grp) < 2
        assert len(h2grp) < 2
        if len(h1grp) == 0 and len(h2grp) == 0:
          grps.append({h1, h2}) # new group
        elif len(h1grp) == 1 and len(h2grp) == 0:
          grps[h1grp[0]].add(h2) # add to existing group
        elif len(h1grp) == 0 and len(h2grp) == 1:
          grps[h2grp[0]].add(h1) # add to existing group
        elif len(h1grp) == 1 and len(h2grp) == 1 and h1grp[0] != h2grp[0]:
          # merge groups
          grps[h1grp[0]].update(grps[h2grp[0]])
          grps.pop(h2grp[0])
    if(len(grps) == 0):
      break # no groups could be made
    # Test Hind/He for these groups
    hindheOK = [HindHe([countsmat[h] for h in g]) <= expHindHe for g in grps]
    if all(hindheOK):
      break
    currP = currP / 10
  return([grps, currP])

def AdjustHapAssignByAlAssociations(grps, hapAssign):
  '''Adjust hapAssign if necessary so that each group that was made based on
  negative associations between alleles is in just one haplotype group.'''
  nLoc = len(hapAssign)
  for grp in grps:
    numPerHA = [sum([g in ha for g in grp]) for ha in hapAssign]
    assert sum(numPerHA) == len(grp)
    haInGrp = [i for i in range(nLoc) if numPerHA[i] > 0]
    if len(haInGrp) == 1:
      continue # no rearrangement needed
    # go to the isolocus where most of these are, or a random one.
    maxPerHA = max(numPerHA)
    matchmax = [i for i in haInGrp if numPerHA[i] == maxPerHA]
    if len(matchmax) == 1:
      targetLoc = matchmax[0]
    else:
      targetLoc = choice(matchmax, size = 1)[0]
    # rearrange haplotypes
    [hapAssign[i].remove(h) for i in haInGrp for h in grp if h in hapAssign[i]]
    hapAssign[targetLoc].extend(grp)

  return hapAssign

def HindHeExcess(hindhe, expHindHe):
  '''Return the mean amount by which per-isolocus Hind/He exceeds the expected
  value.'''
  return mean([max([0, h - expHindHe]) for h in hindhe if h != None])

def AnnealLocus(countsmat, NMmat, seqlen, expHindHe, base = 0.5, maxreps = 100,
                T0 = 0.1, rho = 0.95, corrstartP = 0.1, logcon = None):
  '''Perform simulated annealing on one group of haplotypes to split it into isoloci.
  logcon is a file connection that is already open, for logging progress.'''
  hapAssign = InitHapAssign(NMmat) # initial assignment of haplotypes to isoloci
  hindhe = HindHeByIsolocus(countsmat, hapAssign)
  NM_mean = MeanNMperLoc(NMmat, hapAssign)
  if logcon != None:
    logcon.write("Initial Hind/He: {}\n".format(" ".join([str(h) for h in hindhe])))
    logcon.write("Initial average NM: {}\n".format(NM_mean))
  # if already fixed or within Hind/He expectations at each isolocus, don't do simulated annealing
  if all([h == None or h < expHindHe for h in hindhe]):
    return hapAssign
  # get groups based on allele correlations, and adjust hapAssign if needed
  corrgrps, corrP = GroupByAlAssociations(countsmat, expHindHe, startP = corrstartP)
  hapAssign = AdjustHapAssignByAlAssociations(corrgrps, hapAssign)
  # get the mean amount by which each Hind/He is greater than expectations
  hindhe_mean = HindHeExcess(hindhe, expHindHe)

  # number of swaps to attempt per temperature
  # roughly allow each allele to get moved to each isolocus
  swapspertemp = math.factorial(len(NMmat)) * len(NMmat[0]) * 2
  Ti = T0 # current temperature

  # set aside objects to hold best solution
  hapAssign_best = deepcopy(hapAssign)
  hindhe_mean_best = deepcopy(hindhe_mean)
  NM_mean_best = deepcopy(NM_mean)

  # are we likely to explore much of the solution space?  If so, track solutions.
  nsol = len(NMmat) ** len(NMmat[0]) # number of possible solutions
  tracksol = swapspertemp * maxreps >= nsol # true if we could explore whole space
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
      hapAssign_new = SwapHap(NMmat, hapAssign, seqlen, corrgrps, corrP, base = base)
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
        hapAssign_best = hapAssign_new
        didswap = False
        if logcon != None:
          logcon.write("All isoloci fixed.\n")
        break
      if (not tracksol) or (not explored[haInd]): # Mean Hind/He minus max expected Hind/He
        hindhe_mean_new = HindHeExcess(hindhe_new, expHindHe)
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
        # update best solution if necessary
        if hindhe_mean < hindhe_mean_best or (hindhe_mean == hindhe_mean_best and NM_mean < NM_mean_best):
          hapAssign_best = deepcopy(hapAssign)
          hindhe_mean_best = deepcopy(hindhe_mean)
          NM_mean_best = deepcopy(NM_mean)
        #if logcon != None:
        #  logcon.write("Temperature: {}, Current Hind/He: {}\n".format(Ti, " ".join([str(h) for h in hindhe])))
    if not didswap:
      if hindhe_mean == hindhe_mean_best and NM_mean == NM_mean_best:
        break # no swaps happened, optimal solution found, algorithm is done
      else: # restart algorithm at best solution
        hapAssign = deepcopy(hapAssign_best)
        hindhe_mean = deepcopy(hindhe_mean_best)
        NM_mean = deepcopy(NM_mean_best)
    Ti = Ti * rho

  if logcon != None:
    logcon.write("Final Hind/He: {}\n".format(" ".join([str(h) for h in hindhe])))
    logcon.write("Final average NM: {}\n".format(MeanNMperLoc(NMmat, hapAssign)))
    logcon.write("Final temperature: {}\n".format(Ti))

  return hapAssign_best

def FindNeighbors(hapAssign, corrgrps, tabu):
  '''For the Tabu Search algorithm, find neighboring solutions to the current
  one that have not recently been examined.  corrgrps is a list of groups of
  haplotypes that can be swapped together, including all individual haplotypes
  that were not placed in groups.'''
  nLoc = len(hapAssign)
  assert sum([len(grp) for grp in corrgrps]) == sum([len(ha) for ha in hapAssign])
  haList = []
  for grp in corrgrps:
    for loc in range(nLoc):
      if not all([g in hapAssign[loc] for g in grp]):
        hapAssign_new = deepcopy(hapAssign)
        [[ha.remove(g) for g in grp if g in ha] for ha in hapAssign_new]
        hapAssign_new[loc].extend(grp)
        haInd = IndexHapAssign(hapAssign_new)
        if haInd not in tabu:
          haList.append(hapAssign_new)
  return haList

def TabuLocus(countsmat, NMmat, expHindHe, \
reps = 25, maxTabu = 5, corrstartP = 0.01, logcon = None):
  '''A Tabu Search to try to optimize first Hind/He, then NM from reference,
  while keeping together groups of alleles that are negatively associated.'''
  nHap = len(countsmat)
  nLoc = len(NMmat)
  hapAssign = InitHapAssign(NMmat) # initial assignment of haplotypes to isoloci
  hindhe = HindHeByIsolocus(countsmat, hapAssign) # doesn't get updated
  NM_mean = MeanNMperLoc(NMmat, hapAssign)        # doesn't get updated
  if logcon != None:
    logcon.write("Initial Hind/He: {}\n".format(" ".join([str(h) for h in hindhe])))
    logcon.write("Initial average NM: {}\n".format(NM_mean))
  # if already fixed or within Hind/He expectations at each isolocus, don't do search
  if all([h == None or h < expHindHe for h in hindhe]):
    return hapAssign

  # get groups based on allele correlations, and adjust hapAssign if needed
  corrgrps, corrP = GroupByAlAssociations(countsmat, expHindHe, startP = corrstartP)
  hapAssign = AdjustHapAssignByAlAssociations(corrgrps, hapAssign)
  hindhe = HindHeByIsolocus(countsmat, hapAssign)
  # if everything ok after adjusting by corrgrps, don't do search
  if all([h == None or h < expHindHe for h in hindhe]):
    return hapAssign
  # expand corrgrps to add individual alleles
  if len(corrgrps) == 0:
    corrgrps = [{h} for h in range(nHap)]
  else:
    corrgrps.extend([{h} for h in range(nHap) if all([h not in grp for grp in corrgrps])])
  nGrp = len(corrgrps)

  # set up tabu list
  tabu = [0 for i in range(maxTabu)]
  tabu[0] = IndexHapAssign(hapAssign)
  tabuIndex = 1

  # set aside objects to hold best solution
  hapAssign_best = deepcopy(hapAssign)
  hindhe_excess_best = HindHeExcess(hindhe, expHindHe)
  NM_mean_best = MeanNMperLoc(NMmat, hapAssign)
  best_rep = 0 # rep where best solution was found

  # tabu search algorithm
  for rep in range(reps):
    # get all non-tabu neighbors and choose the best
    neighbors = FindNeighbors(hapAssign, corrgrps, tabu)
    if len(neighbors) == 0:
      break
    hindhe_neighbors = [HindHeByIsolocus(countsmat, ha) for ha in neighbors]
    hindhe_excess_neighbors = [HindHeExcess(hh, expHindHe) for hh in hindhe_neighbors]
    min_excess = min(hindhe_excess_neighbors)
    min_neighbors = [neighbors[i] for i in range(len(neighbors)) \
                     if hindhe_excess_neighbors[i] == min_excess]
    if len(min_neighbors) == 1:
      hapAssign = min_neighbors[0]
    else:
      NM_neighbors = [MeanNMperLoc(NMmat, ha) for ha in min_neighbors]
      min_NM = min(NM_neighbors)
      min_neighbors = [min_neighbors[i] for i in range(len(min_neighbors)) \
                       if NM_neighbors[i] == min_NM]
      if len(min_neighbors) == 1:
        hapAssign = min_neighbors[0]
      else:
        hapAssign = min_neighbors[int(choice(range(len(min_neighbors)), size = 1))]
    # update the tabu list
    tabu[tabuIndex % maxTabu] = IndexHapAssign(hapAssign)
    tabuIndex += 1
    # update the best solution if appropriate
    min_NM = MeanNMperLoc(NMmat, hapAssign)
    if min_excess < hindhe_excess_best or (min_excess == hindhe_excess_best and min_NM < NM_mean_best):
      hapAssign_best = deepcopy(hapAssign)
      hindhe_excess_best = min_excess
      NM_mean_best = min_NM
      best_rep = rep

  if logcon != None:
    logcon.write("Final Hind/He: {}\n".format(" ".join([str(h) for h in HindHeByIsolocus(countsmat, hapAssign_best)])))
    logcon.write("Final average NM: {}\n".format(NM_mean_best))
    logcon.write("Rep where best solution found: {}\n".format(best_rep))

  return hapAssign_best
