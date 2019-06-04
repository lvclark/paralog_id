#!/usr/bin/env python3
import isoloci_fun
import csv
import math

## Script to process tag depth and sort haplotypes into isoloci ##

# variables to set up as command line args later
maxisoloci = 2  # how many subgenomes are there
ploidy = 2      # expected ploidy after sorting
alignfile = "../marker_CSV/190525twoalign_Chr1Chr2.csv" # alignment locations
depthfile = "../marker_CSV/190523diploid_Chr1Chr2.csv"  # read depth
logfile = "../log/190531hindhelogT10_5xSwaps.txt"

# maximum tolerable Hind/He: halfway between this and the next ploidy, on a log scale
p2 = ploidy * 2
maxHindHe = math.exp((math.log((ploidy - 1)/ploidy) + math.log((p2 - 1)/p2))/2)
expHindHe = (ploidy - 1)/ploidy

def ProcessRowGroup(alignrows, depthrows, nisoloci, thresh, expHindHe, logcon):
  '''Process two matching groups of rows showing alignment and depth for a
  group of tags corresponding to one set of alignment locations.'''
  # write marker being analyzed to log
  logcon.write(" ".join(alignrows[0][:nisoloci]) + "\n")
  assert len(depthrows) > 1
  
  depths = [[int(d) for d in row[1:]] for row in depthrows] # integer depths
  # filter out any alleles without depth
  packed = [(dep, ar) for dep, ar in zip(depths, alignrows) if not all([d < 2 for d in dep])]
  if len(packed) < 2:
    logcon.write("Insufficient read depth.\n")
    return None
  depths, alignrows = zip(*packed)
  # check if the marker should be split
  hindheUnsplit = isoloci_fun.HindHe(depths)
  logcon.write("Hind/He without splitting: {}\n".format(hindheUnsplit))
  if hindheUnsplit <= thresh:
    return None

  # number of mutations from reference locations
  NM = [[int(row[nisoloci + 1 + i]) for row in alignrows] for i in range(nisoloci)]
  seqlen = max([len(row[0]) for row in depthrows]) # maximum tag length
  
  hapAssign = isoloci_fun.AnnealLocus(depths, NM, seqlen, expHindHe, T0 = 0.10, logcon = logcon)
  return None

# files need to already have tags in the same order, grouped by alignment location
# loop through markers
try:
  depthcon = open(depthfile, newline = '', mode = 'r')
  aligncon = open(alignfile, newline = '', mode = 'r')
  logcon = open(logfile, mode = 'w')
  depthreader = csv.reader(depthcon)
  alignreader = csv.reader(aligncon)
  
  # header info from tag file
  header = next(depthreader)
  samples = header[1:]
  
  currdepthrows = [next(depthreader)]
  curralignrows = [next(alignreader)]
  
  rowcount = 0 # for testing, to prevent going thru whole file
  
  for row in alignreader:
    newalignrow = row
    newdepthrow = next(depthreader)
    tag = newdepthrow[0] # tag sequence
    assert newalignrow[maxisoloci] == tag
    if newalignrow[:maxisoloci] == curralignrows[0][:maxisoloci]: # same alignment
      currdepthrows.append(newdepthrow)
      curralignrows.append(newalignrow)
    else: # new alignment; process last one and start new
      ProcessRowGroup(curralignrows, currdepthrows, maxisoloci, maxHindHe, expHindHe, logcon)
      currdepthrows = [newdepthrow]
      curralignrows = [newalignrow]
    rowcount += 1
    # if rowcount > 5000: # for testing only
    #   break
    if rowcount % 1000 == 0:
      print(rowcount)
  ProcessRowGroup(curralignrows, currdepthrows, maxisoloci, maxHindHe, expHindHe, logcon)
finally:
  depthcon.close()
  aligncon.close()
  logcon.close()
