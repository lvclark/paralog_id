#!/usr/bin/env python3
import isoloci_fun
import csv
import math

## Script to process tag depth and sort haplotypes into isoloci ##
parser = argparse.ArgumentParser(description =
'''Process the output of one chunk from process_sam_multi.py.  Estimate
Hind/He both to filter markers and to determine which groups of tags could be
adjusted in terms of assignment of tags to isoloci.  For groups of tags needing
adjustment, analyze allele correlations to identify groups that putatively
belong to the same isolocus, then perform tabu search to find groups of tags
within Hind/He expections that minimize number of mutations from the reference
sequence.  Output read depth and assignment of alleles to loci, for import by
polyRAD.''',
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("alignfile", nargs = '?',
                    help = "Path to alignment file output by process_sam_multi.py.")
parser.add_argument("depthfile", nargs = '?',
                    help = "Path to read depth file output by process_sam_multi.py.")
parser.add_argument("out", nargs = '?',
                    help = "Base file name for output.") # tack on number from input files if present
parser.add_argument("--ploidy", "-p", nargs = '?', type = int, default = 2,
                    help = "Expected ploidy after splitting isoloci.")
parser.add_argument("--inbreeding", "-f", nargs = '?', type = float, default = 0.0,
                    help = "Inbreeding coefficient, ranging from 0 to 1.")
parser.add_argument("--logfile", "-l", nargs = '?', defaults = "",
                    help = "Optional path to file where log should be written.")

args = parser.args()
alignfile = args.alignfile
depthfile = args.depthfile
outbase = args.out
ploidy = args.ploidy
inbreeding = args.inbreeding
logfile = args.logfile

# maximum tolerable Hind/He: halfway between this and the next ploidy, on a log scale
p2 = ploidy * 2
maxHindHe = math.exp((math.log((ploidy - 1)/ploidy) + math.log((p2 - 1)/p2) + 2 * math.log(1 - inbreeding))/2)
expHindHe = (ploidy - 1)/ploidy * (1 - inbreeding)

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
  #seqlen = max([len(row[0]) for row in depthrows]) # maximum tag length

  #hapAssign = isoloci_fun.AnnealLocus(depths, NM, seqlen, expHindHe, T0 = 0.10, rho = 0.95, logcon = logcon)
  hapAssign = isoloci_fun.TabuLocus(depths, NM, expHindHe, logcon = logcon)
  return None

# files need to already have tags in the same order, grouped by alignment location
# loop through markers
try:
  depthcon = open(depthfile, newline = '', mode = 'r')
  aligncon = open(alignfile, newline = '', mode = 'r')
  if logfile = "":
    logcon = None
  else:
    logcon = open(logfile, mode = 'w')
  depthreader = csv.reader(depthcon)
  alignreader = csv.reader(aligncon)

  # header info from files
  depthheader = next(depthreader)
  samples = depthheader[1:]
  alignheader = next(alignreader)
  maxisoloci = sum([h.startswith("Alignment") for h in alignheader])
  if maxisoloci == 0:
    raise Exception("Columns starting with 'Alignment' not found in " + alignfile)

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
    #    break
    if rowcount % 1000 == 0:
      print(rowcount)
  ProcessRowGroup(curralignrows, currdepthrows, maxisoloci, maxHindHe, expHindHe, logcon)
finally:
  depthcon.close()
  aligncon.close()
  logcon.close()
