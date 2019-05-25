#!/usr/bin/env python3
import isoloci_fun
import csv

## Script to process tag depth and sort haplotypes into isoloci ##

# variables to set up as command line args later
maxisoloci = 2  # how many subgenomes are there
ploidy = 2      # expected ploidy after sorting
alignfile = "../marker_CSV/190517twoalign.csv"         # alignment locations
depthfile = "../marker_CSV/190523diploid_Chr1Chr2.csv" # read depth

# import depth
with open(depthfile, newline = '') as depthcon:
  depthreader = csv.reader(depthcon)
  header = next(depthreader)
  samples = header[1:]
  tags = []
  depths = []
  for row in depthreader:
    tags.append(row[0])
    depths.append([int(d) for d in row[1:]])
    
# import alignement data
with open(alignfile, newline = '') as aligncon:
  alignreader = csv.reader(aligncon)
  markernames = []
  markerNM = []
  markertagindex = []
  for row in alignreader:
    tag = row[maxisoloci]
    if tag not in tags:
      continue
    thisalign = tuple(row[:maxisoloci])
    thisNM = [int(nm) for nm in row[maxisoloci + 1:]]
    if len(markernames) > 0 and markernames[-1] == thisalign:
      # same marker as last tag read
      markertagindex[-1].append(tags.index(tag)) # index of tag in depth matrix
      [markerNM[-1][i].append(thisNM[i]) for i in range(maxisoloci)]
    else: # new marker
      markernames.append(thisalign)
      markertagindex.append([tags.index(tag)])
      markerNM.append([[nm] for nm in thisNM])

for i in range(5):
  print(markernames[i])
  print([tags[j] for j in markertagindex[i]])
  print(markerNM[i])
