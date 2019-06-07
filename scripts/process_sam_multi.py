#!/usr/bin/env python3

# Process a SAM file of alignment of M. sacchariflorus tags to the M. sinensis
# reference.  Output tags that align once, and tags that align twice.  Include
# NM, the number of mutational steps from the reference for each alignment
# location.
# Bowtie2 was run with -k 8 to generate the SAM file.

import csv
import re
import argparse

parser = argparse.ArgumentParser(description =
'''Process a SAM file reporting multiple alignments, and output one or more CSV
files that list tag sequences organized by sets of aligment locations.  A
TagTaxaDist file output by TASSEL is also processed to extract read depth for
the same tags.''')
parser.add_argument("sam", nargs = 1, help = "Path to SAM file.")
parser.add_argument("ttd", nargs = 1, help = "Path to TagTaxaDist file.")
parser.add_argument("out", nargs = 1, help = "Base file name for output.")
parser.add_argument("--subgenomes", "-g", nargs = 1, type = int, default = 2,
                    help = "Number of subgenomes, i.e. maximum number of alignments expected per tag.")
parser.add_argument("--chunks", "-c", nargs = 1, type = int, default = 1,
                    help = "Number of files to split the output into.")
args = parser.parse_args()
mysam = args.sam
myttd = args.ttd
outbase = args.out
maxalign = args.subgenomes
nchunks = args.chunks

#mysam = "D:/TASSELGBS_Msa/190517aligned_tags_multi.sam"
#maxalign = 2 # maximum number of alignments allowed per tag
onealign_file = "marker_CSV/190517onealign.csv"
twoalign_file = "marker_CSV/190517twoalign.csv"

# dictionaries to store tags aligning once and twice
onealign = dict()
twoalign = dict()
count_one = 0
count_two = 0

# lists to store marker data for the current tag
these_mnames = [] # marker names
these_NM = [] # number of mutational steps between tag and reference

with open(mysam, mode = "r") as samcon:
  for line in samcon:
    if line[0] == "@": # header lines
      continue
    row = line.split()

    flag = int(row[1])
    if flag & 4 == 4: # no alignment
        continue

    # extract various information about the alignment
    tagseq = row[0][7:]
    chrom = row[2]
    pos = int(row[3])
    cigar = row[5]
    if flag & 16 == 16:
      strand = "bot"
      deletions = sum([int(x[:-1]) for x in re.findall('\d+D', cigar)])
      insertions = sum([int(x[:-1]) for x in re.findall('\d+I', cigar)])
      pos = pos + len(tagseq) - insertions + deletions - 1
    else:
      strand = "top"
    mname = "{}-{:0>{width}}-{}".format(chrom, pos, strand, width=9)
    if row[17].startswith("NM:i:"):
      NM = int(row[17][5:])
    else:
      assert row[16].startswith("NM:i:")
      NM = int(row[16][5:])

    secondary = flag & 256 == 256
    if secondary: # add alignments to the list
      these_mnames.append(mname)
      these_NM.append(NM)
    else: # start a new tag
      # put the last tag into single alignment dict if applicable
      if(len(these_mnames) == 1):
        mname_out = these_mnames[0]
        if mname_out not in onealign.keys():
          onealign[mname_out] = []
        onealign[mname_out].append((lasttagseq, these_NM[0]))
        count_one += 1
        if count_one % 5000 == 0:
          print("{} tags aligning once".format(count_one))
        # if count_one == 5000:
        #   break ### temporary for testing
      # put the last tag into double alignment dict if applicable
      if(len(these_mnames) == 2):
        these_mnames, these_NM = zip(*sorted(zip(these_mnames, these_NM)))
        if these_mnames not in twoalign.keys():
          twoalign[these_mnames] = []
        twoalign[these_mnames].append((lasttagseq, these_NM))
        count_two += 1
        if count_two % 5000 == 0:
          print("{} tags aligning twice".format(count_two))
      # reset variables for new tag
      lasttagseq = tagseq
      these_mnames = [mname]
      these_NM = [NM]

# write the last marker to dict if necessary
# put the last tag into single alignment dict if applicable
if(len(these_mnames) == 1):
  mname_out = these_mnames[0]
  if mname_out not in onealign.keys():
    onealign[mname_out] = []
  onealign[mname_out].append((lasttagseq, these_NM[0]))
  count_one += 1
  if count_one % 500 == 0:
    print("{} tags aligning once".format(count_one))
# put the last tag into double alignment dict if applicable
if(len(these_mnames) == 2):
  these_mnames, these_NM = zip(*sorted(zip(these_mnames, these_NM)))
  if these_mnames not in twoalign.keys():
    twoalign[these_mnames] = []
  twoalign[these_mnames].append((lasttagseq, these_NM))
  count_two += 1
  if count_two % 500 == 0:
    print("{} tags aligning twice".format(count_two))

# filter to polymorphic sites
onealign = {k:v for k, v in onealign.items() if len(v) > 1}
twoalign = {k:v for k, v in twoalign.items() if len(v) > 1}

# write to files
with open(onealign_file, mode = "w", newline = '') as outcon:
  mywriter = csv.writer(outcon)
  for mrkr in sorted(onealign.keys()):
    for tagtup in onealign[mrkr]:
      mywriter.writerow([mrkr] + list(tagtup))

with open(twoalign_file, mode = "w", newline = '') as outcon:
  mywriter = csv.writer(outcon)
  for mrkrs in sorted(twoalign.keys()):
    for tagtup in twoalign[mrkrs]:
      mywriter.writerow(list(mrkrs) + [tagtup[0]] + list(tagtup[1]))
