#!/usr/bin/env python3

# Process a SAM file of alignment of M. sacchariflorus tags to the M. sinensis
# reference.  Output tags that align once, and tags that align twice.  Include
# NM, the number of mutational steps from the reference for each alignment
# location.
# Bowtie2 was run with -k 8 to generate the SAM file.

import csv
import re
import argparse
import bisect

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
parser.add_argument("--min_ind_with_reads", "-m", nargs = 1, type = int,
                    default = 100,
                    help = "Minimum number of samples with sequencing reads needed to retain a tag set.")
args = parser.parse_args()
mysam = args.sam
myttd = args.ttd
outbase = args.out
maxalign = args.subgenomes
nchunks = args.chunks
min_ind_with_reads = args.min_ind_with_reads

#mysam = "D:/TASSELGBS_Msa/190517aligned_tags_multi.sam"
# dictionary to store alignments
aligndict = dict()
count = 0

# lists to store marker data for the current tag
these_mnames = [] # marker names
these_NM = [] # number of mutational steps between tag and reference

# subroutine for updating the alignment dictionary
def update_aligndict(these_mnames, these_NM, lasttagseq):
  these_mnames, these_NM = zip(*sorted(zip(these_mnames, these_NM)))
  if these_mnames not in aligndict.keys():
    aligndict[these_mnames] = []
  aligndict[these_mnames].append((lasttagseq, these_NM))

# read the sam file and process to aligndict
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
      # put the last tag into alignment dict if applicable
      if(len(these_mnames) <= maxalign):
        update_aligndict(these_mnames, these_NM, lasttagseq)
        count += 1
        if count % 5000 == 0:
          print("{} tags aligning".format(count))
      # reset variables for new tag
      lasttagseq = tagseq
      these_mnames = [mname]
      these_NM = [NM]

# write the last marker to dict if necessary
if(len(these_mnames) <= maxalign):
  update_aligndict(these_mnames, these_NM, lasttagseq)
  count += 1
print("{} tags aligning".format(count))

# filter to polymorphic sites
aligndict = {k:v for k, v in aligndict.items() if len(v) > 1}

## Part two: reading TagTaxaDist, filtering, and export

def extractTTD(tags, ttdfile):
  "Extract a given set of tags from the TTD file"
  # set up output matrix for TagTaxaDist
  ttd_mat = [[tag] for tag in tags]
  # sort tags for binary search
  sorted_tags, tagindex = zip(*sorted(zip(tags, range(len(tags)))))
  with open(ttdfile, mode = 'r') as mycon:
    header = next(mycon).split()[1:]
    for line in mycon:
      row = line.split()
      tag = row[0]
      ti = bisect.bisect(sorted_tags, tag)
      if sorted_tags[ti] != tag:
        continue # skip if this is not a tag we wanted to keep
      ti2 = tagindex[ti]
      ttd_mat[ti2].extend([int(d) for d in row[1:]])
  return header, ttd_mat

def keepMarker(depths, min_ind_with_reads):
  "Determine whether to keep a marker, based on missing data rate."
  nind = len(depths[0])
  ind_with_reads = [any([d[i] > 0 for d in depths]) for i in range(nind)]
  return sum(ind_with_reads) >= min_ind_with_reads

# divide into chunks and process TTD
print("Sorting alignment locations and processing TagTaxaDist")
all_markers = sorted(aligndict.keys())
m_per_chunk = len(all_markers) % nchunks
for chnk in range(nchunks):
  endm = len(these_markers) if chnk == nchunks - 1 else (chnk + 1) * m_per_chunk
  these_markers = all_markers[chnk * m_per_chunk : endm]
  nm = len(these_markers)
  tag_table = [] # to hold tag info in spreadsheet-like format
  # to index rows that correspond to given markers, same order as these_markers
  table_row_per_marker = [[-1] for i in range(len(these_markers))]
  curr_row = 0
  # loop to build table of tag info
  for mi in range(nm):
    m = these_markers[mi]
    # pad out marker names to maximum
    npad = maxalign - len(m)
    m_exp = list(m) + ['' for i in range(npad)]
    # get tags and NM
    m_tags = [tup[0] for tup in aligndict[m]]
    m_NM = [list(tup[1]) + ['' for i in range(npad)] for tup in aligndict[m]]
    nt = len(m_tags) # number of tags for this marker
    assert len(m_NM) == nt
    # add to table
    tag_table.extend([m_exp + m_tags[i] + m_NM[i] for i in range(nt)])
    table_row_per_marker[mi] = range(curr_row, curr_row + nt)
    curr_row += nt
  # extract all tag sequences
  these_tags = [tt[maxalign] for tt in tag_table]
  # extract read depth for these tags
  ttdheader, this_ttd = extractTTD(these_tags, myttd)
  # determine which markers (rows) to keep
  markers_kept = [keepMarker([this_ttd[ti] for ti in table_row_per_marker[mi]],
                             min_ind_with_reads) for mi in range(nm)]
  rows_kept = [ti for ti in table_row_per_marker[mi] for mi in range(nm) \
               if markers_kept[mi]]
  # subset tag table and depth table
  tag_table = [tag_table[ti] for ti in rows_kept]
  this_ttd = [this_ttd[ti] for ti in rows_kept]
  # export to files

# write to files
# with open(onealign_file, mode = "w", newline = '') as outcon:
#   mywriter = csv.writer(outcon)
#   for mrkr in sorted(onealign.keys()):
#     for tagtup in onealign[mrkr]:
#       mywriter.writerow([mrkr] + list(tagtup))
#
# with open(twoalign_file, mode = "w", newline = '') as outcon:
#   mywriter = csv.writer(outcon)
#   for mrkrs in sorted(twoalign.keys()):
#     for tagtup in twoalign[mrkrs]:
#       mywriter.writerow(list(mrkrs) + [tagtup[0]] + list(tagtup[1]))
