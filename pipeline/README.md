# Pipeline for sorting tags into isoloci

This pipeline is intended for genotype calling from reduced-representation
DNA sequencing data, such as GBS or RADseq, where the reference genome contains
one or more whole genome duplications, for example in an allopolyploid organism.
Most ancestral loci correspond to two or more isoloci (paralogous loci).
Because conventional alignment software may not align each read to the correct
paralog, this pipeline uses the distribution of read depth in individuals and
populations, along with alignment results, to assign sequence tags to the
correct isolocus, or filter loci out that cannot be adequately sorted.

Python 3 is required for the Python scripts.

## Step 1: TASSEL-GBSv2 pipeline

The [TASSEL-GBSv2](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline)
pipeline is used to find all unique sequence tags in a set of FASTQ files, and
tally the depth of every tag in every individual in the dataset.

Run `GBSSeqtoTagDBPlugin` and `TagExportToFastqPlugin` as normal. Additionally,
run `GetTagTaxaDistFromDBPlugin` to get a tab-delimited text file with the
depth of each tag in each individual.

## Step 2: Sequence alignment

Under default parameters, Bowtie2 and BWA only return one alignment for each
sequence read.  However, for this pipeline we want multiple alignments for
each read, under the assumption that the best alignment may not be correct,
or there may be multiple best alignments with equal alignment scores.

With [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), align
the FASTQ output by `TagExportToFastqPlugin` to your reference genome using
`-k` mode.  Because `-k` mode does not guarantee that the best alignments are
reported, but rather reports the first `k` alignments that it finds, I recommend
setting `k` to something higher than the number of anticipated isoloci.  For
example, in an allotetraploid, I expect two isoloci, but I need to set `-k 3`
or higher to confirm that there are not additional valid alignments.

The command might look like:

```
bowtie2 --very-sensitive -k 3 -x my_reference -U tagexport.fq -S myalign.sam
```

## Step 3: Grouping tags by sets of alignments

The `process_sam_multi.py` script parses the TagTaxaDist file from Step 1 and
the SAM file from Step 2.  Tags that aligned more times than the expected
number of isoloci are discarded.  The remaining tags are placed into groups
based on alignment location or unordered set of alignment locations.  This step
also allows filtering of loci based on a minimum number of individuals with
reads (*i.e.* filtering based on the missing data rate), as well as subsetting
of the read depth matrix by individual.  The dataset may optionally be split
into chunks for parallel processing and/or memory conservation downstream.
For each chunk, two files are generated: (1) an alignment file, containing the
tag sequence, alignment locations up to the maximum allowed based on the number
of expected isoloci, and the edit distance (`NM`) for each alignment, and (2) a
subsetted `TagTaxaDist` file with read depth for the same tags.

Required positional arguments:

* `sam`: Path to the SAM file output by Step 2.
* `ttd`: Path to the `TagTaxaDist` file output by Step 1.
* `out`: Base file name for output.

Optional, named arguments:

* `--subgenomes`: The number of expected isoloci for each group of tags, *i.e.*
the number of subgenomes within an allopolyploid.  Default 2.
* `--chunks`: Number of chunks into which to split the output.  Default 1.
* `--min_ind_with_reads`: For a given group of tags aligning to the same
location or set of locations, this many individuals must have more than zero
reads for the tag group to be retained for further analysis.  Default 100.
* `--samples`: Path to a file containing names of samples to retain, matching
the sample names that went into TASSEL.  The file should simply contain one
sample name per line and no other data.  If your dataset is a mixture of
ploidies, only process samples for one ploidy at a time.

Example use:

``` bash
python process_sam_multi.py myalign.sam mytagtaxadist.txt myout --chunks 10
```

## Step 4: Sorting alleles into isoloci

(Note that dollar signs indicate LaTeX notation; this may not render in GitHub
but will render in other Markdown previewers such as RStudio.)

The `process_isoloci.py` script takes the output from Step 3 (one chunk at a
time) and attempts to sort groups of tags into the appropriate number of
isoloci, based on how many alignments were reported.

A statistic called $H_{ind}/H_E$ (Hind/He) is used to determine whether a group
of tags is behaving like a normal locus, or whether it likely contains tags
from multiple isoloci.  $H_{ind}$ is the probability of sampling two different
tags, with replacement, from all the reads in that tag group for a given
individual, averaged across individuals.  $H_E$ is the expected heterozygosity,
based on estimated allele frequencies, assuming all tags belong to one locus.
The expected value for $H_{ind}/H_E$ is
$\frac{ploidy - 1}{ploidy} * F$, where $F$ is the inbreeding coefficient.  Higher
values indicate the presence of paralogous tags in the group.  The algorithm
strives to get $H_{ind}/H_E$ at or below the expected value for each isolocus,
and discards isoloci that exceed a higher threshold (halfway between the
expected value for $ploidy$ and $2 * ploidy$, on a log scale).

Initial assignments are made based on `NM` values; tags are assigned to
whichever isolocus they better aligned, or a random isolocus chosen among
equally good alignments.  Negative allele associations are then detected
similarly to the method of
[Clark and Schreier (2017)](https://doi.org/10.1111/1755-0998.12639), and
allele assignments are adjusted if necessary to make sure that sets of
negatively associated alleles are assigned to the same isolocus.  Lastly, if
any isolocus exceeds the expected value of $H_{ind}/H_E$, a
tabu search is performed, moving alleles to different isoloci for 25 cycles.
The best solution found in the tabu search is then returned as the optimal
allele assignment.  The 'best' solution is considered to be the one that
exceeds the expected $H_{ind}/H_E$ by the minimum amount, or in the case of
a tie, has the fewest number of mutations between alleles and the reference
genome.

The script outputs a single CSV file containing alignment information for each
tag (corrected according to the above algorithm) as well as read depth for every
individual.

Required positional arguments:

* `alignfile`: Path to alignment file output by Step 3.
* `depthfile`: Path to depth file output by Step 3.

Optional, named arguments:

* `--out`: File name for output.  Generated from input files if not provided.
* `--ploidy`: Expected ploidy after splitting isoloci.  Default 2.
* `--inbreeding`: Inbreeding coefficient $F$.  Default 0.  Can range from 0 to 1.
* `--logfile`: Path to file for writing log.  By default, no log will be written.
For every set of alignment positions, the log contains information about initial
and final $H_{ind}/H_E$ and average `NM` values.

For the script to run, `isoloci_fun.py` must be in the same directory as
`process_isoloci.py`.

Example use:

``` bash
python process_isoloci.py myout_4_align.csv myout_4_depth.csv --inbreeding 0.2
```

## Step 5: Genotype calling with polyRAD

I still need to make an import function, but the file output by Step 4 contains
all information required by polyRAD.
