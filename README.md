# paralog_id
Exploration and testing of methods to identify paralogs in GBS data

This repository includes scripts, data, and drafts pertaining to methods for
sorting GBS/RAD tags into paralogous loci.  The dataset being used is
from *Miscanthus sacchariflorus*, which contains both diploid and tetraploid
individuals.  Groups of tags were identified that aligned to one location
in the *Sorghum bicolor* v3 reference genome, but two locations in the
*Miscanthus sinensis* v7 reference genome, where those two locations are
on appropriate chromosomes given the known synteny between *S. bicolor* and
*M. sinensis*.

## File descriptions

### CSV files

`marker_CSV/190513miscanthus_sorghum_match.csv` contains TagDigger output listing every
tag location from the *M. sinensis* reference and how it matches, if at all,
to the *S. bicolor* reference.  SAM files for TagDigger were generated by
Bowtie2 as part of the TASSEL-GBSv2 pipeline and are available at the
[Illinois Data Bank](https://doi.org/). **Need to add DOI once published**

`marker_CSV/190513paralogs.csv` is filtered from the above file, to only contain pairs
of *M. sinensis* tag locations that align to the same *S. bicolor* location,
and only with appropriate synteny.

`marker_CSV/190515paralog_tags.csv` contains tag sequences for all tag locations listed
in the above file.

### Code

`scripts/hindhe_by_ploidy.R` imports *M. sacchariflorus* read depths from a VCF and
plots the distribution of Hind/He by ploidy.  This is an older figure used in
presentations but not included in the manuscript.

`scripts/filter_markers.R` imports `marker_CSV/190513miscanthus_sorghum_match.csv`
and filters it down to loci that are clearly paralogous, outputting
`marker_CSV/190513paralogs.csv`.

`scripts/get_tag_seq.R` imports `marker_CSV/190513paralogs.csv` and a large
SAM file containing alignments of *M. sacchariflorus* tags to the
*M. sinensis* reference, then exports the sequences of tags for all
selected markers to `marker_CSV/190515paralog_tags.csv`.

`scripts/import_tagtaxadist.R` imports `marker_CSV/190515paralog_tags.csv` as
well as a list of ploidy by accession and a large TagTaxaDist file containing
read counts output by TASSEL.  It imports read counts just for the desired
tags and just in diploid and tetraploid individuals, then exports the read
depth matrices to `workspaces/190515counts_matrices.RData`.

`scripts/H_stats_sorghum_miscanthus.R` imports `workspaces/190515counts_matrices.RData`
and `marker_CSV/190515paralog_tags.csv`, then 
estimates $H_{ind}$ and $H_E$ for all markers.
The goal is to demonstrate that Hind/He can clearly distinguish one-copy and
two-copy loci, *i.e.* those aligned to *Miscanthus* vs. those aligned to
*Sorghum*.  Figures are plotted to compare the distribution of Hind/He
when tags are aligned to *Miscanthus* vs. *Sorghum*.

`scripts/optimize_temperature.R` examines the output of `process_isoloci.py` to
evaluate the effectiveness of the tabu search algorithm for optimizing Hind/He.
Figures are plotted to compare Hind/He before and after the tabu search.

`scripts/get_inbreeding.R` tests code to estimate inbreeding from preliminary
Hind/He distributions, and also plots Hind/He vs. ploidy, depth, and proportion
_M. sinensis_ ancestry for the manuscript.

`scripts/isoloci_fun_test.py` contains code for testing individual functions
within the variant calling pipeline.

`scripts/snps_v_haps.R` imports variants from a VCF, with and without phasing
SNPs into haplotypes.  The distribution and variance of Hind/He is then compared
between SNPs and haplotypes.  A figure is generated for visualizing the
distribution.

`scripts/variance_and_bias.R` uses simulated data to explore the impact of
population and techinical parameters on the variance and bias of Hind/He.

### Documentation

`doc/using_hindhe.Rmd` contains a brief exploration of the Hind/He statistic
using an example dataset provided with polyRAD.

`doc/proof_heterozygosity_for_isolocus_detection.docx` is a mathematical proof
demonstrating that the average value for $H_{ind}$ is always expected to be
$\frac{ploidy - 1}{ploidy}H_E$.  It also describes the algorithm implemented
in the pipeline.

`doc/next_steps.md` describes additional work not done yet.
