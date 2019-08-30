# Thoughts on a general protocol for sorting paralogs

Each tag is aligned to the reference using Bowtie2, but instead of only returning
the top hit, two (or more) top hits are returned.  Many paralogs are going to be
highly similar or indistiguishable between subgenomes, so we need to accept that
many tags are not properly aligned.  Additionally, a difference between the
two subgenomes in the reference may not be fixed in the population, causing
tags without that variant to align to the wrong subgenome.

Tags are grouped based on the unordered set of top hits.  Tags that align more
times than expected are discarded, at least for now.  For tags that align once,
we can look at Hind/He to see if it makes sense to try to split them.

Preliminary grouping of tags into isoloci is based on which alignment to the
reference has fewer apparent mutations.  If a tag aligns equally well to two
locations, it is assigned to one at random.

Negative associations similar to those used in polysat are estimated using
Kendall's Tau on depth ratios.  To ensure that depth ratios are not negatively
associated between all pairs of alleles, the depth of the other allele in
question is omitted from the denominator of the depth ratio of a given allele.
A p-value threshold is lowered until groups are identified that don't violate
Hind/HE expectations.

An algorithm like simulated annealing or MCMC is used to shuffle tags between the
two or more isoloci.  The goal is to minimize Hind/HE for both/all isoloci.
The chance of moving a tag from one isolocus to another is influenced somewhat
by similarity to the reference.  After Hind/HE is below the expected value for
all isoloci, the algorithm tries to find solutions where haplotypes are most
similar to the reference for their assigned isolocus.  If there were groups
identified with Kendall's Tau, the algorithm tries to keep them together.

Perhaps this is all done in a Python script before import to polyRAD.  It will
necessitate looping through loci, which will be slow as molasses in R unless
Rcpp is used very heavily.

If one tag belongs to multiple isoloci, either the whole marker has to be
discarded, or polyRAD needs better models for dealing with allopolyploidy.
This will certainly happen, so there should be mechanisms in place for
recognizing it.

Will need to cite polyDog (doi: 10.1186/1471-2156-16-S2-S4)
