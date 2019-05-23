# Thoughts on a general protocol for sorting paralogs

He is now based on depth ratios rather than raw depth, so the Word document should
be updated.

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
locations, it is assigned to one at random.  Negative associations similar to
those used in polysat don't seem to be helpful.

An algorithm like simulated annealing or MCMC is used to shuffle tags between the
two or more isoloci.  The goal is to minimize Hind/HE for both/all isoloci.
The chance of moving a tag from one isolocus to another is influenced somewhat
by similarity to the reference.

Perhaps this is all done in a Python script before import to polyRAD.  It will
necessitate looping through loci, which will be slow as molasses in R unless
Rcpp is used very heavily.

If one tag belongs to multiple isoloci, either the whole marker has to be
discarded, or polyRAD needs better models for dealing with allopolyploidy.
This will certainly happen, so there should be mechanisms in place for
recognizing it.
