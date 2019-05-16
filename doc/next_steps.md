# Thoughts on a general protocol for sorting paralogs

Each tag is aligned to the reference using Bowtie2, but instead of only returning
the top hit, two (or more) top hits are returned.  Many paralogs are going to be
highly similar or indistiguishable between subgenomes, so we need to accept that
many tags are not properly aligned.

Tags are grouped based on the unordered set of top hits.

Similarly to how polysat looks for negative associations between the presence
and absence of alleles in order to group them into isoloci, I would like to look
for negative associations between read depth ratios (as estimated in polyRAD when
it builds the RADdata object) to make a preliminary grouping of tags into isoloci.

The preliminary grouping should also be influenced by sequence similarity.  How
much should depend on whether allelic associations or sequence similarity gives
cleaner results.  Maybe this is similar to Bayesian stats where posterior probs
are influenced sometimes more by the priors and sometimes more by the likelihoods.

An algorithm like simulated annealing or MCMC is used to shuffle tags between the
two or more isoloci.  The goal is to minimize Hind/HE for both/all isoloci.

Perhaps this is all done in a Python script before import to polyRAD.  It will
necessitate looping through loci, which will be slow as molasses in R unless
Rcpp is used very heavily.
