# Notes on Hind/He in a mapping population

In the diversity panel version, we have He * (1-F), which is the probability
that two alleles drawn from an individual without replacement are different
from each other.  This is mulitiplied by (ploidy - 1)/ploidy to get the
probability of drawing two different alleles with replacement.

In a mapping population, we still need the probability that two alleles in
a progeny drawn without replacement are different from each other.  In a diploid
F1, this is the same as the probability that one allele drawn from parent 1 and
one allele drawn from parent 2 are different from each other.

In a polyploid, do we also need the probability that two alleles drawn from one
parent are different from each other?

Cross       | Prob. parents diff | Prob alleles in progeny diff
-------------------------------------------------------------------------
AB x AB     |                0.5 |                          0.5
AABB x AABB |                0.5 |                          5/9

1/4 AA * 0
1/2 AB * 1
1/4 AA * 0

1/36  AAAA * 0
8/36  AAAB * 1/2
18/36 AABB * 2/3
8/36  ABBB * 1/2
1/36  BBBB * 0

= 4/36 + 12/36 + 4/36 = 5/9
