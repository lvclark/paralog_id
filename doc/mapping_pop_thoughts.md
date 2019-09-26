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

What we need is:

prob that two alleles drawn come from parent 1 * prob of two alleles from parent 1 being the same +
same for parent 2 +
prob that two alleles drawn come from different parents * prob of different alleles from two parents

If we do this with replacement, so ploidy doesn't matter:

AA x BB: A = 0.25 * 1 + 0.25 * 0 + 0.5 * 0 = 0.25; B = 0.25; Hind = 1 - 0.25 - 0.25 = 0.5
AA x AB: A = 0.25 * 1 + 0.25 * 0.5 + 0.5 * 0.5 = 0.5; B = 0.25 * 0 + 0.25 * 0.5 + 0.5 * 0 = 0.125
  Hind = 1 - 0.5 - 0.125 = 0.375
  
Nah, probably still needs to be ploidy-aware.
Prob that a pair of alleles, with replacment, both came from parent 1, parent 2, or different
parents, depends on ploidy, backcrossing, and inbreeding.

With replacement
| Ploidy | Generation | Both al. from p. 1 | Both al. from p. 2 | Al from diff. parents |
-----------------------------------------------------------------------------------------
| 2x     |         F1 |                  0 |                  0 |                     1 |
| 2x     |         F2 |              0.25  |              0.25  |                  0.5  |
| 2x     |         F3 |              0.375 |              0.375 |                  0.25 |
| 2x     |        BC1 |              0.5   |              0     |                  0.5  |
| 4x     |         F1 |                1/6 |                1/6 |                   4/6 |
| 4x     |

Maybe recursive functions to get gamete proportions, then calc from there.