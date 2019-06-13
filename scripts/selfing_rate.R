Rcpp::sourceCpp("C:/users/lvclark/RPackages/ploidyverseVcf/src/multiallele_utils.cpp")

# confirm that observed heterozygosity decreases at the same rate in a selfing tetraploid pop as selfing diploid pop

sm42 <- selfingMatrix(4, 2)
sm22 <- selfingMatrix(2, 2)
sm42
sm22

qmat <- function(ploidy){
  # heterozygosity by genotype, for two alleles only
  cp = 0:ploidy
  q <- cp/ploidy * (ploidy - cp)/(ploidy - 1) * 2
  return(matrix(q, nrow = ploidy + 1, ncol = 1))
}

# probability of sampling two different alleles (no replacement) from genotype
q2 <- qmat(2)
q4 <- qmat(4)

sm22 %*% q2 # goes down by half
sm42 %*% q4 # goes down by 5/6
sm42 %*% q4 / q4

sm62 <- selfingMatrix(6, 2)
q6 <- qmat(6)
q6

sm62 %*% q6
sm62 %*% q6 / q6 # goes down by 9/10

selfingMatrix(8, 2) %*% qmat(8) / qmat(8)

# so...
# 2 = 1/2
# 4 = 5/6
# 6 = 9/10
# 8 = 13/14
# (2 * ploidy - 3) / (2 * ploidy - 2)
