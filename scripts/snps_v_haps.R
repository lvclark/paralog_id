library(VariantAnnotation)
library(polyRAD)

# How do SNPs vs. haplotypes affect variance in Hind/He estimate?

myvcf <- "~/DOE Msi study\\Seq\\TASSEL-GBS\\PstI_GBSv2_170602/170608Msi_PstI_genotypes.vcf.bgz"

myparam <- ScanVcfParam(fixed = "ALT", info = NA, geno = "AD",
                        which = GRanges("01", IRanges(1, 1e8)))

mydata1 <- VCF2RADdata(myvcf, svparam = myparam, phaseSNPs = TRUE, yieldSize = NA_integer_)

mydata2 <- VCF2RADdata(myvcf, svparam = myparam, phaseSNPs = FALSE, yieldSize = NA_integer_)

hh1 <- HindHe(mydata1)
hh2 <- HindHe(mydata2)

hh1loc <- colMeans(hh1, na.rm = TRUE)
hh2loc <- colMeans(hh2, na.rm = TRUE)

hist(hh1loc, col = "lightgrey")
hist(hh2loc, col = "lightgrey")

median(hh1loc) # 0.263
median(hh2loc) # 0.276

var(hh1loc) # 0.0267
var(hh2loc) # 0.0370

tab1 <- table(mydata1$alleles2loc)

plot(as.vector(tab1), hh1loc)

plot(tapply(hh1loc, tab1, var))

mean(hh1loc)
mean(hh2loc)

hist(rowMeans(hh1, na.rm = TRUE), breaks = 30)
hist(rowMeans(hh2, na.rm = TRUE), breaks = 30)

var(as.vector(hh1), na.rm = TRUE) # 1.37
var(as.vector(hh2), na.rm = TRUE) # 3.58
