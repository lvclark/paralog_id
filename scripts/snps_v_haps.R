library(VariantAnnotation)
library(polyRAD)

# How do SNPs vs. haplotypes affect variance in Hind/He estimate?

myvcf <- "~/DOE Msa study\\Seq\\GBSv2_180110/180208Msa_filtered.vcf.bgz"

myparam <- ScanVcfParam(fixed = "ALT", info = NA, geno = "AD",
                        which = GRanges("01", IRanges(1, 1e8)))

mydata1 <- VCF2RADdata(myvcf, svparam = myparam, phaseSNPs = TRUE, yieldSize = NA_integer_)

mydata2 <- VCF2RADdata(myvcf, svparam = myparam, phaseSNPs = FALSE, yieldSize = NA_integer_)

diploids <- readLines("marker_CSV/diploids.txt")
tetraploids <- readLines("marker_CSV/tetraploids.txt")

hh1_2x <- HindHe(SubsetByTaxon(mydata1, diploids))
hh2_2x <- HindHe(SubsetByTaxon(mydata2, diploids))

hh1_4x <- HindHe(SubsetByTaxon(mydata1, tetraploids))
hh2_4x <- HindHe(SubsetByTaxon(mydata2, tetraploids))

hh1loc_2x <- colMeans(hh1_2x, na.rm = TRUE)
hh2loc_2x <- colMeans(hh2_2x, na.rm = TRUE)
hh1loc_4x <- colMeans(hh1_4x, na.rm = TRUE)
hh2loc_4x <- colMeans(hh2_4x, na.rm = TRUE)

hist(hh1loc, col = "lightgrey")
hist(hh2loc, col = "lightgrey")

median(hh1loc_2x, na.rm = TRUE) # 0.372
median(hh2loc_2x, na.rm = TRUE) # 0.398

median(hh1loc_4x, na.rm = TRUE) # 0.535
median(hh2loc_4x, na.rm = TRUE) # 0.539

var(hh1loc_2x, na.rm = TRUE) # 0.0624
var(hh2loc_2x, na.rm = TRUE) # 0.0759

var(hh1loc_4x, na.rm = TRUE) # 0.0359
var(hh2loc_4x, na.rm = TRUE) # 0.0424

tab1 <- table(mydata1$alleles2loc)

plot(as.vector(tab1), hh1loc_2x)

plot(tapply(hh1loc_2x, tab1, var))

mean(hh1loc)
mean(hh2loc)

hist(rowMeans(hh1_2x, na.rm = TRUE), breaks = 30)
hist(rowMeans(hh2_2x, na.rm = TRUE), breaks = 30)

var(as.vector(hh1_2x), na.rm = TRUE) # 10.76
var(as.vector(hh2_2x), na.rm = TRUE) # 19.90

# identify regions with collapsed paralogs?
library(ggplot2)

ggplot(mapping = aes(x = mydata1$locTable$Pos, y = hh1loc_2x)) +
  geom_point() +
  geom_smooth()

# density plots of Hind/He for manuscript
df <- data.frame(HindHe = c(hh1loc_2x, hh2loc_2x, hh1loc_4x, hh2loc_4x),
                 Ploidy = rep(c("Diploids", "Tetraploids"),
                              times = c(length(hh1loc_2x) + length(hh2loc_2x),
                                        length(hh1loc_4x) + length(hh2loc_4x))),
                 Marker = c("Haplotype", "SNP")[rep(c(1,2,1,2),
                                                    times = c(length(hh1loc_2x), length(hh2loc_2x),
                                                              length(hh1loc_4x), length(hh2loc_4x)))])
expvals <- data.frame(Ploidy = c("Diploids", "Tetraploids"),
                      Val = c(1/2, 3/4))

tiff("191219snp_vs_hap.tiff", width = 6.5 * 300, height = 4 * 300, res = 300,
     compression = "lzw")
ggplot(df, aes(x = HindHe, color = Marker)) +
  geom_density() +
  facet_wrap(~ Ploidy) +
  geom_vline(aes(xintercept = Val), data = expvals, lty = 2) +
  labs(color = "Marker type")
dev.off()