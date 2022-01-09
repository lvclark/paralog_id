library(VariantAnnotation)
library(polyRAD)
library(GenomicRanges)
library(ggplot2)
library(dplyr)

# How do SNPs vs. haplotypes affect variance in Hind/He estimate?

myvcf <- "large_datasets/180208Msa_filtered.vcf.bgz"

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

# tiff("191219snp_vs_hap.tiff", width = 6.5 * 300, height = 4 * 300, res = 300,
#      compression = "lzw")
ggplot(df, aes(x = HindHe, color = Marker)) +
  geom_density() +
  facet_wrap(~ Ploidy) +
  geom_vline(aes(xintercept = Val), data = expvals, lty = 2) +
  labs(color = "Marker type")
# dev.off()

#cairo_pdf("SuppFig1_SNPs_vs_haps.pdf", width = 3.35, height = 4.5)
tiff("SuppFig1_SNPs_vs_haps.tiff", width = 3.35 * 300, height = 4.5 * 300, res = 300,
     compression = "lzw")
ggplot(df, aes(x = HindHe, color = Marker)) +
  geom_density() +
  facet_wrap(~ Ploidy, nrow = 2) +
  geom_vline(aes(xintercept = Val), data = expvals, lty = 2) +
  labs(color = "Marker type", x = expression(H[ind] / H[E])) +
  theme(legend.position = "bottom")
dev.off()

# Distance to nearest gene ####
gff0 <- rtracklayer::import("~/Genomes/Miscanthus reference/Msinensis_497_v7.1.gene_exons.gff3.gz")
gff1 <- gff0[gff0$type == "gene"]

gff1

snpGR <- GRanges("Chr01", IRanges(mydata2$locTable$Pos, width = 1))

gndist <- distanceToNearest(snpGR, gff1, ignore.strand = TRUE)
gndist
hist(log1p(mcols(gndist)$distance))
median(mcols(gndist)$distance) # 873 bp
summary(mcols(gndist)$distance)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.0      0.0    872.5   9484.4   7045.8 597751.0 

gndistDF <- data.frame(Locus = rep(GetLoci(mydata2), times = 2),
                       Ploidy = rep(c("Diploid", "Tetraploid"), each = nLoci(mydata2)),
                       HindHe = c(hh2loc_2x, hh2loc_4x),
                       Distance_to_gene = rep(mcols(gndist)$distance, times = 2))
gndistDF$Distance_class <- NA_character_
gndistDF$Distance_class[gndistDF$Distance_to_gene == 0] <- "In gene"
gndistDF$Distance_class[gndistDF$Distance_to_gene > 0 & gndistDF$Distance_to_gene <= 5e3] <- "Within 5kb of gene"
gndistDF$Distance_class[gndistDF$Distance_to_gene > 5e3 & gndistDF$Distance_to_gene <= 30e3] <- "Within 30kb of gene"
gndistDF$Distance_class[gndistDF$Distance_to_gene > 30e3] <- "Further than 30kb from gene"
any(is.na(gndistDF$Distance_class))
table(gndistDF$Distance_class)
# Further than 30kb from gene                     In gene         Within 30kb of gene          Within 5kb of gene 
#                        2050                        8220                        4038                        6608 

gndistDF$Distance_class <- factor(gndistDF$Distance_class,
                                  levels = c("In gene", "Within 5kb of gene",
                                             "Within 30kb of gene", "Further than 30kb from gene"))

mean(gndistDF$HindHe <= 2, na.rm = TRUE) # 100%

gndistDF %>%
  ggplot(aes(x = HindHe, fill = Distance_class)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Ploidy)

gndistDF %>%
  ggplot(aes(x = HindHe, fill = Distance_to_gene == 0)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Ploidy)
