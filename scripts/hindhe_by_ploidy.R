library(polyRAD)
library(ggplot2)

# Look at mix of diploids and tetraploids to see if distinguishable ####
library(VariantAnnotation)
myRAD <- VCF2RADdata("large_datasets/180208Msa_filtered.vcf.bgz",
                     svparam = ScanVcfParam(fixed = "ALT", info = NA, geno = "AD",
                                            which = GRanges("05", IRanges(1, 112e6))),
                     yieldSize = NA)

myHindHe <- HindHe(myRAD)

TotDepthT <- rowSums(myRAD$locDepth)

myHindHeByInd <- rowMeans(myHindHe, na.rm = TRUE)

accessions <- read.csv("marker_csv/all_accession_names.csv",
                       stringsAsFactors = FALSE)
ploidies <- accessions$Ploidy[match(names(myHindHeByInd), accessions$Accession)]

ggplot(mapping = aes(x = TotDepthT, y = myHindHeByInd, color = ploidies)) +
  geom_point() +
  scale_x_log10()

hist(myHindHeByInd, breaks = 50)
# bimodal distribution, not perfect but could be very helpful

tiff("MsaPloidyID.tiff", width = 1800, height = 1200, res = 300, compression = "lzw",
     pointsize = 16)
ggplot(data.frame(Depth = TotDepthT, HindHe = myHindHeByInd, Ploidy = ploidies)[ploidies %in% c("2x", "3x", "4x"),],
       mapping = aes(x = Depth, y = HindHe, color = Ploidy)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~ Ploidy) +
  geom_hline(data = data.frame(Ploidy = c("2x", "3x", "4x"),
                               ExpHindHe = c(1/2, 2/3, 3/4)),
             mapping = aes(yintercept = ExpHindHe), lty = 2)
dev.off()
