# Identify markers where we clearly have two paralogs in miscanthus for one
# marker location in sorghum

markers <- read.csv("190513miscanthus_sorghum_match.csv",
                    stringsAsFactors = FALSE)
head(markers)

markers <- markers[markers$Sbicolor_3 != "",]
markers <- markers[order(markers$Sbicolor_3),]

sorghum_nhits <- table(markers$Sbicolor_3)
head(sorghum_nhits)
twohits <- names(sorghum_nhits)[sorghum_nhits == 2]
  # --> 25,382 markers in sorghum match two locations in miscanthus

markers <- markers[markers$Sbicolor_3 %in% twohits,]

# chromosome matches that make sense
okmatches <- list(Chr01 = c("Chr01", "Chr02"),
                  Chr02 = c("Chr03", "Chr04"),
                  Chr03 = c("Chr05", "Chr06"),
                  Chr04 = c("Chr07", "Chr08"),
                  Chr05 = c("Chr09", "Chr10"),
                  Chr06 = c("Chr11", "Chr12"),
                  Chr07 = c("Chr07", "Chr13"),
                  Chr08 = c("Chr14", "Chr15"),
                  Chr09 = c("Chr16", "Chr17"),
                  Chr10 = c("Chr18", "Chr19"))

# check that the matches make sense
keep <- logical(length(twohits))
names(keep) <- twohits
for(m in twohits){
  sbchrom <- sub("-.*$", "", m)
  mismark <- markers$Query[markers$Sbicolor_3 == m]
  mischrom <- sub("-.*$", "", mismark)
  keep[m] <- all(okmatches[[sbchrom]] %in% mischrom)
}

sum(keep)

markers <- markers[markers$Sbicolor_3 %in% twohits[keep],]

head(markers)
markers[10001:10030,]
markers[20001:20030,]
markers[30001:30030,]
markers[35001:35030,]

write.csv(markers, file = "190513paralogs.csv", row.names = FALSE)
