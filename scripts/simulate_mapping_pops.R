library(polyRAD)
library(ggplot2)
library(dplyr)

# cross types to simulate
# values are allele copy numbers for parents 1 and 2

cross_types <-
  list(Diploid = list(Testcross = c(1,0),
                       F2 = c(1, 1)),
        Tetraploid = list(Simplex_x_Nulliplex = c(1, 0),
                          Duplex_x_Nulliplex = c(2, 0),
                          Triplex_x_Nulliplex = c(3, 0),
                          Simplex_x_Simplex = c(1, 1),
                          Duplex_x_Simplex = c(2, 1),
                          Triplex_x_Simplex = c(3, 1),
                          Duplex_x_Duplex = c(2, 2)))

# function to get an estimate of hind/he and its variance, given biological and technical parameters.
# Cross type, in terms of the genotype of parent 1 and parent 2
# nsam = number of samples
# depth = read depth (one value; will be uniform across samples and loci)
# overdispersion = overdispersion parameter
HHByParamMapping <- function(cross_type = c(1, 0), nsam = 500, depth = 20, overdispersion = 20,
                      ploidy = 2, nloc = 5000, n.gen.backcrossing = 0, n.gen.selfing = 0){
  alleles2loc <- rep(1:nloc, each = 2)
  parent1 <- rep(c(cross_type[1], ploidy - cross_type[1]), times = nloc)
  parent2 <- rep(c(cross_type[2], ploidy - cross_type[2]), times = nloc)
  
  geno <- SimGenotypesMapping(parent1, parent2, alleles2loc, nsam, ploidy,
                              n.gen.backcrossing, n.gen.selfing)
  geno <- rbind(parent1, parent2, geno)
  alleleDepth <- SimAlleleDepth(matrix(depth, nrow = nsam + 2, ncol = nloc,
                                       dimnames = list(NULL, as.character(1:nloc))),
                                geno, alleles2loc, overdispersion,
                                contamRate = 0, errorRate = 0.001)
  rownames(alleleDepth) <- c("parent1", "parent2", paste0("sam", 1:nsam))
  simrad <- RADdata(alleleDepth, alleles2loc,
                    locTable = data.frame(row.names = paste0("loc", 1:nloc)),
                    possiblePloidies = list(as.integer(ploidy)),
                    contamRate = 0.001,
                    alleleNucleotides = rep(c("A", "G"), times = nloc))
  simrad <- SetDonorParent(simrad, "parent1")
  simrad <- SetRecurrentParent(simrad, "parent2")
  
  hh <- HindHeMapping(simrad, n.gen.backcrossing, n.gen.selfing)
  hhByLoc <- colMeans(hh, na.rm = TRUE)
  return(list(mean = mean(hhByLoc, na.rm = TRUE),
              var = var(hhByLoc, na.rm = TRUE),
              nonNA = mean(!is.na(hhByLoc)),
              values = hhByLoc))
}

test <- HHByParamMapping()
test
hist(test$values)

# set up data frames for output
tottests5 <- sum(lengths(cross_types))
nloc5 <- 5000L

testres5 <- data.frame(Ploidy = rep(names(cross_types), times = lengths(cross_types)),
                       Cross_type = unlist(lapply(cross_types, names)),
                       Mean = numeric(tottests5),
                       Variance = numeric(tottests5),
                       Prop_estimated = numeric(tottests5))

locres5 <- data.frame(Ploidy = rep(testres5$Ploidy, each = nloc5),
                      Cross_type = rep(testres5$Cross_type, each = nloc5),
                      Estimate = rep(NA_real_, tottests5 * nloc5))

for(i in seq_len(tottests5)){
  pld <- testres5$Ploidy[i]
  ct <- testres5$Cross_type[i]
  if(pld == "Diploid") p <- 2L
  if(pld == "Tetraploid") p <- 4L
  
  res <- HHByParamMapping(cross_types[[pld]][[ct]], nsam = 500, depth = 20,
                          overdispersion = 20, ploidy = p, nloc = nloc5, 
                          n.gen.backcrossing = 0L, n.gen.selfing = 0L)
  testres5$Mean[i] <- res$mean
  testres5$Variance[i] <- res$var
  testres5$Prop_estimated[i] <- res$nonNA
  
  theserows <- which(locres5$Ploidy == pld & locres5$Cross_type == ct)
  theserows <- head(theserows, length(res$values))
  locres5$Estimate[theserows] <- res$values
}

#save(testres5, locres5, file = "workspaces/sim_mapping_pops_2022-01-29.RData")

locres5$Cross_type <- factor(locres5$Cross_type, levels = unique(locres5$Cross_type))

locres5 %>% filter(Ploidy == "Diploid") %>%
ggplot(mapping = aes(x = Estimate)) +
  geom_histogram() +
  facet_wrap(~ Cross_type)

locres5 %>% filter(Ploidy == "Tetraploid") %>%
  ggplot(mapping = aes(x = Estimate)) +
  geom_histogram() +
  facet_wrap(~ Cross_type)

ggplot(locres5, aes(x = Estimate, color = Cross_type)) +
  geom_density() +
  facet_wrap(~ Ploidy, scales = "free_x") +
  scale_color_manual(values = dittoSeq::dittoColors(1))

testres5
sqrt(testres5$Variance) # standard deviation of estimate

# Expected het for each cross type
locres5$Expected_heterozygosity <- numeric(nrow(locres5))
locres5$Expected_heterozygosity[locres5$Ploidy == "Diploid"] <- 0.5
for(crs in names(cross_types$Tetraploid)){
  parents <- cross_types$Tetraploid[[crs]]
  gamprob1 <- polyRAD:::.gameteProb(polyRAD:::.makeGametes(parents[1], 4), 4)
  gamprob2 <- polyRAD:::.gameteProb(polyRAD:::.makeGametes(parents[2], 4), 4)
  out <- 0
  for(i in 0:2){
    for(j in 0:2){
      progprob <- gamprob1[i+1,] * gamprob2[j+1,]
      geno <- i + j
      het <- (1 - (geno / 4) ^ 2 - ((4 - geno) / 4) ^ 2) * 0.75
      out <- out + progprob * het
    }
  }
  locres5$Expected_heterozygosity[locres5$Cross_type == crs] <- out
}


ggplot(locres5, aes(x = Estimate)) +
  geom_density() +
  facet_wrap(~ paste(Ploidy, Cross_type), scales = "free_x")

# cairo_pdf("Fig7_mapping.pdf", width = 6.7, height = 3.5)
# tiff("Fig7_mapping.tiff", res = 300, width = 6.7 * 300, height = 3.5 * 300,
#      compression = "lzw")
ggplot(locres5, aes(y = Estimate, x = Cross_type, fill = as.factor(Expected_heterozygosity))) +
  geom_violin(draw_quantiles = 0.5) +
  facet_grid(rows = vars(Ploidy) , scales = "free_y", space = "free_y") +
  coord_flip() +
  scale_fill_viridis_d() +
  labs(fill = "Expected heterozygosity", x = "Cross type")
# dev.off()
