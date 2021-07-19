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
                                geno, alleles2loc, overdispersion)
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

#save(testres5, locres5, file = "workspaces/sim_mapping_pops.RData")

locres5$Cross_type <- factor(locres5$Cross_type, levels = unique(locres5$Cross_type))

locres5 %>% filter(Ploidy == "Diploid") %>%
ggplot(mapping = aes(x = Estimate)) +
  geom_histogram() +
  facet_wrap(~ Cross_type)

locres5 %>% filter(Ploidy == "Tetraploid") %>%
  ggplot(mapping = aes(x = Estimate)) +
  geom_histogram() +
  facet_wrap(~ Cross_type)

#cairo_pdf("Fig5_mapping.pdf", width = 6.7, height = 3.5)
ggplot(locres5, aes(x = Estimate, color = Cross_type)) +
  geom_density() +
  facet_wrap(~ Ploidy, scales = "free_x") +
  scale_color_manual(values = dittoSeq::dittoColors(1))
#dev.off()

testres5
sqrt(testres5$Variance) # standard deviation of estimate
