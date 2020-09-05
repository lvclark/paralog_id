# Estimate variance and bias in Hind/He
library(polyRAD) # v1.3, Sept. 5 2020
library(ggplot2)

# function to get an estimate of hind/he and its variance, given biological and technical parameters.
# MAF = minor allele frequency.  Biallelic loci only.
# nsam = number of samples
# depth = read depth (one value; will be uniform across samples and loci)
# overdispersion = overdispersion parameter
# inbreeding = inbreeding coefficient
HHByParam <- function(MAF = 0.05, nsam = 500, depth = 20, overdispersion = 20,
                      inbreeding = 0, ploidy = 2, nloc = 5000){
  alleleFreq <- rep(c(1 - MAF, MAF), times = nloc)
  alleles2loc <- rep(1:nloc, each = 2)
  
  geno <- SimGenotypes(alleleFreq, alleles2loc, nsam, inbreeding, ploidy)
  alleleDepth <- SimAlleleDepth(matrix(depth, nrow = nsam, ncol = nloc,
                                       dimnames = list(NULL, as.character(1:nloc))),
                                geno, alleles2loc, overdispersion)
  rownames(alleleDepth) <- paste0("sam", 1:nsam)
  simrad <- RADdata(alleleDepth, alleles2loc,
                    locTable = data.frame(row.names = paste0("loc", 1:nloc)),
                    possiblePloidies = list(as.integer(ploidy)),
                    contamRate = 0.001,
                    alleleNucleotides = rep(c("A", "G"), times = nloc))
  
  hh <- HindHe(simrad)
  hhByLoc <- colMeans(hh, na.rm = TRUE)
  return(list(mean = mean(hhByLoc, na.rm = TRUE),
              var = var(hhByLoc, na.rm = TRUE),
              nonNA = mean(!is.na(hhByLoc)),
              values = hhByLoc))
}

test <- HHByParam()
hist(test$values)
sqrt(test$var)

# run tests for variance
ploidies <- c(2L, 4L)
freqs <- c(0.01, 0.05, 0.1, 0.2)
samsizes <- c(100, 500, 1000)
depths <- c(2, 5, 10, 20, 50, 100, 200)

tottests <- length(ploidies) * length(freqs) * length(samsizes) * length(depths)

testres <- data.frame(Ploidy = integer(tottests),
                      MAF = numeric(tottests),
                      N_sam = numeric(tottests),
                      Depth = numeric(tottests),
                      Mean = numeric(tottests),
                      Variance = numeric(tottests),
                      Prop_estimated = numeric(tottests))
currrow <- 1

for(p in ploidies){
  for(f in freqs){
    for(s in samsizes){
      for(d in depths){
        res <- HHByParam(MAF = f, nsam = s, depth = d, ploidy = p,
                         overdispersion = 20, inbreeding = 0, nloc = 5000)
        testres$Ploidy[currrow] <- p
        testres$MAF[currrow] <- f
        testres$N_sam[currrow] <- s
        testres$Depth[currrow] <- d
        testres$Mean[currrow] <- res$mean
        testres$Variance[currrow] <- res$var
        testres$Prop_estimated[currrow] <- res$nonNA
        currrow <- currrow + 1
        if(currrow %% 10 == 0) print(currrow)
      }
    }
  }
}

#save(testres, file = "workspaces/variance_estimates.RData")

testres$PloidyText <- ifelse(testres$Ploidy == 2, "Diploid", "Tetraploid")

# plot results
ggplot(testres[testres$Depth < 200,],
       aes(x = Depth, y = sqrt(Variance), color = as.factor(MAF), group = MAF)) +
  geom_line() +
  facet_grid(PloidyText ~ N_sam, labeller = labeller(N_sam = function(x) paste("N =", x))) +
  scale_x_continuous(breaks = depths, trans = "log2") +
  labs(y = "Standard deviation of estimate", color = "MAF", x = "Read depth")

# bias -- overestimated for rare alleles and small sample size, since so many NA
ggplot(testres[testres$Depth < 200,],
       aes(x = Depth, y = Mean, color = as.factor(MAF), group = MAF)) +
  geom_line() +
  facet_grid(PloidyText ~ N_sam, labeller = labeller(N_sam = function(x) paste("N =", x))) +
  scale_x_continuous(breaks = depths, trans = "log2") +
  labs(y = "Mean estimate", color = "MAF", x = "Read depth")
