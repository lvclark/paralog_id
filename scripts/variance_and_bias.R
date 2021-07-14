# Estimate variance and bias in Hind/He
library(polyRAD) # v1.3, Sept. 5 2020
library(ggplot2)
library(dplyr)

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
load("workspaces/variance_estimates.RData")

testres$PloidyText <- ifelse(testres$Ploidy == 2, "Diploid", "Tetraploid")

# plot results
ggplot(testres[testres$Depth < 200,],
       aes(x = Depth, y = sqrt(Variance), color = as.factor(MAF), group = MAF)) +
  geom_line() +
  facet_grid(PloidyText ~ N_sam, labeller = labeller(N_sam = function(x) paste("N =", x))) +
  scale_x_continuous(breaks = depths, trans = "log2") +
  labs(y = "Standard deviation of estimate", color = "MAF", x = "Read depth")
ggsave("~/NSF polyRAD/Year 4 report/fig1a.png",
       width = 6.5, height = 3.5)

# bias -- overestimated for rare alleles and small sample size, since so many NA
ggplot(testres[testres$Depth < 200,],
       aes(x = Depth, y = Mean, color = as.factor(MAF), group = MAF)) +
  geom_line() +
  facet_grid(PloidyText ~ N_sam, labeller = labeller(N_sam = function(x) paste("N =", x))) +
  scale_x_continuous(breaks = depths, trans = "log2") +
  labs(y = "Mean estimate", color = "MAF", x = "Read depth")
ggsave("~/NSF polyRAD/Year 4 report/fig1b.png",
       width = 6.5, height = 3.5)

# Explore overdispersion.  polyRAD 1.4 June 2021.
ods <- 5:20
freqs2 <- c(0.01, 0.05)
inbreed <- seq(0, 1, by = 0.1)

tottests2 <- length(ploidies) * length(ods) * length(freqs2) * length(inbreed)
testres2 <- data.frame(Ploidy = integer(tottests2),
                      MAF = numeric(tottests2),
                      Overdispersion = numeric(tottests2),
                      Inbreeding = numeric(tottests2),
                      Mean = numeric(tottests2),
                      Variance = numeric(tottests2),
                      Prop_estimated = numeric(tottests2))
currrow <- 1

for(p in ploidies){
  for(f in freqs2){
    for(o in ods){
      for(i in inbreed){
        res <- HHByParam(MAF = f, nsam = 500L, depth = 20L, ploidy = p,
                         overdispersion = o, inbreeding = i, nloc = 20000)
        testres2$Ploidy[currrow] <- p
        testres2$MAF[currrow] <- f
        testres2$Overdispersion[currrow] <- o
        testres2$Inbreeding[currrow] <- i
        testres2$Mean[currrow] <- res$mean
        testres2$Variance[currrow] <- res$var
        testres2$Prop_estimated[currrow] <- res$nonNA
        currrow <- currrow + 1L
        if(currrow %% 10 == 0) print(currrow)
      }
    }
  }
}

testres2$PloidyText <- ifelse(testres2$Ploidy == 2, "Diploid", "Tetraploid")

#save(testres2, file = "workspaces/variance_estimates_overdispersion_inbreeding.RData")
load("workspaces/variance_estimates_overdispersion_inbreeding.RData")

# Do this fig for hexaploids too since there is interest from sweetpotato
p <- 6L

tottests3 <- length(ods) * length(freqs2) * length(inbreed)
testres3 <- data.frame(Ploidy = 6L,
                       MAF = numeric(tottests3),
                       Overdispersion = numeric(tottests3),
                       Inbreeding = numeric(tottests3),
                       Mean = numeric(tottests3),
                       Variance = numeric(tottests3),
                       Prop_estimated = numeric(tottests3),
                       PloidyText = "Hexaploid")
currrow <- 1

for(f in freqs2){
  for(o in ods){
    for(i in inbreed){
      res <- HHByParam(MAF = f, nsam = 500L, depth = 20L, ploidy = p,
                       overdispersion = o, inbreeding = i, nloc = 20000)
      testres3$MAF[currrow] <- f
      testres3$Overdispersion[currrow] <- o
      testres3$Inbreeding[currrow] <- i
      testres3$Mean[currrow] <- res$mean
      testres3$Variance[currrow] <- res$var
      testres3$Prop_estimated[currrow] <- res$nonNA
      currrow <- currrow + 1L
      if(currrow %% 10 == 0) print(currrow)
    }
  }
}

#save(testres3, file = "workspaces/variance_estimates_overdispersion_inbreeding_hexaploid.RData")

testres2 <- rbind(testres2, testres3)
testres2$PloidyText <- factor(testres2$PloidyText, levels = c("Diploid", "Tetraploid", "Hexaploid"))

ggplot(testres2, aes(x = Overdispersion, y = Mean, color = Inbreeding, group = Inbreeding)) +
  geom_line() +
  facet_grid(PloidyText ~ MAF, scales = "free_y") +
  scale_color_viridis_c()

# version of figure for tutorial
png("overdispersion_inbreeding.png", width = 900, height = 600)
testres2 %>% filter(Overdispersion <= 20) %>%
ggplot(mapping = aes(x = Overdispersion, y = Mean, color = Inbreeding, group = paste(Inbreeding, MAF))) +
  geom_smooth(aes(lty = as.factor(MAF))) + ## Note that this is geom_smooth and not geom_line
  facet_grid(~ PloidyText, scales = "free_y") +
  scale_color_viridis_c() +
  labs(lty = "MAF", y = "Mean Hind/He") +
  geom_text(data = testres2[testres2$Overdispersion == 20 & testres2$MAF == 0.01,],
            mapping = aes(label = Inbreeding, y = Mean + 0.02),
            x = 20) +
  scale_x_continuous(minor_breaks = 5:20) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1))
dev.off()

## If we test a range of overdispersion and inbreeding, users can look at the chart for 
## a better way to estimate inbreeding if they know their overdispersion.

# Null alleles

# NAF is the frequency of a null allele
HHByParam_Null <- function(MAF = 0.05, NAF = 0.05, nsam = 500, depth = 20, overdispersion = 20,
                      inbreeding = 0, ploidy = 2, nloc = 5000){
  alleleFreq <- rep(c(1 - MAF - NAF, MAF, NAF), times = nloc)
  alleles2loc <- rep(1:nloc, each = 3)
  
  geno <- SimGenotypes(alleleFreq, alleles2loc, nsam, inbreeding, ploidy)
  alleleDepth <- SimAlleleDepth(matrix(depth, nrow = nsam, ncol = nloc,
                                       dimnames = list(NULL, as.character(1:nloc))),
                                geno, alleles2loc, overdispersion)
  rownames(alleleDepth) <- paste0("sam", 1:nsam)
  alleleDepth <- alleleDepth[,-seq(3, nloc * 3, by = 3)] # remove null alleles
  alleles2loc <- alleles2loc[-seq(3, nloc * 3, by = 3)]
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

tottest4 <- length(ploidies) * length(freqs2) * length(freqs)

testres4 <- data.frame(Ploidy = integer(tottest4),
                       MAF = numeric(tottest4),
                       NAF = numeric(tottest4),
                       Mean = numeric(tottest4),
                       Variance = numeric(tottest4),
                       Prop_estimated = numeric(tottest4))

currrow <- 1

for(p in ploidies){
  for(f1 in freqs2){
    for(f2 in freqs){
      res <- HHByParam_Null(MAF = f1, NAF = f2, nsam = 500L, depth = 20L,
                            overdispersion = 20, inbreeding = 0, ploidy = p, nloc = 5000L)
      testres4$Ploidy[currrow] <- p
      testres4$MAF[currrow] <- f1
      testres4$NAF[currrow] <- f2
      testres4$Mean[currrow] <- res$mean
      testres4$Variance[currrow] <- res$var
      testres4$Prop_estimated[currrow] <- res$nonNA
      currrow <- currrow + 1L
      if(currrow %% 10 == 0) print(currrow)
    }
  }
}

#save(testres4, file = "workspaces/variance_estimates_null_alleles.RData")
load("workspaces/variance_estimates_null_alleles.RData")

testres4$PloidyText <- ifelse(testres4$Ploidy == 2, "Diploid", "Tetraploid")

ggplot(testres4, aes(x = NAF, y = Mean)) +
  geom_line(aes(lty = as.factor(MAF))) +
  geom_point() +
  facet_grid(~ PloidyText) +
  labs(lty = "MAF", y = "Mean Hind/He", x = "Null allele frequency")

ggplot(testres4, aes(x = NAF, y = Variance)) +
  geom_line(aes(lty = as.factor(MAF))) +
  geom_point() +
  facet_grid(~ PloidyText)
