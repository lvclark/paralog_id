# Estimate variance and bias in Hind/He
library(polyRAD) # v1.6, January 2022
library(ggplot2)
library(dplyr)
library(gridExtra)

# function to get an estimate of hind/he and its variance, given biological and technical parameters.
# MAF = minor allele frequency.  Biallelic loci only.
# nsam = number of samples
# depth = read depth (one value; will be uniform across samples and loci)
# overdispersion = overdispersion parameter
# inbreeding = inbreeding coefficient
HHByParam <- function(MAF = 0.05, nsam = 500, depth = 20, overdispersion = 20,
                      inbreeding = 0, ploidy = 2, nloc = 5000,
                      contamRate = 0, errorRate = 0){
  alleleFreq <- rep(c(1 - MAF, MAF), times = nloc)
  alleles2loc <- rep(1:nloc, each = 2)
  
  geno <- SimGenotypes(alleleFreq, alleles2loc, nsam, inbreeding, ploidy)
  alleleDepth <- SimAlleleDepth(matrix(depth, nrow = nsam, ncol = nloc,
                                       dimnames = list(NULL, as.character(1:nloc))),
                                geno, alleles2loc, overdispersion, contamRate,
                                errorRate)
  rownames(alleleDepth) <- paste0("sam", 1:nsam)
  simrad <- RADdata(alleleDepth, alleles2loc,
                    locTable = data.frame(row.names = paste0("loc", 1:nloc)),
                    possiblePloidies = list(as.integer(ploidy)),
                    contamRate = contamRate,
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
samsizes <- c(50, 100, 500, 1000)
depths <- c(2, 5, 10, 20, 50, 100, 200)
contams <- c(0, 0.001, 0.01)
errors <- c(0, 0.001)

tottests <- length(ploidies) * length(freqs) * length(samsizes) * length(depths) *
  length(contams) * length(errors)

testres <- data.frame(Ploidy = integer(tottests),
                      MAF = numeric(tottests),
                      N_sam = numeric(tottests),
                      ContamRate = numeric(tottests),
                      ErrorRate = numeric(tottests),
                      Depth = numeric(tottests),
                      Mean = numeric(tottests),
                      Variance = numeric(tottests),
                      Prop_estimated = numeric(tottests))
currrow <- 1

for(p in ploidies){
  for(f in freqs){
    for(s in samsizes){
      for(d in depths){
        for(c in contams){
          for(e in errors){
            res <- HHByParam(MAF = f, nsam = s, depth = d, ploidy = p,
                             overdispersion = 20, inbreeding = 0, nloc = 5000,
                             contamRate = c, errorRate = e)
            testres$Ploidy[currrow] <- p
            testres$MAF[currrow] <- f
            testres$N_sam[currrow] <- s
            testres$Depth[currrow] <- d
            testres$ContamRate[currrow] <- c
            testres$ErrorRate[currrow] <- e
            testres$Mean[currrow] <- res$mean
            testres$Variance[currrow] <- res$var
            testres$Prop_estimated[currrow] <- res$nonNA
            currrow <- currrow + 1
            if(currrow %% 10 == 0) print(currrow)
          }
        }
      }
    }
  }
}

#save(testres, file = "workspaces/variance_estimates_2022-01-16.RData")
load("workspaces/variance_estimates_2022-01-16.RData")

testres$PloidyText <- ifelse(testres$Ploidy == 2, "Diploid", "Tetraploid")

# plot results
p1 <- ggplot(testres[testres$Depth < 200 & testres$ContamRate == 0 & testres$ErrorRate == 0,],
       aes(x = Depth, y = sqrt(Variance), color = as.factor(MAF), group = MAF)) +
  geom_line() +
  facet_grid(PloidyText ~ N_sam, labeller = labeller(N_sam = function(x) paste("N =", x))) +
  scale_x_continuous(breaks = depths, trans = "log2") +
  labs(y = "Standard deviation of estimate", color = "MAF", x = "Read depth")
p1
# ggsave("~/NSF polyRAD/Year 4 report/fig1a.png",
#        width = 6.5, height = 3.5)

# bias -- overestimated for rare alleles and small sample size, since so many NA
p2 <- ggplot(testres[testres$Depth < 200 & testres$ContamRate == 0 & testres$ErrorRate == 0,],
       aes(x = Depth, y = Mean, color = as.factor(MAF), group = MAF)) +
  geom_line() +
  facet_grid(PloidyText ~ N_sam, labeller = labeller(N_sam = function(x) paste("N =", x))) +
  scale_x_continuous(breaks = depths, trans = "log2") +
  labs(y = "Mean estimate", color = "MAF", x = "Read depth")
p2
# ggsave("~/NSF polyRAD/Year 4 report/fig1b.png",
#        width = 6.5, height = 3.5)

#cairo_pdf("Fig4_samplesize_depth_maf.pdf", width = 6.7, height = 6)
grid.arrange(arrangeGrob(p1 + ggtitle("A"), p2 + ggtitle("B")))
#dev.off()


# Explore overdispersion, inbreeding, and sequencing error.  polyRAD 1.6 Jan 2022.
ods <- 5:20
freqs2 <- c(0.01, 0.05)
inbreed <- seq(0, 1, by = 0.1)

tottests2 <- length(ploidies) * length(ods) * length(freqs2) * length(inbreed) * length(errors)
testres2 <- data.frame(Ploidy = integer(tottests2),
                      MAF = numeric(tottests2),
                      Overdispersion = numeric(tottests2),
                      Inbreeding = numeric(tottests2),
                      ErrorRate = numeric(tottests2),
                      Mean = numeric(tottests2),
                      Variance = numeric(tottests2),
                      Prop_estimated = numeric(tottests2))
currrow <- 1

for(p in ploidies){
  for(f in freqs2){
    for(o in ods){
      for(i in inbreed){
        for(e in errors){
          res <- HHByParam(MAF = f, nsam = 500L, depth = 20L, ploidy = p,
                           overdispersion = o, inbreeding = i, nloc = 20000,
                           errorRate = e)
          testres2$Ploidy[currrow] <- p
          testres2$MAF[currrow] <- f
          testres2$Overdispersion[currrow] <- o
          testres2$Inbreeding[currrow] <- i
          testres2$ErrorRate[currrow] <- e
          testres2$Mean[currrow] <- res$mean
          testres2$Variance[currrow] <- res$var
          testres2$Prop_estimated[currrow] <- res$nonNA
          currrow <- currrow + 1L
          if(currrow %% 10 == 0) print(currrow)
        }
      }
    }
  }
}

testres2$PloidyText <- ifelse(testres2$Ploidy == 2, "Diploid", "Tetraploid")

#save(testres2, file = "workspaces/variance_estimates_overdispersion_inbreeding_2022-01-17.RData")
load("workspaces/variance_estimates_overdispersion_inbreeding_2022-01-17.RData")

# figure for manuscript
#cairo_pdf("SuppFig2_overdispersion_inbreeding.pdf", width = 6.7, height = 3.5)
# tiff("SuppFig2_overdispersion_inbreeding.tiff", width = 6.7 * 300, height = 3.5 * 300,
#      res = 300, compression = "lzw")
testres2 %>%
  filter(ErrorRate == 0) %>%
ggplot(mapping = aes(x = Overdispersion, y = Mean, color = Inbreeding, group = paste(Inbreeding, MAF))) +
  geom_line(aes(lty = as.factor(MAF))) +
  facet_grid(~ PloidyText, scales = "free_y") +
  scale_color_viridis_c() +
  labs(lty = "MAF", y = "Mean estimate")
# dev.off()

testres2 %>%
  filter(ErrorRate == 0.001) %>%
  ggplot(mapping = aes(x = Overdispersion, y = Mean, color = Inbreeding, group = paste(Inbreeding, MAF))) +
  geom_line(aes(lty = as.factor(MAF))) +
  facet_grid(~ PloidyText, scales = "free_y") +
  scale_color_viridis_c() +
  labs(lty = "MAF", y = "Mean estimate")

testres2 %>%
  filter(Overdispersion == 20) %>%
  ggplot(aes(x = Inbreeding, y = Mean, linetype = as.factor(MAF),
             color = as.factor(ErrorRate), shape = as.factor(MAF))) +
  geom_line() +
  geom_point() +
  facet_grid(~ PloidyText, scales = "free_y") +
  labs(lty = "MAF", shape = "MAF",
       y = "Mean estimate", color = "Sequencing error rate") +
  scale_color_manual(values = c("darkgreen", "dodgerblue"))

testres2 %>%
  filter(Overdispersion == 20) %>%
  mutate(MAFtext = paste("MAF = ", MAF)) %>%
  ggplot(aes(x = Inbreeding, y = Mean, linetype = as.factor(ErrorRate),
             color = PloidyText, shape = as.factor(ErrorRate))) +
  geom_line() +
  geom_point() +
  facet_grid(~ MAFtext) +
  labs(lty = "Sequencing error rate", shape = "Sequencing error rate",
       y = "Mean estimate", color = "Ploidy") +
  scale_color_manual(values = c("darkgreen", "dodgerblue"))

# Do this fig for hexaploids too since there is interest from sweetpotato
p <- 6L

tottests3 <- length(ods) * length(freqs2) * length(inbreed) * length(errors)
testres3 <- data.frame(Ploidy = 6L,
                       MAF = numeric(tottests3),
                       Overdispersion = numeric(tottests3),
                       Inbreeding = numeric(tottests3),
                       ErrorRate = numeric(tottests3),
                       Mean = numeric(tottests3),
                       Variance = numeric(tottests3),
                       Prop_estimated = numeric(tottests3),
                       PloidyText = "Hexaploid")
currrow <- 1

for(f in freqs2){
  for(o in ods){
    for(i in inbreed){
      for(e in errors){
        res <- HHByParam(MAF = f, nsam = 500L, depth = 20L, ploidy = p,
                         overdispersion = o, inbreeding = i, nloc = 20000,
                         errorRate = e)
        testres3$MAF[currrow] <- f
        testres3$Overdispersion[currrow] <- o
        testres3$Inbreeding[currrow] <- i
        testres2$ErrorRate[currrow] <- e
        testres3$Mean[currrow] <- res$mean
        testres3$Variance[currrow] <- res$var
        testres3$Prop_estimated[currrow] <- res$nonNA
        currrow <- currrow + 1L
        if(currrow %% 10 == 0) print(currrow)
      }
    }
  }
}

save(testres3, file = "workspaces/variance_estimates_overdispersion_inbreeding_hexaploid_2022-01-22.RData")

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

#cairo_pdf("SuppFig3_null_alleles.pdf", width = 3.35, height = 5.2)
# tiff("SuppFig3_null_alleles.tiff", width = 3.35 * 300, height = 5.2 * 300,
#      res = 300, compression = "lzw")
ggplot(testres4, aes(x = NAF, y = Mean)) +
  geom_line(aes(lty = as.factor(MAF))) +
  geom_point() +
  facet_wrap(~ PloidyText, nrow = 2, scales = "free_y") +
  labs(lty = "MAF", y = "Mean estimate", x = "Null allele frequency") +
  theme(legend.position = "bottom")
# dev.off()

ggplot(testres4, aes(x = NAF, y = Variance)) +
  geom_line(aes(lty = as.factor(MAF))) +
  geom_point() +
  facet_grid(~ PloidyText)

# Effects of contamination and sequencing error ####
testres %>%
  filter(N_sam == 500 & MAF == 0.05) %>%
  ggplot(aes(x = Depth, y = Variance, color = as.character(ErrorRate))) +
  geom_line() +
  facet_grid(PloidyText ~ ContamRate)

testres %>%
  filter(N_sam == 500 & MAF == 0.05) %>%
  ggplot(aes(x = Depth, y = Mean, color = as.character(ErrorRate))) +
  geom_line() +
  facet_grid(PloidyText ~ ContamRate)

# As expected, contamination and seq error bias slightly upwards
# Variance doesn't seem to be affected

testres$MAFtext <- paste("MAF =", testres$MAF)

testres %>%
  filter(N_sam == 500, Depth == 20) %>%
  ggplot(aes(x = ContamRate, y = Mean)) +
  geom_line(aes(linetype = as.character(ErrorRate))) +
  geom_point(pch = 1) +
  facet_grid(PloidyText ~ MAFtext, scales = "free_y") +
  labs(x = "Contamination rate", y = "Mean estimate",
       linetype = "Sequencing error rate") +
  scale_x_continuous(breaks = seq(0, 0.01, by = 0.001),
                     minor_breaks = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Slight increase with high contamination rate
# Sequencing error causes more inflation at lower MAF

testres %>%
  filter(N_sam == 500, MAF == 0.05) %>%
  ggplot(aes(x = ContamRate, y = Mean)) +
  geom_line(aes(linetype = as.character(ErrorRate))) +
  geom_point(pch = 1) +
  facet_grid(PloidyText ~ Depth, scales = "free_y") +
  labs(x = "Contamination rate", y = "Mean estimate",
       linetype = "Sequencing error rate") +
  scale_x_continuous(breaks = seq(0, 0.01, by = 0.001),
                     minor_breaks = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Estimates are higher at depth below 10, but we already knew that

testres %>%
  filter(MAF == 0.05, Depth == 20) %>%
  ggplot(aes(x = ContamRate, y = Mean)) +
  geom_line(aes(linetype = as.character(ErrorRate))) +
  geom_point(pch = 1) +
  facet_grid(PloidyText ~ N_sam, scales = "free_y") +
  labs(x = "Contamination rate", y = "Mean estimate",
       linetype = "Sequencing error rate") +
  scale_x_continuous(breaks = seq(0, 0.01, by = 0.001),
                     minor_breaks = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Estimates are higher at sample size below 500, but we already knew that

# Test error effects across inbreeding and MAF values

