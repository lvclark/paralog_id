# Can we find negative correlations among read depth ratios of alleles that
# belong to the same isolocus?
Rcpp::sourceCpp("src/simpson.cpp")

load("workspaces/190515counts_matrices.RData")
tagtab <- read.csv("marker_CSV/190515paralog_tags.csv", stringsAsFactors = FALSE)

# Set up some objects as in H_stats_sorghum_miscanthus.R
alleles2locM <- match(tagtab$Miscanthus, unique(tagtab$Miscanthus))
alleles2locS <- match(tagtab$Sorghum, unique(tagtab$Sorghum))
locDepth_dip_Sorg <- t(rowsum(t(diploid_mat), alleles2locS))

str(locDepth_dip_Sorg)
str(diploid_mat)
locDepth_dip_Sorg[1:10,1:10]

# get depth ratios, based on sorghum reference
depthRatio_dip_S <- diploid_mat / locDepth_dip_Sorg[,as.character(alleles2locS)]
str(depthRatio_dip_S)

depthRatio_dip_S[1:10,alleles2locS == 6]

# function for distance between sequences ####
seqdist <- function(seq1, seq2){
  return(sum(strsplit(seq1, "")[[1]] != strsplit(seq2, "")[[1]]))
}

# exploration - get info on a locus of choice ####
sorgloc <- 6 # locus to explore from sorghum
theseal <- which(alleles2locS == sorgloc) # columns containing these alleles
miscloc <- unique(alleles2locM[theseal]) # corresponding locus numbers in Miscanthus
if(length(miscloc) != 2) stop("Check Miscanthus loci.")
miscal1 <- which(alleles2locM == miscloc[1]) # columns containing these alleles in Miscanthus
miscal2 <- which(alleles2locM == miscloc[2])
tagseq <- colnames(diploid_mat)[theseal]

alnames <- c(paste("MiscAl1", miscal1, sep = "."),
             paste("MiscAl2", miscal2, sep = "."))
cormat <- matrix(NA_real_, nrow = length(theseal), ncol = length(theseal),
                 dimnames = list(alnames, alnames))
seqmat <- cormat
cormatBin <- cormat
cormatC <- cormat
# get correlations
for(i in 1:(length(theseal) - 1)){
  for(j in (i+1):length(theseal)){
    thiscor <- cor.test(depthRatio_dip_S[,theseal[i]],
                   depthRatio_dip_S[,theseal[j]],
                   alternative = "less")
    cormat[i,j] <- thiscor$p.value
    cormat[j,i] <- thiscor$p.value
    cormatC[i,j] <- thiscor$estimate
    cormatC[j,i] <- thiscor$estimate
    
    thisseqdist <- seqdist(tagseq[i], tagseq[j])
    seqmat[i,j] <- thisseqdist
    seqmat[j,i] <- thisseqdist
    
    # thiscorBin <- fisher.test(diploid_mat[,theseal[i]] > 0,
    #                           diploid_mat[,theseal[j]] > 0,
    #                           alternative = "less")
  }
}

image(cormat)
heatmap(cormat)
myclust <- hclust(as.dist(cormat))
plot(myclust)
grps <- cutree(myclust, k = 2)

# how does sequence similarity match up by group?
cat(tagseq[grps == 1], sep = "\n")
cat(tagseq[grps == 2], sep = "\n")

myclustseq <- hclust(as.dist(seqmat))
plot(myclustseq)
cutree(myclustseq, k = 2)

image(cormatC)
heatmap(cormatC)

# OR, rather than starting with clustering, use alignment as initial
# configuration for sorting algorithm.

# how did this perform in terms of Hind/He
h_orig <- HReadDepth(diploid_mat[,theseal], alleles2locM[theseal] - miscloc[1] + 1, 2)
h_cor <- HReadDepth(diploid_mat[,theseal], grps, 2)

colMeans(h_orig[[1]], na.rm=TRUE)/h_orig[[2]]
colMeans(h_cor[[1]], na.rm=TRUE)/h_cor[[2]]

# Use depth ratios, but don't include the locus for comparison in the denominator ####
# countsmat is a depth matrix for one locus
# al1 and al2 are numeric indices for the two alleles of interest
exclusiveDepthRatios <- function(countsmat, al1, al2){
  tot1 <- rowSums(countsmat[,-al1, drop = FALSE])
  tot2 <- rowSums(countsmat[,-al2, drop = FALSE])
  rat1 <- countsmat[,al1]/tot2
  rat2 <- countsmat[,al2]/tot1
  
  # fill in zeros where the other allele took up all reads
  rat1[is.nan(rat1) & !is.nan(rat2)] <- 0
  rat2[is.nan(rat2) & !is.nan(rat1)] <- 0
  
  return(list(rat1, rat2))
}

testrat <- exclusiveDepthRatios(diploid_mat[,theseal], 1, 2)
image(cbind(testrat[[1]], testrat[[2]]))
summary(testrat[[1]])
summary(testrat[[2]])

head(diploid_mat[,theseal])
cor.test(testrat[[1]], testrat[[2]])

cormat <- matrix(NA_real_, nrow = length(theseal), ncol = length(theseal),
                 dimnames = list(alnames, alnames))
cormatP <- cormat
for(i in 1:(length(theseal) - 1)){
  for(j in (i+1):length(theseal)){
    thisrat <- exclusiveDepthRatios(diploid_mat[,theseal], i, j)
    testout <- cor.test(thisrat[[1]], thisrat[[2]], alternative = "less", method = "kendall")
    cormat[i,j] <- cormat[j,i] <- testout$estimate
    cormatP[i,j] <- cormatP[j,i] <- testout$p.value
  }
}

image(cormat)
cormat
image(cormatP)
heatmap(cormatP)
sum(cormatP[lower.tri(cormatP)] < 0.05)
which(cormatP < 0.05, arr.ind = TRUE)
# for locus 6
alnames[c(1, 7, 8, 11)]
alnames[c(6, 13)]
image(diploid_mat[,theseal[c(1, 7, 8, 11)]])
image(diploid_mat[,theseal[c(1, 11, 8, 7)]])
cormatP[c(1, 7, 8, 11), c(1, 7, 8, 11)]
image(diploid_mat[,theseal[c(6, 13)]])
hist(diploid_mat[,theseal[6]])
hist(diploid_mat[,theseal[13]])

testrat <- exclusiveDepthRatios(diploid_mat[,theseal], 6, 13)
hist(testrat[[1]])
hist(testrat[[2]])
cormat[6,13]
cormatP[6,13]
cor.test(testrat[[1]], testrat[[2]], alternative = "less")

testrat <- exclusiveDepthRatios(diploid_mat[,theseal], 7, 8)
plot(testrat[[1]], testrat[[2]])
cor.test(testrat[[1]], testrat[[2]], alternative = "less")
cor.test(testrat[[1]], testrat[[2]], alternative = "two.sided")
cor.test(testrat[[1]], testrat[[2]], alternative = "less", method = "kendall") # actually much more significant
cor.test(testrat[[1]], testrat[[2]], alternative = "less", method = "spearman") # warns about ties

seqmat[c(1, 7, 8, 11),c(1, 7, 8, 11)]

# for locus 33
algrp <- c(1, 2, 6, 9, 10, 12, 13, 14, 15, 20)
cormatP[algrp, algrp]
image(diploid_mat[,theseal[algrp]])
image(diploid_mat[,theseal[c(1,9,14)]])
image(diploid_mat[,theseal[-algrp]]) # no common alleles left over

h33 <- HReadDepth(diploid_mat[,theseal], rep(1, length(theseal)), 1)
mean(h33[[1]], na.rm = TRUE)/h33[[2]] # 0.80; too high
h33a <- HReadDepth(diploid_mat[,theseal[algrp]], rep(1, length(algrp)), 1)
mean(h33a[[1]], na.rm = TRUE)/h33a[[2]] # still 0.80

# Perhaps keep lowering the p-value threshold, and find one where the
# groups don't violate Hind/He expectations.
pthresh = 1 # -log10 p-value at which to group alleles
repeat{
  whicharr <- which(cormatP < 10 ^ -pthresh, arr.ind = TRUE)
  if(nrow(whicharr) == 0){
    # no groups can be made
    grps <- list()
    break
  }
  grps <- list(unname(whicharr[1,]))
  for(i in 2:nrow(whicharr)){
    if(i > nrow(whicharr)) break # for cases of only one row
    matches1 <- which(sapply(grps, function(x) whicharr[i,1] %in% x))
    matches2 <- which(sapply(grps, function(x) whicharr[i,2] %in% x))
    if(length(matches1) > 1 || length(matches2) > 1){
      stop("Too many matches.")
    }
    if(length(matches1) == 0 && length(matches2) == 0){
      # new group
      grps[[length(grps) + 1]] <- unname(whicharr[i,])
    } else if(length(matches1) == 0 && length(matches2) == 1){
      # add allele to existing group
      grps[[matches2]] <- c(grps[[matches2]], whicharr[i,1])
    } else if(length(matches1) == 1 && length(matches2) == 0){
      grps[[matches1]] <- c(grps[[matches1]], whicharr[i,2])
    } else if(length(matches1) == 1 && length(matches2) == 1){
      if(matches1 != matches2){
        # merge groups
        grps[[matches1]] <- c(grps[[matches1]], grps[[matches2]])
        grps <- grps[-matches2]
      }
    }
  }
  grps <- lapply(grps, sort)
  hindhe <- sapply(grps,
                   function(x){
                     h <- HReadDepth(diploid_mat[,theseal[x]],
                                     rep(1, length(x)), 1)
                     mean(h[[1]], na.rm = TRUE)/h[[2]]
                   })
  if(all(hindhe < 0.6)) break
  pthresh <- pthresh + 1
}

grps
hindhe
pthresh
alnames[grps[[1]]]
alnames[grps[[2]]]
cat(tagseq[grps[[1]]], sep = "\n")
cat(tagseq[grps[[2]]], sep = "\n")
image(diploid_mat[,theseal[grps[[1]]]])
image(diploid_mat[,theseal[grps[[2]]]])
cormatP[grps[[1]], grps[[1]]]
