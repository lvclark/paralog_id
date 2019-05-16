# Can we find negative correlations among read depth ratios of alleles that
# belong to the same isolocus?

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
# get correlations
for(i in 1:(length(theseal) - 1)){
  for(j in (i+1):length(theseal)){
    thiscor <- cor.test(depthRatio_dip_S[,theseal[i]],
                   depthRatio_dip_S[,theseal[j]])$p.value
    cormat[i,j] <- thiscor
    cormat[j,i] <- thiscor
    
    thisseqdist <- seqdist(tagseq[i], tagseq[j])
    seqmat[i,j] <- thisseqdist
    seqmat[j,i] <- thisseqdist
  }
}

image(cormat)
myclust <- hclust(as.dist(cormat))
plot(myclust)
grps <- cutree(myclust, k = 2)

# how does sequence similarity match up by group?
cat(tagseq[grps == 1], sep = "\n")
cat(tagseq[grps == 2], sep = "\n")

myclustseq <- hclust(as.dist(seqmat))
plot(myclustseq)
cutree(myclustseq, k = 2)
