# Import read depth for tags from TagTaxaDist, only for tags from paralogs
# for experiment to ID paralogs

# tags to evaluate, and which loci they belong to
tagtab <- read.csv("marker_CSV/190515paralog_tags.csv", stringsAsFactors = FALSE)

# table indicating ploidy
accessions <- read.csv("~/DOE Msa study/Seq/GBSv2_180110/all_accession_names.csv",
                       stringsAsFactors = FALSE)
head(accessions)
diploids <- accessions$Accession[accessions$Ploidy == "2x"]    # 356
tetraploids <- accessions$Accession[accessions$Ploidy == "4x"] # 268

# find accessions in tagtaxadist file
ttdfile <- "D:/TASSELGBS_Msa/180813tagtaxadist.txt"
header <- strsplit(readLines(ttdfile, 1), split = "\t")[[1]]
header[1:50]
all(diploids %in% header)
all(tetraploids %in% header)

diploid_index <- match(diploids, header)
tetraploid_index <- match(tetraploids, header)

# set up matrices for output
diploid_mat <- matrix(0L, nrow = length(diploids), ncol = nrow(tagtab),
                      dimnames = list(diploids, tagtab$TagSeq))
tetraploid_mat <- matrix(0L, nrow = length(tetraploids), ncol = nrow(tagtab),
                      dimnames = list(tetraploids, tagtab$TagSeq))

# read through tagtaxadist file
mywhich <- list(character(), integer())
mywhich <- mywhich[c(1, rep(2, length(header) - 1))]

ttdcon <- file(ttdfile, open = "r")
depthdat <- scan(ttdcon, what = mywhich, nmax = 10000, sep = "\t", skip = 1)
tagsread <- 0L
while(length(depthdat[[1]])){
  # identify tags that we wanted to keep
  keep <- which(depthdat[[1]] %in% colnames(diploid_mat))
  if(length(keep) == 0) next
  thesetags <- depthdat[[1]][keep]
  
  # extract the counts for these tags for diploids and tetraploids
  dip_counts <- lapply(depthdat[diploid_index], function(x) x[keep])
  tet_counts <- lapply(depthdat[tetraploid_index], function(x) x[keep])
  
  # fill in matrices
  diploid_mat[,thesetags] <- matrix(unlist(dip_counts), nrow = length(diploids),
                                    ncol = length(thesetags), byrow = TRUE)
  tetraploid_mat[,thesetags] <- matrix(unlist(tet_counts), nrow = length(tetraploids),
                                       ncol = length(thesetags), byrow = TRUE)
  
  tagsread <- tagsread + length(thesetags)
  cat(paste(tagsread, "tags read"), sep = "\n")
  
  # read the next chunk
  depthdat <- scan(ttdcon, what = mywhich, nmax = 10000, sep = "\t")
}
close(ttdcon)

# check data
totdepthTags <- colSums(diploid_mat) + colSums(tetraploid_mat)
summary(totdepthTags)
mean(totdepthTags == 0) # 0.05%, probably just a few that were seqeuncing errors in 3x samples or blanks
hist(totdepthTags)
hist(log(totdepthTags))

totdepthTaxa <- c(rowSums(diploid_mat), rowSums(tetraploid_mat))
summary(totdepthTaxa)
hist(totdepthTaxa, breaks = 50)

#save(diploid_mat, tetraploid_mat, file = "workspaces/190515counts_matrices.RData")
