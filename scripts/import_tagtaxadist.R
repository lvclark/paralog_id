# Import read depth for tags from TagTaxaDist, only for tags from paralogs
# for experiment to ID paralogs

# tags to evaluate, and which loci they belong to
tagtab <- read.csv("marker_CSV/190515paralog_tags.csv", stringsAsFactors = FALSE)

# table indicating ploidy
accessions <- read.csv("marker_CSV/all_accession_names.csv",
                       stringsAsFactors = FALSE)
head(accessions)
diploids <- accessions$Accession[accessions$Ploidy == "2x"]    # 356
tetraploids <- accessions$Accession[accessions$Ploidy == "4x"] # 268

writeLines(diploids, con = "marker_CSV/diploids.txt")
writeLines(tetraploids, con = "marker_CSV/tetraploids.txt")

readTagTaxaDist <- function(tags, diploids, tetraploids,
                            ttdfile = "D:/TASSELGBS_Msa/180813tagtaxadist.txt"){
  # find accessions in tagtaxadist file
  header <- strsplit(readLines(ttdfile, 1), split = "\t")[[1]]
  header[1:50]
  all(diploids %in% header)
  all(tetraploids %in% header)
  
  diploid_index <- match(diploids, header)
  tetraploid_index <- match(tetraploids, header)
  
  # set up matrices for output
  diploid_mat <- matrix(0L, nrow = length(diploids), ncol = length(tags),
                        dimnames = list(diploids, tags))
  tetraploid_mat <- matrix(0L, nrow = length(tetraploids), ncol = length(tags),
                           dimnames = list(tetraploids, tags))
  # read through tagtaxadist file
  mywhich <- list(character(), integer())
  mywhich <- mywhich[c(1, rep(2, length(header) - 1))]
  
  ttdcon <- file(ttdfile, open = "r")
  depthdat <- scan(ttdcon, what = mywhich, nmax = 10000, sep = "\t", skip = 1)
  tagsread <- 0L
  
  while(length(depthdat[[1]])){
    # identify tags that we wanted to keep
    keep <- which(depthdat[[1]] %in% tags)
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
    message(paste(tagsread, "tags read"))
    
    # read the next chunk
    depthdat <- scan(ttdcon, what = mywhich, nmax = 10000, sep = "\t")
  }
  close(ttdcon)
  
  return(list(diploid_mat, tetraploid_mat))
}

mats <- readTagTaxaDist(tagtab$TagSeq, diploids = diploids, tetraploids = tetraploids)
diploid_mat <- tagtab[[1]]
tetraploid_mat <- tagtab[[2]]
rm(mats)

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

# Instead import tags that aligned twice in the Miscanthus genome ####
twoalign <- read.csv("marker_CSV/190517twoalign.csv", stringsAsFactors = FALSE,
                     header = FALSE)
twoalign_Chr1_2 <- twoalign[grepl("Chr01", twoalign$V1) &
                              grepl("Chr02", twoalign$V2),]
mats <- readTagTaxaDist(twoalign_Chr1_2$V3, diploids = diploids, tetraploids = tetraploids)
mats[[1]][1:10,1:10]

#write.table(twoalign_Chr1_2, file = "marker_CSV/190525twoalign_Chr1Chr2.csv", col.names = FALSE,
#          row.names = FALSE, sep = ",")
#write.csv(t(mats[[1]]), file = "marker_CSV/190523diploid_Chr1Chr2.csv")
#write.csv(t(mats[[2]]), file = "marker_CSV/190523tetraploid_Chr1Chr2.csv")

table(table(paste(twoalign_Chr1_2$V1, twoalign_Chr1_2$V2)))
