# Get tag sequences for the paralogous markers in order to identify them in
# TagTaxaDist table and start import to polyRAD.
library(fastmatch)
library(Biostrings)

paralogs <- read.csv("marker_CSV/190513paralogs.csv", stringsAsFactors = FALSE)
head(paralogs)
tail(paralogs)
# eliminate a few from sorghum scaffolds
paralogs <- paralogs[startsWith(paralogs$Sbicolor_3, "Chr"),]
# eliminate miscanthus markers that match two sorghum positions
misctab <- table(paralogs$Query)
doubled <- names(misctab)[misctab > 1]
paralogs <- paralogs[!paralogs$Query %in% doubled,]

# function for getting total insertion or deletion size
IndelSize <- function(gregexprout, cigar){
  outlen <- 0L
  if(gregexprout[1] != -1){
    for(k in 1:length(gregexprout)){
      size <- substring(cigar, gregexprout[k],
                        gregexprout[k] + attr(gregexprout, "match.length")[k] - 2)
      outlen <- outlen + as.integer(size)
    }
  }
  return(outlen)
}

# list to hold tag sequences for each Miscanthus tag location
tagseq <- list()
length(tagseq) <- nrow(paralogs)
names(tagseq) <- paralogs$Query

tagsadded <- 0

# read through SAM file to find sequences
samfile <- "large_datasets/180110alignedtags.sam"
samcon <- file(samfile, open = "r")
while(length(samlines <- readLines(samcon, 10000))){
  # remove header lines
  samlines <- samlines[!grepl("^@", samlines)]
  if(length(samlines) == 0) next
  # remove unaligned sequence
  samlines <- samlines[!grepl("^tagSeq=[ACGT]+\t4\t", samlines)]
  if(length(samlines) == 0) next
  
  # split aligned lines and get marker info
  splitlines <- strsplit(samlines, split = "\t")
  botstrand <- sapply(splitlines, function(x) x[2]) == "16"
  botstrandN <- which(botstrand)
  positions <- as.integer(sapply(splitlines, function(x) x[4]))
  taglen <- nchar(sapply(splitlines, function(x) x[10]))
  cigar <- sapply(splitlines, function(x) x[6])
  # position correction for bottom strand, to match tagdigger
  insertpos <- gregexpr("[[:digit:]]+I", cigar[botstrand])
  deletepos <- gregexpr("[[:digit:]]+D", cigar[botstrand])
  for(i in 1:length(botstrandN)){
    j <- botstrandN[i] # j is index in all markers, i is index in bottom strand markers
    insertlen <- IndelSize(insertpos[[i]], cigar[j])
    deletelen <- IndelSize(deletepos[[i]], cigar[j])
    positions[j] <- positions[j] + taglen[j] - insertlen + deletelen - 1L
  }
  if(!is.integer(positions)) stop("positions somehow changed from integer.")
  
  # generate marker names
  markernames <- paste(sapply(splitlines, function(x) x[3]), # chromosome
                       formatC(positions,
                               width = 9, flag = "0"), # all are 9 digits wide
                       ifelse(botstrand, "bot", "top"),
                       sep = "-")
  # subset to markers from our paralogs list
  keepmrkr <- markernames %fin% paralogs$Query
  splitlines <- splitlines[keepmrkr]
  markernames <- markernames[keepmrkr]
  if(length(markernames) == 0) next
  
  # get tag sequences
  thesetags <- sapply(splitlines, function(x) x[10])
  # reverse complement for bottom strand
  botstrand2 <- grep("bot$", markernames)
  thesetags[botstrand2] <- 
    as.character(reverseComplement(DNAStringSet(thesetags[botstrand2])))
  
  # add sequences to set
  for(i in 1:length(markernames)){
    m <- markernames[i]
    tagseq[[m]] <- c(tagseq[[m]], thesetags[i])
  }
  
  # report results
  tagsadded <- tagsadded + length(markernames)
  cat(paste(tagsadded, "tags added to list"), sep = "\n")
}
close(samcon)


mylen <- sapply(tagseq, length)
summary(mylen)

# debug
mean(mylen == 0)
notags <- names(tagseq)[mylen == 0]

# tagseq[["Chr02-003541244-top"]]
# tagseq[["Chr02-095323814-top"]]
# 
# samcon <- file(samfile, open = "r")
# while(length(samlines <- readLines(samcon, 10000))){
# #  if(any(grepl("Chr01\t5877848", samlines))) break
# #  if(any(grepl("Chr02\t3541244", samlines))) break
#   if(any(grepl("Chr05\t32686449", samlines))) break
# }
# close(samcon)
# 
# grep("Chr02\t3541244", samlines, value = TRUE)
# grep("Chr05\t32686449", samlines, value = TRUE)
# grep("Chr05\t32686449", samlines)

# build a table of tag sequences, with marker names in miscanthus and sorghum
tagtab <- data.frame(TagSeq = unlist(tagseq),
                     Miscanthus = rep(names(tagseq), times = mylen),
                     stringsAsFactors = FALSE)
tagtab$Sorghum <- paralogs$Sbicolor_3[fmatch(tagtab$Miscanthus, paralogs$Query)]

head(tagtab)

#write.csv(tagtab, file = "marker_CSV/190515paralog_tags.csv", row.names = FALSE)
