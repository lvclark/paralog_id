# Explore the dataset of loci with two alignment locations in the Miscanthus
# reference, to determine some good methods to use for the simulated annlealing
# algorithm to sort into isoloci.

twoalign <- read.csv("marker_CSV/190517twoalign.csv", stringsAsFactors = FALSE,
                     header = FALSE)
head(twoalign, 50)
head(twoalign[grepl("Chr01", twoalign$V1) &
                grepl("Chr02", twoalign$V2),], 50)

# pick a marker to look at
exampleMrkr <- twoalign[2158:2189,]
exampleMrkr

# get a statistic to use for transition probabilities from one paralog to other
oneToTwo <- 0.5 - (exampleMrkr$V5 - exampleMrkr$V4)/80
data.frame(exampleMrkr[,-3], oneToTwo)
hist(oneToTwo/sum(oneToTwo))
hist(oneToTwo)

# how many mismatches can we expect?
max(abs(twoalign$V4 - twoalign$V5)) # diff of 19 between them
max(c(twoalign$V4, twoalign$V5))    # max of 21/80
hist(abs(twoalign$V4 - twoalign$V5))

# get things for testing in Python
paste(exampleMrkr$V4, collapse = ", ")
paste(exampleMrkr$V5, collapse = ", ")
