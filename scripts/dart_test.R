# Test Hind/He on Jeff Endelman's DArT markers
library(polyRAD)
library(ggplot2)
library(viridis)

load("workspaces/updog_results.rda")

ans[[1]]

nLoc <- length(ans) # 678 loci
nInd <- length(ans[[1]]$postmean) # 89 individuals

wmgeno <- matrix(NA, nrow = nInd, ncol = nLoc)
alleleDepths <- matrix(NA_integer_, nrow = nInd, ncol = nLoc * 2)

for(i in 1:nLoc){
  wmgeno[,i] <- ans[[i]]$postmean
  alleleDepths[,(2 * i - 1)] <- ans[[i]]$input$refvec
  alleleDepths[, 2 * i] <- ans[[i]]$input$sizevec - ans[[i]]$input$refvec
}

# look at posterior mean genotypes to try to id parents
mydist <- dist(wmgeno)
mypcoa <- cmdscale(mydist)
plot(mypcoa[,1], mypcoa[,2]) # three groups, so I'm not really sure what is going on here.

rownames(alleleDepths) <- paste0("Sam", 1:nInd)
colnames(alleleDepths) <- paste(rep(paste0("Loc", 1:nLoc), each = 2), c("A", "G"), sep = "_")

# set up RADdata object
myRAD <- RADdata(alleleDepths, rep(1:nLoc, each = 2), 
                 data.frame(row.names = paste0("Loc", 1:nLoc)), list(4L),
                 0.001, rep(c("A", "G"), times = nLoc))
myRAD

# Hind/He
myhh <- HindHe(myRAD)

hh_byInd <- rowMeans(myhh, na.rm = TRUE)
hist(hh_byInd)

ggplot(mapping = aes(x = mypcoa[,1], y = mypcoa[,2], col = hh_byInd)) +
  geom_point() +
  scale_color_viridis()
# one group more heterozygous than others

hh_byLoc <- colMeans(myhh, na.rm = TRUE)
hist(hh_byLoc, breaks = 30, xlab = "Hind/He", main = "Histogram of Hind/He across loci")
abline(v = 0.75, col = "blue", lwd = 2)

# Classify markers based on Hind/He
markers_out <- data.frame(marker_number = 1:nLoc,
                          HindHe = hh_byLoc,
                          Class = character(nLoc),
                          stringsAsFactors = FALSE)
markers_out$Class[hh_byLoc < 0.6] <- "low"
markers_out$Class[hh_byLoc > 0.9] <- "high"
markers_out$Class[hh_byLoc >= 0.6 & hh_byLoc <= 0.9] <- "ok"

table(markers_out$Class)

# plot markers in different categories
plotmarker <- function(i, rad = myRAD){
  plot(rad$alleleDepth[,(2 * i - 1)], rad$alleleDepth[,2 * i],
       xlab = "Reference counts", ylab = "Alternative counts",
       main = i)
  abline(a = 0, b = 1/3)
  abline(a = 0, b = 1)
  abline(a = 0, b = 3)
  abline(h = 0)
  abline(v = 0)
}

high_locs <- sample(which(markers_out$Class == "high"), 9)
low_locs <- sample(which(markers_out$Class == "low"), 9)
ok_locs <- sample(which(markers_out$Class == "ok"), 9)

par(mfrow = c(3, 3))
for(L in high_locs){
  plotmarker(L)
}
for(L in low_locs){
  plotmarker(L)
}
for(L in ok_locs){
  plotmarker(L)
}

really_low <- sample(which(markers_out$HindHe < 0.3), 9)
par(mfrow = c(3, 3))
for(L in really_low){
  plotmarker(L)
}

write.csv(markers_out, file = "hindhe_by_marker.csv", row.names = FALSE)
