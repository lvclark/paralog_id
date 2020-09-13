# look at log output from process_isoloci.py and see what temperatures worked well.
library(ggplot2)
library(viridis)

# function to process log file and extract information on markers that underwent
# simulated annealing.
extractLog <- function(file, finalTemp = FALSE, tabu = FALSE){
  if(finalTemp && tabu){
    stop("Cannot have both finalTemp and tabu")
  }
  logout <- readLines(file)
  initHindHeLines <- grep("^Initial Hind/He:", logout)
  finalHindHeLines <- initHindHeLines + 2
  # filter out where isoloci were fixed
  keep <- grepl("^Final Hind/He:", logout[finalHindHeLines])
  initHindHeLines <- initHindHeLines[keep]
  finalHindHeLines <- finalHindHeLines[keep]
  
  markerLines <- initHindHeLines - 1
  initNMLines <- initHindHeLines + 1
  finalNMLines <- initHindHeLines + 3
  if(finalTemp){
    finalTempLines <- initHindHeLines + 4
  }
  if(tabu){
    repLines <- initHindHeLines + 4
  }
  
  # extract numeric values
  init_NM <- as.numeric(sub("Initial average NM: ", "", logout[initNMLines]))
  final_NM <- as.numeric(sub("Final average NM: ", "", logout[finalNMLines]))
  
  initHindHeSplit <- strsplit(logout[initHindHeLines], " ")
  finalHindHeSplit <- strsplit(logout[finalHindHeLines], " ")
  init1_HindHe <- as.numeric(sapply(initHindHeSplit,
                                    function(x){
                                      out <- x[3]
                                      out[out == "None"] <- NA
                                      return(out)}))
  init2_HindHe <- as.numeric(sapply(initHindHeSplit,
                                    function(x){
                                      out <- x[4]
                                      out[out == "None"] <- NA
                                      return(out)}))
  final1_HindHe <- as.numeric(sapply(finalHindHeSplit,
                                    function(x){
                                      out <- x[3]
                                      out[out == "None"] <- NA
                                      return(out)}))
  final2_HindHe <- as.numeric(sapply(finalHindHeSplit,
                                    function(x){
                                      out <- x[4]
                                      out[out == "None"] <- NA
                                      return(out)}))
  
  # get maximum initial and final Hind/He
  initMax_HindHe <- apply(cbind(init1_HindHe, init2_HindHe), 1, max, na.rm = TRUE)
  finalMax_HindHe <- apply(cbind(final1_HindHe, final2_HindHe), 1, max, na.rm = TRUE)
  
  # construct data frame
  outdf <- data.frame(row.names = logout[markerLines],
                      Init_NM = init_NM,
                      Final_NM = final_NM,
                      Init1_HindHe = init1_HindHe,
                      Init2_HindHe = init2_HindHe,
                      InitMax_HindHe = initMax_HindHe,
                      Final1_HindHe = final1_HindHe,
                      Final2_HindHe = final2_HindHe,
                      FinalMax_HindHe = finalMax_HindHe)
  if(finalTemp){
    outdf$Final_Temp <- as.numeric(sub("Final temperature: ", "", logout[finalTempLines]))
  }
  if(tabu){
    outdf$Reps <- as.integer(sub("Rep where best solution found: ", "", logout[repLines]))
  }
  
  return(outdf)
}

# Tabu search Dec. 2019 ####
dfTabu2x <- extractLog("log/Msa_1_diploids.log", tabu = TRUE)
dfTabu4x <- extractLog("log/Msa_1_tetraploids.log", tabu = TRUE)

dfTabu2x <- dfTabu2x[is.finite(dfTabu2x$FinalMax_HindHe),]
dfTabu4x <- dfTabu4x[is.finite(dfTabu4x$FinalMax_HindHe),]

ggplot(dfTabu2x, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Reps)) +
  geom_point() +
  geom_hline(yintercept = 1/2 * 0.85, lty = 2) +
  geom_vline(xintercept = 1/2 * 0.85, lty = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  labs(x = expression("Initial" ~ H[ind]/H[E] ~ "(max out of two loci)"),
       y = expression("Final" ~ H[ind]/H[E] ~ "(max out of two loci)"),
       col = "Tabu reps")

ggplot(dfTabu4x, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Reps)) +
  geom_point() +
  geom_hline(yintercept = 3/4 * 0.85, lty = 2) +
  geom_vline(xintercept = 3/4 * 0.85, lty = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  labs(x = expression("Initial" ~ H[ind]/H[E] ~ "(max out of two loci)"),
       y = expression("Final" ~ H[ind]/H[E] ~ "(max out of two loci)"),
       col = "Tabu reps")

dfTabuAll <- rbind(dfTabu2x, dfTabu4x)
dfTabuAll$Ploidy <- rep(c("Diploids", "Tetraploids"), times = c(nrow(dfTabu2x), nrow(dfTabu4x)))

expvals <- data.frame(Ploidy = c("Diploids", "Tetraploids"),
                      Val = c(1/2 * 0.85, 3/4 * 0.85))

tiff("191219tabu_output.tiff", width = 6.5 * 300, height = 6.5 * 300, res = 300,
     compression = "lzw")
ggplot(dfTabuAll, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Reps)) +
  geom_point() +
  geom_hline(mapping = aes(yintercept = Val), data = expvals, lty = 2) +
  geom_vline(mapping = aes(xintercept = Val), data = expvals, lty = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  labs(x = expression("Initial" ~ H[ind]/H[E] ~ "(max out of two loci)"),
       y = expression("Final" ~ H[ind]/H[E] ~ "(max out of two loci)"),
       col = "Tabu reps") +
  facet_wrap(~ Ploidy, scales = "free_x", nrow = 2, ncol = 1)
dev.off()

success <- dfTabuAll$FinalMax_HindHe <= ifelse(dfTabuAll$Ploidy == "Diploids", 1/2 * 0.85, 3/4 * 0.85)
mean(success,
     na.rm = TRUE)
mean(dfTabuAll$Reps[success] == 0)
hist(dfTabuAll$Reps[success])
table(dfTabuAll$Reps[success])
