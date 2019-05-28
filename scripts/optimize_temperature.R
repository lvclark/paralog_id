# look at log output from process_isoloci.py and see what temperatures worked well.
library(ggplot2)
library(viridis)

# function to process log file and extract information on markers that underwent
# simulated annealing.
extractLog <- function(file){
  logout <- readLines(file)
  initHindHeLines <- grep("^Initial Hind/He:", logout)
  finalHindHeLines <- initHindHeLines + 2
  # filter out where isoloci were fixed
  keep <- grepl("^Final Hind/He:", logout[finalHindHeLines])
  initHindHeLines <- initHindHeLines[keep]
  finalHindHeLines <- finalHindHeLines[keep]
  
  markerLines <- initHindHeLines - 2
  initNMLines <- initHindHeLines + 1
  finalNMLines <- initHindHeLines + 3
  unsplitHindHeLines <- initHindHeLines - 1
  
  # extract numeric values
  unsplit_HindHe <- as.numeric(sub("Hind/He without splitting: ", "", logout[unsplitHindHeLines]))
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
                      Unsplit_HindHe = unsplit_HindHe,
                      Init_NM = init_NM,
                      Final_NM = final_NM,
                      Init1_HindHe = init1_HindHe,
                      Init2_HindHe = init2_HindHe,
                      InitMax_HindHe = initMax_HindHe,
                      Final1_HindHe = final1_HindHe,
                      Final2_HindHe = final2_HindHe,
                      FinalMax_HindHe = finalMax_HindHe)
  
  return(outdf)
}

# extract values and plot ####
df50 <- extractLog("log/190528hindhelogT50.txt")
head(df50)

ggplot(df50, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Final_NM - Init_NM)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))

ggplot(df50) +
  geom_density(aes(x = FinalMax_HindHe), col = "blue") +
  geom_density(aes(x = InitMax_HindHe), col = "red")

ggplot(df50[df50$InitMax_HindHe > 0.5,]) +
  geom_density(aes(x = FinalMax_HindHe), col = "blue") +
  geom_density(aes(x = InitMax_HindHe), col = "red") +
  geom_vline(xintercept = 0.5)

# Notes from T0 = 0.50:
# Doesn't fix everything but many are considerably improved, mean is moved to a reasonable place.
# If initial Hind/He for isoloci is not higher than expected, do not bother shuffling them around.

# try with T0 = 0.25
df25 <- extractLog("log/190528hindhelogT25.txt")

ggplot(df25, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Final_NM - Init_NM)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))

ggplot(df25[df25$InitMax_HindHe > 0.5,]) +
  geom_density(aes(x = FinalMax_HindHe), col = "blue") +
  geom_density(aes(x = InitMax_HindHe), col = "red") +
  geom_vline(xintercept = 0.5) 

# T0 = 0.25 is the same as or slightly better than T0 = 0.5

# try with T0 = 0.1
df10 <- extractLog("log/190528hindhelogT10.txt")

ggplot(df10, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Final_NM - Init_NM)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))

ggplot(df10[df10$InitMax_HindHe > 0.5,]) +
  geom_density(aes(x = FinalMax_HindHe), col = "blue") +
  geom_density(aes(x = InitMax_HindHe), col = "red") +
  geom_vline(xintercept = 0.5) 

# compare quantitatively
median(df50$FinalMax_HindHe)
median(df25$FinalMax_HindHe)
median(df10$FinalMax_HindHe)
mean(df50$FinalMax_HindHe)
mean(df25$FinalMax_HindHe)
mean(df10$FinalMax_HindHe)
sd(df50$FinalMax_HindHe)
sd(df25$FinalMax_HindHe)
sd(df10$FinalMax_HindHe)

t.test(df50$FinalMax_HindHe, df10$FinalMax_HindHe)

hist(log(df10$Final_NM / df10$Init_NM))
median(df50$Final_NM / df50$Init_NM)
median(df25$Final_NM / df25$Init_NM)
median(df10$Final_NM / df10$Init_NM)

plot(df50$Init_NM, df25$Init_NM)
