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
  
  markerLines <- initHindHeLines - 2
  initNMLines <- initHindHeLines + 1
  finalNMLines <- initHindHeLines + 3
  unsplitHindHeLines <- initHindHeLines - 1
  if(finalTemp){
    finalTempLines <- initHindHeLines + 4
  }
  if(tabu){
    repLines <- initHindHeLines + 4
  }
  
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
  if(finalTemp){
    outdf$Final_Temp <- as.numeric(sub("Final temperature: ", "", logout[finalTempLines]))
  }
  if(tabu){
    outdf$Reps <- as.integer(sub("Rep where best solution found: ", "", logout[repLines]))
  }
  
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

# try with T0 = 0.05
df05 <- extractLog("log/190528hindhelogT05.txt")

ggplot(df05, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Final_NM - Init_NM)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))

ggplot(df05[df05$InitMax_HindHe > 0.5,]) +
  geom_density(aes(x = FinalMax_HindHe), col = "blue") +
  geom_density(aes(x = InitMax_HindHe), col = "red") +
  geom_vline(xintercept = 0.5)

# compare quantitatively
median(df50$FinalMax_HindHe)
median(df25$FinalMax_HindHe)
median(df10$FinalMax_HindHe)
median(df05$FinalMax_HindHe)
mean(df50$FinalMax_HindHe)
mean(df25$FinalMax_HindHe)
mean(df10$FinalMax_HindHe)
mean(df05$FinalMax_HindHe)
sd(df50$FinalMax_HindHe)
sd(df25$FinalMax_HindHe)
sd(df10$FinalMax_HindHe)
sd(df05$FinalMax_HindHe)

sd(df10$FinalMax_HindHe[df10$InitMax_HindHe > 0.5])
sd(df05$FinalMax_HindHe[df05$InitMax_HindHe > 0.5])

t.test(df50$FinalMax_HindHe, df05$FinalMax_HindHe)

hist(log(df10$Final_NM / df10$Init_NM))
median(df50$Final_NM / df50$Init_NM)
median(df25$Final_NM / df25$Init_NM)
median(df10$Final_NM / df10$Init_NM)

plot(df50$Init_NM, df25$Init_NM)

mean(df50$FinalMax_HindHe > 0.4 & df50$FinalMax_HindHe < 0.6)
mean(df25$FinalMax_HindHe > 0.4 & df25$FinalMax_HindHe < 0.6)
mean(df10$FinalMax_HindHe > 0.4 & df10$FinalMax_HindHe < 0.6)
mean(df05$FinalMax_HindHe > 0.4 & df05$FinalMax_HindHe < 0.6)

# quantitative differences are not big, but T0 = 0.10 seems to best bring values to around 0.5

# Re-run incorporating final temperature, and seeing if reps per temp changes things
df10 <- extractLog("log/190529hindhelogT10.txt", finalTemp = TRUE)
df10$Num_Temps <- log(df10$Final_Temp/0.1, base = 0.95) # number of temperatures attempted
head(df10)
summary(df10$Num_Temps) #(NM is wrong in these two)

ggplot(df10, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Num_Temps)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))
# see that more effort has been put into the ones that didn't get to the desired value

# five times as many attempted swaps per temperature
df10_5x <- extractLog("log/190529hindhelogT10_5xSwaps.txt", finalTemp = TRUE)
df10_5x$Num_Temps <- log(df10_5x$Final_Temp/0.1, base = 0.95)
# There are 47 in this one and 44 in the last one, I assume from random chance of 
# how haplotypes were initially assigned.

ggplot(df10_5x, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Num_Temps)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))

median(df10_5x$Num_Temps) # 59 - took longer to converge (probably b/c more chances for a swap)
median(df10$Num_Temps)    # 44
median(df10_5x$Final_Temp) # 0.005
median(df10$Final_Temp)    # 0.010

mean(df10$FinalMax_HindHe < 0.6)    # 50%
mean(df10_5x$FinalMax_HindHe < 0.6) # 66% - slightly better but should try with more markers.

# how many haplotypes per group do we typically have?
twoalignChr1Chr2 <- read.table("marker_CSV/190525twoalign_Chr1Chr2.csv", sep = ",", header = FALSE)
table(table(paste(twoalignChr1Chr2$V1, twoalignChr1Chr2$V2)))
# generally few, but up to over 100
log2(1e7)

# when a record of the best solution is kept
df10_5x_keepbest <- extractLog("log/190531hindhelogT10_5xSwaps.txt", finalTemp = TRUE)
df10_5x_keepbest$Num_Temps <- log(df10_5x_keepbest$Final_Temp/0.1, base = 0.95)

ggplot(df10_5x_keepbest, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Num_Temps)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))

identical(rownames(df10_5x), rownames(df10_5x_keepbest))
mean(rownames(df10_5x) %in% rownames(df10_5x_keepbest))
commonsites <- rownames(df10_5x)[rownames(df10_5x) %in% rownames(df10_5x_keepbest)]

plot(df10_5x[commonsites,"FinalMax_HindHe"],
     df10_5x_keepbest[commonsites,"FinalMax_HindHe"])
abline(b = 1, a = 0, col = "red")

# the actual value we are trying to optimize
init_excess_hindhe <- as.matrix(df10_5x_keepbest[,c("Init1_HindHe", "Init2_HindHe")]) - 0.5
init_excess_hindhe[init_excess_hindhe < 0] <- 0
df10_5x_keepbest$InitExcess_HindHe <- rowMeans(init_excess_hindhe, na.rm = TRUE)

final_excess_hindhe <- as.matrix(df10_5x_keepbest[,c("Final1_HindHe", "Final2_HindHe")]) - 0.5
final_excess_hindhe[final_excess_hindhe < 0] <- 0
df10_5x_keepbest$FinalExcess_HindHe <- rowMeans(final_excess_hindhe, na.rm = TRUE)

init_excess_hindhe <- as.matrix(df10_5x[,c("Init1_HindHe", "Init2_HindHe")]) - 0.5
init_excess_hindhe[init_excess_hindhe < 0] <- 0
df10_5x$InitExcess_HindHe <- rowMeans(init_excess_hindhe, na.rm = TRUE)

final_excess_hindhe <- as.matrix(df10_5x[,c("Final1_HindHe", "Final2_HindHe")]) - 0.5
final_excess_hindhe[final_excess_hindhe < 0] <- 0
df10_5x$FinalExcess_HindHe <- rowMeans(final_excess_hindhe, na.rm = TRUE)

plot(df10_5x_keepbest$InitExcess_HindHe, df10_5x_keepbest$FinalExcess_HindHe)
abline(b = 1, a = 0, col = "red")
grid()

plot(df10_5x$InitExcess_HindHe, df10_5x$FinalExcess_HindHe)
abline(b = 1, a = 0, col = "red")
grid()

plot(df10_5x[commonsites,"FinalExcess_HindHe"],
     df10_5x_keepbest[commonsites,"FinalExcess_HindHe"])
abline(b = 1, a = 0, col = "red")

plot(df10_5x[commonsites,"Final_Temp"],
     df10_5x_keepbest[commonsites,"Final_Temp"])
abline(b = 1, a = 0, col = "red")

# number of haplotypes vs. whether a group goes thru simulated annealing?
nhap_by_group <- table(paste(twoalignChr1Chr2$V1[1:5000], twoalignChr1Chr2$V2[1:5000]))
hist(nhap_by_group)
length(nhap_by_group)
did_anneal = names(nhap_by_group) %in% rownames(df10_5x_keepbest)
summary(as.integer(nhap_by_group[did_anneal]))
summary(as.integer(nhap_by_group[!did_anneal]))
# tends to be ones where there were ~30 haplotypes

df10_5x_keepbest$nHap <- as.integer(nhap_by_group[rownames(df10_5x_keepbest)])
ggplot(df10_5x_keepbest, aes(y = FinalExcess_HindHe, x = nHap, color = InitExcess_HindHe)) +
  geom_point() +
  scale_color_viridis()

# Adding allele correlations in to the algorithm ####
# this first version had a mistake but let's take a look anyway
dfCorr <- extractLog("log/190603hindhelogT10_5xSwaps_UseCorr.txt", finalTemp = TRUE)
head(dfCorr)
dfCorr$Num_Temps <- log(dfCorr$Final_Temp/0.1, base = 0.95)

ggplot(dfCorr, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Num_Temps)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))

init_excess_hindhe <- as.matrix(dfCorr[,c("Init1_HindHe", "Init2_HindHe")]) - 0.5
init_excess_hindhe[init_excess_hindhe < 0] <- 0
dfCorr$InitExcess_HindHe <- rowMeans(init_excess_hindhe, na.rm = TRUE)

final_excess_hindhe <- as.matrix(dfCorr[,c("Final1_HindHe", "Final2_HindHe")]) - 0.5
final_excess_hindhe[final_excess_hindhe < 0] <- 0
dfCorr$FinalExcess_HindHe <- rowMeans(final_excess_hindhe, na.rm = TRUE)

plot(dfCorr$InitExcess_HindHe, dfCorr$FinalExcess_HindHe)
abline(b = 1, a = 0, col = "red")
grid()

mean(dfCorr$FinalExcess_HindHe == 0)           # 49%
mean(df10_5x_keepbest$FinalExcess_HindHe == 0) # 49%
mean(dfCorr$Num_Temps)           # 63
mean(df10_5x_keepbest$Num_Temps) # 63

# after fixing how correlations estimated
dfCorr1 <- extractLog("log/190604hindhelogT10_5xSwaps_UseCorr.txt", finalTemp = TRUE)
head(dfCorr)
dfCorr1$Num_Temps <- log(dfCorr1$Final_Temp/0.1, base = 0.95)

ggplot(dfCorr1, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Num_Temps)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis() +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1.25))

init_excess_hindhe <- as.matrix(dfCorr1[,c("Init1_HindHe", "Init2_HindHe")]) - 0.5
init_excess_hindhe[init_excess_hindhe < 0] <- 0
dfCorr1$InitExcess_HindHe <- rowMeans(init_excess_hindhe, na.rm = TRUE)

final_excess_hindhe <- as.matrix(dfCorr1[,c("Final1_HindHe", "Final2_HindHe")]) - 0.5
final_excess_hindhe[final_excess_hindhe < 0] <- 0
dfCorr1$FinalExcess_HindHe <- rowMeans(final_excess_hindhe, na.rm = TRUE)

mean(dfCorr1$FinalExcess_HindHe == 0)  # 48%
mean(dfCorr1$Num_Temps)                # 64

plot(dfCorr1$InitExcess_HindHe, dfCorr1$FinalExcess_HindHe)
abline(b = 1, a = 0, col = "red")

# can we speed up run time now?
dfCorr2 <- extractLog("log/190604hindhelogT10_rho90_2xSwaps_UseCorr.txt", finalTemp = TRUE)
dfCorr2$Num_Temps <- log(dfCorr2$Final_Temp/0.1, base = 0.90)

init_excess_hindhe <- as.matrix(dfCorr2[,c("Init1_HindHe", "Init2_HindHe")]) - 0.5
init_excess_hindhe[init_excess_hindhe < 0] <- 0
dfCorr2$InitExcess_HindHe <- rowMeans(init_excess_hindhe, na.rm = TRUE)

final_excess_hindhe <- as.matrix(dfCorr2[,c("Final1_HindHe", "Final2_HindHe")]) - 0.5
final_excess_hindhe[final_excess_hindhe < 0] <- 0
dfCorr2$FinalExcess_HindHe <- rowMeans(final_excess_hindhe, na.rm = TRUE)

mean(dfCorr2$FinalExcess_HindHe == 0)  # 42%
mean(dfCorr2$Num_Temps)                # 31

prop.test(c(sum(dfCorr1$FinalExcess_HindHe == 0), sum(dfCorr2$FinalExcess_HindHe == 0)),
          c(nrow(dfCorr1), nrow(dfCorr2))) # need a bigger sample to know if these are really different

# is it generally the same markers that get the best solution?
Corr1solved <- rownames(dfCorr1)[dfCorr1$FinalExcess_HindHe == 0]
Corr2solved <- rownames(dfCorr2)[dfCorr2$FinalExcess_HindHe == 0]
mean(Corr1solved %in% Corr2solved) # 82%
commonmarkers <- rownames(dfCorr1)[rownames(dfCorr1) %in% rownames(dfCorr2)] # 44

plot(dfCorr1[commonmarkers, "FinalExcess_HindHe"], dfCorr2[commonmarkers, "FinalExcess_HindHe"],
     col = "#00000040")
abline(a = 0, b = 1, col = "red")
# if it was more than zero, it was generally more than zero for both
# quicker algorithm maybe slightly worse, but again we need more data points

# Examine output of tabu search ####
dfTabu <- extractLog("log/190605tabu_log.txt", tabu = TRUE)

init_excess_hindhe <- as.matrix(dfTabu[,c("Init1_HindHe", "Init2_HindHe")]) - 0.5
init_excess_hindhe[init_excess_hindhe < 0] <- 0
dfTabu$InitExcess_HindHe <- rowMeans(init_excess_hindhe, na.rm = TRUE)

final_excess_hindhe <- as.matrix(dfTabu[,c("Final1_HindHe", "Final2_HindHe")]) - 0.5
final_excess_hindhe[final_excess_hindhe < 0] <- 0
dfTabu$FinalExcess_HindHe <- rowMeans(final_excess_hindhe, na.rm = TRUE)

ggplot(dfTabu, aes(x = InitMax_HindHe, y = FinalMax_HindHe, col = Reps)) +
  geom_point() +
  geom_hline(yintercept = 1/2) +
  geom_vline(xintercept = 1/2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis()

ggplot(dfTabu, aes(x = InitExcess_HindHe, y = FinalExcess_HindHe, col = Reps)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis()

mean(dfTabu$FinalExcess_HindHe == 0)  # 70%
table(dfTabu$Reps[dfTabu$FinalExcess_HindHe == 0])
table(dfTabu$Reps[dfTabu$FinalExcess_HindHe < 0.05])
# quite a lot that are zero reps, which I assume means that just using correlations
# worked. don't do tabu if ok after correlations.
# probably 25 reps is an ok cutoff.

mean(dfTabu$FinalExcess_HindHe[dfTabu$Reps > 0] == 0) # 48% if it actually searched
