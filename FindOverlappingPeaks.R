#Micah Lessnick
#peak processing with bioconductor

####Creating GRange Object####
library(GenomicRanges)
library(tidyr)
library(ggplot2)
library(readxl)
library(dplyr)
library(purrr)

csv_A673 <- read_excel("C:/Users/micah/OneDrive/Desktop/Ford Project/Excel project data.xlsx",
sheet = "peaks_A673_FLI_fdr05_idr05_mean")

csv_EWS <- read_excel("C:/Users/micah/OneDrive/Desktop/Ford Project/Excel project data.xlsx",
sheet = "peaks_EWS_502_FLI_fdr05_idr05_m")

csv_SKNMC <- read_excel("C:/Users/micah/OneDrive/Desktop/Ford Project/Excel project data.xlsx",
sheet = "peaks_SK_N_MC_FLI_fdr05_idr05_m")

csv_TC71 <- read_excel("C:/Users/micah/OneDrive/Desktop/Ford Project/Excel project data.xlsx", 
sheet = "peaks_TC71_FLI_fdr05_idr05_mean")

csv_TTC <- read_excel("C:/Users/micah/OneDrive/Desktop/Ford Project/Excel project data.xlsx", 
sheet = "peaks_TTC_466_ERG_fdr05_idr05_m")

#this method should take in a csv and convert it to a grange object
#the csv must have cromosome number, the strand, start, and end locations
  #all other data becomes metadata for the peak and appears on the rhs of the obj
  #in this case, strand is * which means the data is coming from both strands
convToGRange <- function(csv){
  #seqnames refers to the chromosome number
  return(makeGRangesFromDataFrame(csv,
                                  keep.extra.columns = TRUE,
                                  seqnames.field = ("#chrom")))
}

gr_A673 <- convToGRange(csv_A673)
gr_EWS <- convToGRange(csv_EWS)
gr_SKNMC <- convToGRange(csv_SKNMC)
gr_TC71 <- convToGRange(csv_TC71)
gr_TTC <- convToGRange(csv_TTC)

####Finding overlaps####

#subsetByOverlaps can only work on 2 gr objects at a time
  #allows for the invert flag (search for peaks that do NOT overlap)
  #returns full grange object, not hit list (list of indexes that overlap)
oLapFunc <- function(grQ, grS){
  return(subsetByOverlaps(grQ, grS))
}

#Leave TTC out of overlaps for some extra test cases
#Compute overlaps for all other grange objects
oLapFull <- oLapFunc(oLapFunc(oLapFunc(gr_A673, gr_EWS), gr_SKNMC), gr_TC71)

####Finding unique peaks####

#same process as finding overlaps, just set invert flag to true in subsetByOverlaps
noLapFunc <- function(grQ, grS){
  return(subsetByOverlaps(grQ, grS, invert = TRUE))
}

noLapFull <- noLapFunc(noLapFunc(noLapFunc(gr_A673, gr_EWS), gr_SKNMC), gr_TC71)

####Convert granges to dataframes####

oLapDF <- as.data.frame(oLapFull)

noLapDF <- as.data.frame(noLapFull)

####Finish data processing####

#add class column to dataframes
oLapDF <- cbind.data.frame(oLapDF, class = 1)
noLapDF <- cbind.data.frame(noLapDF, class = 0)

#remove unneccessary data columns (start and end, strand, name)
oLapDF <- subset(oLapDF, select = -c(start, end, strand, name))
noLapDF <- subset(noLapDF, select = -c(start, end, strand, name))

#combine into single dataframe
fullDF <- rbind.data.frame(oLapDF, noLapDF)

####Export to csv####
write.csv(fullDF, "C:\\Users\\micah\\OneDrive\\Desktop\\Ford Project\\peakData.csv", row.names = FALSE)

####Generate file for other cell line metadata####
#EWS metadata
#compute overlap
oLapFullEws <- oLapFunc(oLapFunc(oLapFunc(gr_EWS, gr_SKNMC), gr_TC71), gr_A673)
noLapFullEws <- noLapFunc(noLapFunc(noLapFunc(gr_EWS, gr_SKNMC), gr_TC71), gr_A673)

#remove extraneous columns and label classes
oLapEwsDF <- as.data.frame(oLapFullEws) %>%
  cbind.data.frame(class = 1) %>%
  subset(select = -c(start, end, strand, name, seqnames))
noLapEwsDF <- as.data.frame(noLapFullEws) %>%
  cbind.data.frame(class = 0) %>%
  subset(select = -c(start, end, strand, name, seqnames))

#combine into one dataframe
fullEwsDF <- rbind.data.frame(oLapEwsDF, noLapEwsDF)

#export to csv
write.csv(fullEwsDF, "C:\\Users\\micah\\OneDrive\\Desktop\\Ford Project\\peakDataEws.csv", row.names = FALSE)

#SKNMC metadata
#compute overlap
oLapFullSknmc <- oLapFunc(oLapFunc(oLapFunc(gr_SKNMC, gr_TC71), gr_A673), gr_EWS)
noLapFullSknmc <- noLapFunc(noLapFunc(noLapFunc(gr_SKNMC, gr_TC71), gr_A673), gr_EWS)

#remove extraneous columns and label classes
oLapSknmcDF <- as.data.frame(oLapFullSknmc) %>%
  cbind.data.frame(class = 1) %>%
  subset(select = -c(start, end, strand, name, seqnames))
noLapSknmcDF <- as.data.frame(noLapFullSknmc) %>%
  cbind.data.frame(class = 0) %>%
  subset(select = -c(start, end, strand, name, seqnames))

#combine into one dataframe
fullSknmcDF <- rbind.data.frame(oLapSknmcDF, noLapSknmcDF)

#export to csv
write.csv(fullSknmcDF, "C:\\Users\\micah\\OneDrive\\Desktop\\Ford Project\\peakDataSknmc.csv", row.names = FALSE)

#TC71 metadata
#compute overlap
oLapFullTc71 <- oLapFunc(oLapFunc(oLapFunc(gr_TC71, gr_A673), gr_EWS), gr_SKNMC)
noLapFullTc71 <- noLapFunc(noLapFunc(noLapFunc(gr_TC71, gr_A673), gr_EWS), gr_SKNMC)

#remove extraneous columns and label classes
oLapTc71DF <- as.data.frame(oLapFullTc71) %>%
  cbind.data.frame(class = 1) %>%
  subset(select = -c(start, end, strand, name, seqnames))
noLapTc71DF <- as.data.frame(noLapFullTc71) %>%
  cbind.data.frame(class = 0) %>%
  subset(select = -c(start, end, strand, name, seqnames))

#combine into one dataframe
fullTc71DF <- rbind.data.frame(oLapTc71DF, noLapTc71DF)

#export to csv
write.csv(fullTc71DF, "C:\\Users\\micah\\OneDrive\\Desktop\\Ford Project\\peakDataTc71.csv", row.names = FALSE)

####Data exploration of cell lines####
#distribution of log2FC_score for the entire df of the cell lines
par(mfrow=c(2,2))
hist(fullDF$log2FC_score)
hist(fullEwsDF$log2FC_score)
hist(fullSknmcDF$log2FC_score)
hist(fullTc71DF$log2FC_score)

#distribution of log2FC for the entire df of the cell lines
par(mfrow=c(2,2))
hist(fullDF$log2FC)
hist(fullEwsDF$log2FC)
hist(fullSknmcDF$log2FC)
hist(fullTc71DF$log2FC)

#distribution of x.log10pval for the entire df of the cell lines
par(mfrow=c(2,2))
hist(fullDF$X.log10pval)
hist(fullEwsDF$X.log10pval)
hist(fullSknmcDF$X.log10pval)
hist(fullTc71DF$X.log10pval)

#distribution of X.log10FDR for the entire df of the cell lines
par(mfrow=c(2,2))
hist(fullDF$X.log10FDR)
hist(fullEwsDF$X.log10FDR)
hist(fullSknmcDF$X.log10FDR)
hist(fullTc71DF$X.log10FDR)

#distribution of all features for A673
fullDF %>%
  subset(select = -c(seqnames, width, class)) %>%
  gather() %>%
  ggplot(aes(value)) +
    facet_wrap(~ key, scales = "free") +
    geom_density()

#distribution of all features for EWS
fullEwsDF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

#distribution of all features for SKNMC
fullSknmcDF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

#distribution of all features for TC71
fullTc71DF %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_density()

#statistical info about each cell line
fullDF %>%
  subset(select = -c(seqnames)) %>%
  summarise_each(funs(mean, sd))

fullEwsDF %>%
  summarise_each(funs(mean, sd))

fullSknmcDF %>%
  summarise_each(funs(mean, sd))

fullTc71DF %>%
  summarise_each(funs(mean, sd))

####combine and normalize all computed cell line dfs####
fullA673DF <- fullDF %>% 
  subset(select = -c(seqnames))

#zscore normalization function for dataframes
zscoreDFFunc <- function(df){
  dfMean <- colMeans(df)
  dfSd <- apply(df, 2, sd)
  dfTemp <- df
  for (cols in 1:(length(dfMean)-1)){
    df[, cols] <- ((dfTemp[, cols]-dfMean[cols])/dfSd[cols])
  }
  print(head(df))
  return(df)
}

#combine first then normalize:
allCellCNPreN <- rbind.data.frame(fullA673DF, fullEwsDF) %>%
  rbind.data.frame(fullSknmcDF) %>%
  rbind.data.frame(fullTc71DF)

allCellCN <- zscoreDFFunc(allCellCNPreN)
write.csv(allCellCN, "C:\\Users\\micah\\OneDrive\\Desktop\\Ford Project\\peakDataAllCN.csv", row.names = FALSE)

#normalize first then combine
normA673 <- zscoreDFFunc(fullA673DF)
normEws <- zscoreDFFunc(fullEwsDF)
normSknmc <- zscoreDFFunc(fullSknmcDF)
normTc71 <- zscoreDFFunc(fullTc71DF)

allCellNC <- rbind.data.frame(normA673, normEws)%>%
  rbind.data.frame(normSknmc)%>%
  rbind.data.frame(normTc71)

write.csv(allCellNC, "C:\\Users\\micah\\OneDrive\\Desktop\\Ford Project\\peakDataAllNC.csv", row.names = FALSE)

####Continued Data Exploration####
#sample 1k peaks from full dataset (normalized then combined)
set.seed(17)
sampleNC <- slice_sample(allCellNC, n = 1000)

#convert class to character (as opposed to Boolean)
sampleNC$class <- as.character(sampleNC$class)

#colorblind-friendly color pallet
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#graph scatterplot of all points based on FC and FDR, coloring based on class
#FC_score and pval have identical data due to their only difference being scale (which was removed by normalization)
ggplot(sampleNC, aes(x = X.log10FDR, y = log2FC, color = class)) +
  geom_point() +
  geom_smooth(se = FALSE, aes(group = class)) +
  scale_color_manual(values = cbPalette)

#histogram of sample set
ggplot(sampleNC, aes(x = log2FC, fill = class, color = class)) +
  geom_histogram()

####Label 5th cell line####
#compute overlap
oLapFullTtc <- oLapFunc(oLapFunc(oLapFunc(gr_TTC, gr_A673), gr_EWS), gr_SKNMC)
noLapFullTtc <- noLapFunc(noLapFunc(noLapFunc(gr_TTC, gr_A673), gr_EWS), gr_SKNMC)

#remove extraneous columns and label classes
oLapTtcDF <- as.data.frame(oLapFullTtc) %>%
  cbind.data.frame(class = 1) %>%
  subset(select = -c(start, end, strand, name, seqnames))
noLapTtcDF <- as.data.frame(noLapFullTc71) %>%
  cbind.data.frame(class = 0) %>%
  subset(select = -c(start, end, strand, name, seqnames))

#combine into one dataframe
fullTtcDF <- rbind.data.frame(oLapTtcDF, noLapTtcDF)

#normalize dataset
normFullTtcDF <- zscoreDFFunc(fullTtcDF)

#export to csv
write.csv(normFullTtcDF, "C:\\Users\\micah\\OneDrive\\Desktop\\Ford Project\\peakDataTtc.csv", row.names = FALSE)


####Plotting peaks without a class####
#combine all raw datasets into 1 df
rawDf <- as.data.frame(csv_A673) %>%
  rbind.data.frame(as.data.frame(csv_EWS)) %>%
  rbind.data.frame(as.data.frame(csv_SKNMC)) %>%
  rbind.data.frame(as.data.frame(csv_TC71))

#setdiff between full set and raw data by (log2FC, log2FC_score, X.log10pval, X.log10FDR)
#undeterminedPeaks <- setdiff()

#add width and remove excess features from rawDF
#add 1 to width to account for index seen in overlap function
rawDf[, "width"] = (rawDf$end-rawDf$start + 1)

#rename non-matching feature names
rawDf <- rename(rawDf, X.log10pval = "-log10pval") %>%
  rename(X.log10FDR = "-log10FDR")

#remove excess features (because column 1 starts with '#', must remove by index)
rawDf <- rawDf %>%
  subset(select = -c(start, end, strand, name)) %>%
  subset(select = -c(1))

#reorder columns to mirror overlapped data
colOrder <- c("width", "log2FC_score", "log2FC", "X.log10pval", "X.log10FDR")
rawDf <- rawDf[, colOrder]

#compute left join
undetPeaks <- left_join(rawDf, allCellCNPreN, by = c("log2FC_score", "log2FC", "X.log10pval", "X.log10FDR"))

#remove labeled peaks and label and remaining peaks as 2 (for undetermined)
undetPeaks <- filter(undetPeaks, is.na(class)) %>%
  select(-width.y) %>%
  mutate(class = 2)

#rename remaining width.x to width
undetPeaks <- rename(undetPeaks, width = width.x)

#normalize df
undetPeaks <- zscoreDFFunc(undetPeaks)

#plot first 1k of undetermined peaks and 1k determined peaks
set.seed(17)
sampleNC <- slice_sample(allCellNC, n = 1000)
sampleUP <- slice_sample(undetPeaks, n = 1000)

#change class to character so plot will not use gradient color
sampleNC$class <- as.character(sampleNC$class)
sampleUP$class <- as.character(sampleUP$class)

rbind(sampleNC, sampleUP) %>%
  ggplot(aes(x = X.log10FDR, y = log2FC, color = class)) +
  geom_point() +
  geom_smooth(se = FALSE, aes(group = class)) +
  scale_color_manual(values = cbPalette)
