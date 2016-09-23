# Part 1: Reading and analyzing beahvior and physiology data

## load libraries -----
library(dplyr) ## for filtering and selecting rows
library(plyr) ## for renmaing factors
library(ggplot2) ## for awesome plots!
library(reshape2) #@ for melting dataframe

## wrangle the raw wt and fmr1 dataframes ----

## read the data 
setwd("~/Github/qPCR-mouse/Rayna/2013/data")
wt <- read.csv("01_WTbehavephys.csv", header=TRUE, stringsAsFactors = FALSE, na.strings = c("", "ND", "N/A"))
fmr1 <- read.csv("01_FMR1behavephys.csv", header=TRUE, stringsAsFactors = FALSE, na.strings = c("", "ND", "N/A"))

## fixing columns both fmr1 and wt are as equivilent as possible
wt$filename <- as.character(paste(wt$filename, wt$X, sep="_")) ## mgerge 2 columns 

## delete obsolete columns
wt$X <- NULL 
wt$X.1 <- NULL
fmr1$X <- NULL
wt$next. <- NULL
fmr1$next. <- NULL
fmr1$PopSpike.max........mV. <- NA  ## add column
fmr1$PospSpike.Probability <- NA  ## add column
fmr1 <- fmr1[-c(75:80), ] ## removing rows
names(wt)
names(fmr1)


## create ind column with a animal name that matches qpcr (e.g "FMR1 AB" or "BL 19")
wt$ind <- wt$filename
head(wt$ind)
wt$ind <- gsub("[[:blank:]]*bl", "BL", wt$ind) ##removed blank space then changes bl1... to BL 1...
wt$ind <- gsub("D[[:print:]]*", "", wt$ind) ## deletes the D (for Day) and everythign thereafter
wt$ind

fmr1$ind <- fmr1$filename
fmr1$ind <- gsub("RoomTrack_", "", fmr1$ind) ##remove RommTrack_"
fmr1$ind <- gsub("fmr1", "FMR1", fmr1$ind) ## changes to FMR1
fmr1$ind <- gsub("frm1", "FMR1", fmr1$ind) ## because one was typed wrong
fmr1$ind <- gsub("wildtype", "BL", fmr1$ind) ## that one wildtype sample name
fmr1$ind <- gsub("_[ptr][[:print:]]*", "", fmr1$ind) ##deletes the rest of the filename
fmr1$ind

## rename columns
names(wt)[3] <- "TotalTime"
names(wt)[4] <- "TotalPath"
names(wt)[5] <- "Entrances"
names(wt)[6] <- "TimeToFirstEntrance"
names(wt)[7] <- "PathToFirstEntrance"
names(wt)[8] <- "SpeedToFirstEntrance"
names(wt)[9] <- "EntrancePerDistance"
names(wt)[10] <- "TotalShocks"
names(wt)[11] <- "TimeToFirstShock"
names(wt)[12] <- "PathToFirstShock"
names(wt)[13] <- "Speed"
names(wt)[14] <- "SDofSpeed"
names(wt)[15] <- "LineArity"
names(wt)[16] <- "MaxAvoidTime"
names(wt)[17] <- "MaxAvoidPath"
names(wt)[18] <- "TimeToSecondEntrance"
names(wt)[19] <- "PathToSecondEntrance"
names(wt)[20] <- "SpeedToSecondEntranc"
names(wt)[21] <- "TimeTARG"
names(wt)[22] <- "pTimeTARG"
names(wt)[23] <- "pTimeCCW"
names(wt)[24] <- "pTimeOPP"
names(wt)[25] <- "pTimeCW"
names(wt)[26] <- "RayleigLength"
names(wt)[27] <- "RayleigAngle"
names(wt)[28] <- "PolarAveVal"
names(wt)[29] <- "PolarSdVal"
names(wt)[30] <- "PolarMinVal"
names(wt)[31] <- "PolarMinBin"
names(wt)[32] <- "MinLoBin"
names(wt)[33] <- "MinHiBin"
names(wt)[34] <- "PolarMaxVal"
names(wt)[35] <- "PolarMaxBin"
names(wt)[36] <- "MaxLoBin"
names(wt)[37] <- "MaxHiBin"
names(wt)[38] <- "AnnularMinVal"
names(wt)[39] <- "AnnularMinBin"
names(wt)[40] <- "AnnularMaxVal"
names(wt)[41] <- "AnnularMaxBin"
names(wt)[42] <- "AnnularAvg"
names(wt)[43] <- "AnnularSD"
names(wt)[44] <- "AnnularSkewnes"
names(wt)[45] <- "AnnularKurtosis"
names(wt)[46] <- "ind_bad"
names(wt)[47] <- "genotype_bad"
names(wt)[48] <- "APA"
names(wt)[49] <- "TotalLocomotorActivity"
names(wt)[50] <- "LastLocomotorActivity"
names(wt)[51] <- "TotalPunishment"
names(wt)[53] <- "LastPunishment"
names(wt)[54] <- "LastTotalEntrances"
names(wt)[55] <- "T1Retention"
names(wt)[56] <- "T2Retention"
names(wt)[57] <- "IO_Max"
names(wt)[58] <- "LTP_Baseline"
names(wt)[59] <- "LTP_Baseline_SD"
names(wt) # check all good

## do same for fmr1
names(fmr1)[3] <- "TotalTime"
names(fmr1)[4] <- "TotalPath"
names(fmr1)[5] <- "Entrances"
names(fmr1)[6] <- "TimeToFirstEntrance"
names(fmr1)[7] <- "PathToFirstEntrance"
names(fmr1)[8] <- "SpeedToFirstEntrance"
names(fmr1)[9] <- "EntrancePerDistance"
names(fmr1)[10] <- "TotalShocks"
names(fmr1)[11] <- "TimeToFirstShock"
names(fmr1)[12] <- "PathToFirstShock"
names(fmr1)[13] <- "Speed"
names(fmr1)[14] <- "SDofSpeed"
names(fmr1)[15] <- "LineArity"
names(fmr1)[16] <- "MaxAvoidTime"
names(fmr1)[17] <- "MaxAvoidPath"
names(fmr1)[18] <- "TimeToSecondEntrance"
names(fmr1)[19] <- "PathToSecondEntrance"
names(fmr1)[20] <- "SpeedToSecondEntranc"
names(fmr1)[21] <- "TimeTARG"
names(fmr1)[22] <- "pTimeTARG"
names(fmr1)[23] <- "pTimeCCW"
names(fmr1)[24] <- "pTimeOPP"
names(fmr1)[25] <- "pTimeCW"
names(fmr1)[26] <- "RayleigLength"
names(fmr1)[27] <- "RayleigAngle"
names(fmr1)[28] <- "PolarAveVal"
names(fmr1)[29] <- "PolarSdVal"
names(fmr1)[30] <- "PolarMinVal"
names(fmr1)[31] <- "PolarMinBin"
names(fmr1)[32] <- "MinLoBin"
names(fmr1)[33] <- "MinHiBin"
names(fmr1)[34] <- "PolarMaxVal"
names(fmr1)[35] <- "PolarMaxBin"
names(fmr1)[36] <- "MaxLoBin"
names(fmr1)[37] <- "MaxHiBin"
names(fmr1)[38] <- "AnnularMinVal"
names(fmr1)[39] <- "AnnularMinBin"
names(fmr1)[40] <- "AnnularMaxVal"
names(fmr1)[41] <- "AnnularMaxBin"
names(fmr1)[42] <- "AnnularAvg"
names(fmr1)[43] <- "AnnularSD"
names(fmr1)[44] <- "AnnularSkewnes"
names(fmr1)[45] <- "AnnularKurtosis"
names(fmr1)[46] <- "ind_bad"
names(fmr1)[47] <- "genotype_bad"
names(fmr1)[48] <- "APA"
names(fmr1)[49] <- "TotalLocomotorActivity"
names(fmr1)[50] <- "LastLocomotorActivity"
names(fmr1)[51] <- "TotalPunishment"
names(fmr1)[53] <- "LastPunishment"
names(fmr1)[54] <- "LastTotalEntrances"
names(fmr1)[55] <- "T1Retention"
names(fmr1)[56] <- "T2Retention"
names(fmr1)[57] <- "IO_Max"
names(fmr1)[58] <- "LTP_Baseline"
names(fmr1)[59] <- "LTP_Baseline_SD"

## add columns genotype (easier for wt df than for fmr1 df)
wt$genotype <- as.factor("WT")
fmr1$genotype <- ifelse(grepl("fmr1|frm1",fmr1$filename),'FMR1-KO','WT')
fmr1$genotype <- as.factor(fmr1$genotype)

## binding, clearning, filtering, renaming, adding conditions ----

wtfmr1 <- rbind(wt, fmr1) ## combine the datasets into one
wtfmr1 <- filter(wtfmr1, grepl("Room", filename)) ## remove non-data rows
str(wtfmr1)

## making strings factors
wtfmr1$ind <- as.factor(wtfmr1$ind)  
wtfmr1$APA <- as.factor(wtfmr1$APA)  

## making a bunch of columns numeric
wtfmr1[, c(3:35)] <- sapply(wtfmr1[, c(3:35)], as.numeric)
wtfmr1[, c(36:49)] <- sapply(wtfmr1[, c(36:49)], as.numeric)
str(wtfmr1)

## add column session with APA information
wtfmr1$session <- ifelse(grepl("pretraining|pretrain|Hab", wtfmr1$filename), "pretraining", 
                     ifelse(grepl("training1|Train1", wtfmr1$filename), "training1",
                            ifelse(grepl("training2|Train2", wtfmr1$filename), "training2",
                                   ifelse(grepl("training3|Train3", wtfmr1$filename), "training3",
                                          ifelse(grepl("retention|reten|Retest", wtfmr1$filename), "retention", "NA")))))
wtfmr1$session  ## check that all names good with no NAs                                 
wtfmr1$session <- as.factor(wtfmr1$session)  
wtfmr1$session <- factor(wtfmr1$session, levels = c("pretraining", "training1", "training2", "training3", "retention"))

### separate out the summary columns and clean up / wrangle ----
summary <- filter(wtfmr1, session == "retention")
summary <- summary[c(60:61,48,62,49:59,1)] #select and reorder columns, so animal first
names(summary)
str(summary)

## rename APA to match qpcr data
summary$APA <- revalue(summary$APA, c("untrained" = "control")) 
summary$APA <- factor(summary$APA, levels = c("control", "trained"))

## making a bunch of columns numeric
summary$IO_Max <- gsub("%", "", summary$IO_Max) #remove the percent sign
summary[, c(5:15)] <- sapply(summary[, c(5:15)], as.numeric)

## create genoAPA column
summary$genoAPA <- summary$genoAPA <- as.factor(paste(summary$genotype,summary$APA, sep="_"))
summary$genoAPA <- factor(summary$genoAPA, levels = c("WT_control", "FMR1-KO_control", "WT_trained", "FMR1-KO_trained"))
summary$genoAPA

### improve indiviual session dataframe ----
str(wtfmr1)
wtfmr1 <- wtfmr1[c(60:61,62,1:45)]  # removes summary columns

##create APA dataframe to add APA to wtfrm1
APA <-  summary[c(1:3)] 
wtfmr1 <- left_join(wtfmr1,APA) 
wtfmr1 <- wtfmr1[c(1:2,49,3:48)]  # removes summary columns
wtfmr1$APA <- factor(wtfmr1$APA, levels = c("control", "trained"))

## create one level column for genotype*APA
wtfmr1$genoAPA <- wtfmr1$genoAPA <- as.factor(paste(wtfmr1$genotype,wtfmr1$APA, sep="_"))
wtfmr1$genoAPA <- factor(wtfmr1$genoAPA, levels = c("WT_control", "FMR1-KO_control", "WT_trained", "FMR1-KO_trained"))
wtfmr1$genoAPA

### Beahvior ggplots!!! -----

### first, melt the beahvior data to make long for facet wrap with ggplot
wtfmr1_long <- melt(wtfmr1, id=c("ind","genotype", "APA", "session", "filename", "genoAPA"))
wtfmr1_long$value <- as.numeric(wtfmr1_long$value)
str(wtfmr1_long)

## create the color palette
FentonPalette <- c('black','grey50','red','darkorange')

## basic format to behavior of 4 groups by session stat smooth
ggplot(data=wtfmr1, aes(as.numeric(x=session), y=TimeToSecondEntrance, color=genoAPA)) + 
  stat_smooth() + theme_bw() + scale_colour_manual(values=FentonPalette) + 
  scale_x_continuous(name =NULL, 
                     labels=c("1" = "Pretraining", "2" = "Training 1", 
                              "3" = "Training 2", "4" = "Training 3", "5" = "Retention"))

## Plots of time in seconds: saved as 1_TimeMeasures
wtfmr1_long %>% 
  filter(grepl("Time", variable)) %>% 
  filter(!grepl("TotalTime|pTime", variable)) %>% 
  ggplot(aes(as.numeric(x=session), y=value, color=genoAPA)) +
  stat_smooth() + facet_wrap(~variable, scales = "free_y") +
  scale_y_continuous(name="Time (s)") + 
  scale_x_continuous(name =NULL, 
                     labels=c("1" = "Pretraining", "2" = "Training 1", 
                              "3" = "Training 2", "4" = "Training 3", "5" = "Retention")) +
  scale_colour_manual(values=FentonPalette, name="Treatment Group",
                      labels=c("WT Control", "FMR1-KO Control", "WT Trained", "FMR1-KO Trainined")) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), legend.position = c(0.85, 0.25))

## Plots of pTime (probability): saved as 1_TimeProbability
wtfmr1_long %>% filter(grepl("pTime", variable)) %>% 
  ggplot(aes(as.numeric(x=session), y=value, color=genoAPA)) +
  stat_smooth() + facet_wrap(~variable, scales = "free_y") +
  scale_y_continuous(name="Probability") + 
  scale_x_continuous(name =NULL, 
                     labels=c("1" = "Pretraining", "2" = "Training 1", 
                              "3" = "Training 2", "4" = "Training 3", "5" = "Retention")) +
  scale_colour_manual(values=FentonPalette, name="Treatment Group",
                      labels=c("WT Control", "FMR1-KO Control", "WT Trained", "FMR1-KO Trainined")) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

## Plots of path (distance): saved as 1_Path
wtfmr1_long %>% filter(grepl("Path", variable)) %>% 
  ggplot(aes(as.numeric(x=session), y=value, color=genoAPA)) +
  stat_smooth() + facet_wrap(~variable, scales = "free_y") +
  scale_y_continuous(name="Distance (m)") + 
  scale_x_continuous(name =NULL, 
                     labels=c("1" = "Pretraining", "2" = "Training 1", 
                              "3" = "Training 2", "4" = "Training 3", "5" = "Retention")) +
  scale_colour_manual(values=FentonPalette, name="Treatment Group",
                      labels=c("WT Control", "FMR1-KO Control", "WT Trained", "FMR1-KO Trainined")) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), legend.position = c(0.85, 0.25))

## Plots of speed: saved as 1_Speed
wtfmr1_long %>% 
  filter(grepl("Speed|speed", variable))%>% 
  ggplot(aes(as.numeric(x=session), y=value, color=genoAPA)) +
  stat_smooth() + facet_wrap(~variable, scales = "free_y") +
  scale_y_continuous(name="Speed (cm/s)") + 
  scale_x_continuous(name =NULL, 
                     labels=c("1" = "Pretraining", "2" = "Training 1", 
                              "3" = "Training 2", "4" = "Training 3", "5" = "Retention")) +
  scale_colour_manual(values=FentonPalette, name="Treatment Group",
                      labels=c("WT Control", "FMR1-KO Control", "WT Trained", "FMR1-KO Trainined")) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

## Other Plots of speed: Rayleight, Polar, Annular Stuff : Saved as 1-Other
wtfmr1_long %>% 
  filter(grepl("Rayleigh|Polar|Annular", variable)) %>% 
  ggplot(aes(as.numeric(x=session), y=value, color=genoAPA)) +
  stat_smooth() + facet_wrap(~variable, scales = "free_y") +
  scale_y_continuous(name=NULL) + 
  scale_x_continuous(name =NULL, 
                     labels=c("1" = "Pretraining", "2" = "Training 1", 
                              "3" = "Training 2", "4" = "Training 3", "5" = "Retention")) +
  scale_colour_manual(values=FentonPalette, name="Treatment Group",
                      labels=c("WT Control", "FMR1-KO Control", "WT Trained", "FMR1-KO Trainined")) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),legend.position = c(0.775, 0.125))

### Electrophysiology melting and graphing!!! ---- 

## first, melt the data
summary_long <- melt(summary, id=c("ind","genotype", "APA", "session", "filename", "genoAPA"))
summary_long$value <- as.numeric(summary_long$value)
str(summary_long)

## plot all the summary_long data!: saved as 1-SummaryValues
summary_long %>%
  ggplot(aes(x=genoAPA, y=value, color=genoAPA)) +
  geom_boxplot() +  facet_wrap(~variable, scales = "free_y") +
  theme_bw() +
  scale_y_continuous(name=NULL) + 
  scale_x_discrete(name=NULL) + 
  scale_colour_manual(values=FentonPalette, name="Treatment Group",
                      labels=c("WT Control", "FMR1-KO Control", "WT Trained", "FMR1-KO Trainined")) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.9, 0.15),
        axis.text.x = element_blank(), 
        axis.ticks = element_blank())

### Correlation plots!!! ----

## first, make the data a matrix with genoAPAsessionInd as the row names
wtfmr1_matrix <- wtfmr1    #create new dataframe
wtfmr1_matrix$genoAPAsessionInd <- as.factor(paste(wtfmr1_matrix$genoAPA, wtfmr1_matrix$session, wtfmr1_matrix$ind, sep="_")) #create genoAPAsessionInd column
rownames(wtfmr1_matrix) <- wtfmr1_matrix$genoAPAsessionInd     # set $genoAPAsessionInd as rownames
wtfmr1_matrix <- wtfmr1_matrix[-c(1:2,4:6,50:51)] #delete all non-numeric columns
head(wtfmr1_matrix)
str(wtfmr1_matrix)

## next, compute a correlation matrix and melt
wtfmrt_cormat <- round(cor(wtfmr1_matrix),2) # compute correlations
wtfmrt_cormatlong <- melt(wtfmrt_cormat) # melt
head(wtfmrt_cormatlong)

## heatmap NOT clustered!!!
ggplot(data = wtfmrt_cormatlong, aes(x=X1, y=X2, fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


summary_matrix <- summary
summary_matrix$genoAPAsessionInd <- as.factor(paste(summary_matrix$genoAPA, summary_matrix$session, summary_matrix$ind, sep="_")) #create genoAPAsessionInd column
rownames(summary_matrix) <- summary_matrix$genoAPAsessionInd     # set $genoAPAsessionInd as rownames
summary_matrix <- summary_matrix[-c(1:4,16:18)] #delete all non-numeric columns
head(summary_matrix)

### Not used but useful single plots ----

## simple box plot faceted by genotype colored by genotype
ggplot(data = wtfmr1, aes(x = session, y = TimeToFirstEntrance, fill=genoAPA)) +
  geom_boxplot() + facet_wrap( ~ genoAPA)
## simple box plot colored by genotype
ggplot(data = wtfmr1, aes(x = session, y = TimeToSecondEntrance, by=genoAPA, fill=genoAPA)) +
  geom_boxplot()
## line plots for every individual shown individually
ggplot(data=wtfmr1, aes(as.numeric(x=session), y=PathToSecondEntrance, by=ind, color=genoAPA)) +
  geom_line() + facet_wrap( ~ ind)
## line plots for every individual shown by group
ggplot(data=wtfmr1, aes(as.numeric(x=session), y=TimeToSecondEntrance, by=ind, color=genoAPA)) +
  geom_line() + facet_wrap( ~ genoAPA)
## line plots for every individual shown by group colored by individual
ggplot(data=wtfmr1, aes(as.numeric(x=session), y=PathToSecondEntrance, by=ind, color=ind)) +
  geom_line() + facet_wrap(~ genoAPA) 

## Time to 2nd entrance of 4 groups by session stat smooth
ggplot(data=wtfmr1, aes(as.numeric(x=session), y=TimeToSecondEntrance, color=genoAPA)) + 
  stat_smooth() + # plots line with confidence interval
  scale_y_continuous(name="Time to 2nd Entrance (s)") + # renames x axis
  scale_x_continuous(name =NULL, #removed y axis label
                     labels=c("1" = "Pretraining", "2" = "Training 1",  # renmaes x axis labeles
                              "3" = "Training 2", "4" = "Training 3", "5" = "Retention")) +
  scale_colour_manual(values=FentonPalette, name="Treatment Group", # says use color palete, rename legend
                      labels=c("WT Control", "FMR1-KO Control", "WT Trained", "FMR1-KO Trainined")) + 
  theme_bw() # removes grey background

## some statistics ----
lm1 = lm(data=wtfmr1, TimeToSecondEntrance~genoAPA+session+genoAPA:session)
summary(lm1)
xtabs(data=wtfmr1, ~genoAPA+session)
anova(lm1)

