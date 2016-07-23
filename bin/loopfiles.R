## uses read.delim function to read text files
## use "skip=43" for taqman
## use "skip=44" for sybr

taqman <- read.delim("2016-07-20_112033_RA.txt", skip=43, header = TRUE, sep = "\t" )[ ,1:24]
sybr <- read.delim("2016-07-20_115248_ShaynNicole_CREB_PkcZ_Htt.txt", skip=44, header = TRUE, sep = "\t" ) 

file_list <- list.files(pattern = ".txt") #creates a string will all the .xls in a diretory for us to loop through



##taqman exampple
rm(rawdata) # first removed any dataframe called rawdata, if any
rm(temp_dataset) # first removed any dataframe called rawdata, if any

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("rawdata")){
    rawdata <- read.delim(file, skip=43, header = TRUE, sep = "\t" )[ ,1:24]
  }
  
  # if the merged dataset does exist, append to it
  if (exists("rawdata")){
    temp_dataset <- read.delim(file, skip=43, header = TRUE, sep = "\t" )[ ,1:24]
    rawdata<-rbind(rawdata, temp_dataset)
    rm(temp_dataset)
  }
  
}

names(rawdata)
str(rawdata)



##sybr exampple
rm(rawdata) # first removed any dataframe called rawdata, if any
rm(temp_dataset) # first removed any dataframe called rawdata, if any

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("rawdata")){
    rawdata <- read.delim(file, skip=44, header = TRUE, sep = "\t" )[ ,1:24]
  }
  
  # if the merged dataset does exist, append to it
  if (exists("rawdata")){
    temp_dataset <- read.delim(file, skip=44, header = TRUE, sep = "\t" )[ ,1:24]
    rawdata<-rbind(rawdata, temp_dataset)
    rm(temp_dataset)
  }
  
}

names(rawdata)
str(rawdata)