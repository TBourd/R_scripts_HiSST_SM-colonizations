#'HiSST-Assignation


#'This R script allows to identify the HiSST-profile ID of the corresponding S. marcecens isolate(s) or dominant eDNA,
#'by using the sequence type (ST) ID of each HiSST-locus (i.e. bssA, gabR and dhaM) of the isolate(s).
#'
#'The last section of this script allows to add a new HiSST profile to the database.
#'
#'The file 'HiSST_Database.txt' needed for this script is available 
#'at https://github.com/TBourd/R_scripts_for_HiSST_scheme/tree/main/Data.


  
setwd("~/Step2_Assign-ST_infos/") #Set the working directory path of your project
  
HiSST <- read.table("~/v1.4-2022_HiSST_database.txt",header = TRUE) #Import Database 


#'Combine loci-ST profiles of data base in 1 column  
combination_ST <- mapply(paste, sep = "-", HiSST$bssA, HiSST$gabR,HiSST$dhaM)
  
HiSST.cb <- data.frame(HiSST,combination_ST)

## -- Import data --

#### Without input file ####
#ID <- c("a","b","c") #'Enter the ID of the strain(s) : e.g. strain #a = "a"
#query <- data.frame(ID=ID,bssA=c(2,4,3),gabR=c(1,7,1),dhaM=c(1,4,1)) #'Enter ST of each HiSST-locus

#### Import .txt file ####
# These tables are obtained from the R script "dada2_for_HiSST_scheme.R" (lines 242-266).

bssA <- read.table("~/bssA_Samples_and_ASV-ST.txt",header = TRUE)
dhaM <- read.table("~/dhaM_Samples_and_ASV-ST.txt",header = TRUE)
gabR <- read.table("~/gabR_Samples_and_ASV-ST.txt",header = TRUE)

bssA <- data.frame(ID=bssA$Samples, bssA=bssA$ST)
bssA$ID <- gsub("-bssA.*", "", bssA$ID)
bssA$bssA <- gsub("ST", "", bssA$bssA)
dhaM <- data.frame(ID=dhaM$Samples, dhaM=dhaM$ST)
dhaM$ID <- gsub("-dhaM.*", "", dhaM$ID)
dhaM$dhaM <- gsub("ST", "", dhaM$dhaM)
gabR <- data.frame(ID=gabR$Samples, gabR=gabR$ST)
gabR$ID <- gsub("-gabR.*", "", gabR$ID)
gabR$gabR <- gsub("ST", "", gabR$gabR)

library(dplyr)
query <- full_join(bssA, dhaM, by = "ID")
query <- full_join(query, gabR, by = "ID")

query.noNA <- na.omit(query)


#__________________________________________________________________________________#
#### ____ Optional ____ Edit sample names ####

query.noNA$ID =sub("1-BB","BB",query.noNA$ID)
query.noNA$ID =sub("ST", "",query.noNA$ID)
sp.rm <- c("BB1b", "BB28-Sm1", "BB28-Sm2","BWD155-1250-Sm4","BWD162-1313-Sm3","BB10a","BB47","BB48","BB49","BB50") 
query.noNA.rm <- query.noNA[ grep(paste(sp.rm,collapse="|"), query.noNA$ID, invert = TRUE) , ] # Delete sample that contain '001-'


wrong.names <- query.noNA.rm[grepl("00", query.noNA.rm$Samples) 
                         & (grepl("-1", query.noNA.rm$Samples) | grepl("-P", query.noNA.rm$Samples)), 1] # Find sample names that contain characters "00" and ("-1" or "-P")

query.noNA.rm$ID =sub("BD0021-172","BD172-0021",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0035-PLM","BDPLM-0035",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0021-173","BD173-0021",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0021-PLM","WDPLM-0021",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0035-173","WD173-0035",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0035-173","BD173-0035",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0035-172","WD172-0035",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0021-PLM","BDPLM-0021",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0035-PLM","WDPLM-0035",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0021-173","WD173-0021",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0021-172","WD172-0021",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0035-172","BD172-0035",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0028-180","WD180-0028",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0028-180","BD180-0028",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0014-180","WD180-0014",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0014-180","BD180-0014",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0014-155","BD155-0014",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0028-152","BD152-0028",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0014-152","WD152-0014",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD0014-152","BD152-0014",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD0028-152","WD152-0028",query.noNA.rm$ID)

query.noNA.rm$ID =sub("BB1a","BB1",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BB10b","BB10",query.noNA.rm$ID)
query.noNA.rm$ID =sub("1PLM","-PLM",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD1","BD-",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WD1","WD-",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BD-b-2wD","BD1b-2wD",query.noNA.rm$ID)
query.noNA.rm$ID =sub("BDPLM","BD-PLM",query.noNA.rm$ID)
query.noNA.rm$ID =sub("WDPLM","WD-PLM",query.noNA.rm$ID)

query.noNA.rm$ID =sub("PLM","HWSs",query.noNA.rm$ID) #Change PLM for HWS (Hand Washing System Intensive unit)
query.noNA.rm$ID =sub("HWSsS","HWSs",query.noNA.rm$ID) #Change HWSsS for HWSs (Hand Washing System Intensive unit)
query.noNA.rm$ID =sub("HWSsR","HWSr",query.noNA.rm$ID) #Change HWSsR for HWSr (Hand Washing System Intermediate unit)
query.noNA.rm$ID =sub("CF","FK",query.noNA.rm$ID) #Change CF for FK (Family Kitchen)
query.noNA.rm$ID =sub("SA","LR",query.noNA.rm$ID) #Change SA for LR (Lactation Room)

## End


## -- End Import data ----


#'Combine loci-ST profiles of strains in 1 column
combin_query <-setNames(as.data.frame(mapply(paste, sep = "-",query.noNA.rm$bssA,query.noNA.rm$gabR,query.noNA.rm$dhaM)),
                             nm="combination_ST")
row.names(combin_query) <-1:nrow(combin_query)
df_query <- data.frame(ID=query.noNA.rm$ID,combination_ST=combin_query)

#' Assignation of HiSST ID for each query Strains with the data base
assign_HiSST_profile=merge(df_query,HiSST.cb, by=c("combination_ST"),all.x = TRUE)

Strain.ST <- data.frame(ID=assign_HiSST_profile$ID,HiSST=assign_HiSST_profile$HiSST,
                        HiSST.profile=assign_HiSST_profile$combination_ST)
#Strain.ST <- with(Strain.ST,  Strain.ST[order(ID) , ]) #Sort data frame by ID 

write.table(Strain.ST,file="HiSST ID of strains_10-2022.txt", sep="\t",row.names = FALSE, quote = F)

#_______________________________

#### Assign IDs for each eDNA HiSST profile ####

NA.ST <- Strain.ST[rowSums(is.na(Strain.ST)) > 0, ]

df1 <- aggregate(NA.ST[1], NA.ST[3], FUN = function(X) paste(unique(X), collapse=", ")) # Combining duplicated rows

## Assign ST IDs for each sample ##
library(dplyr)

df2 <- data.frame(df1, HiSST = seq(1:nrow(df1))) # Create ID for each unique eDNA HiSST-profile
df.join <- left_join(NA.ST, df2, by = "HiSST.profile")

ST.env <- data.frame(ID = df.join$ID.x, HiSST = paste0("e",df.join$HiSST.y))

## Combine clinical strain and eDNA profiles into the file "Samples_info.txt" ##
Samples_info <- read.table("~/Step2_Assign-ST+infos/Info_Samples.txt",header = TRUE, sep = "\t")

sp.rm <- c("BB47","BB48","BB49","BB50") 
Samples_info <- Samples_info[ grep(paste(sp.rm,collapse="|"), Samples_info$Samples, invert = TRUE) , ]
Samples_info$Samples =sub("BDPLM","BD-PLM",Samples_info$Samples)
Samples_info$Samples =sub("WDPLM","WD-PLM",Samples_info$Samples)
Samples_info$Samples =sub("PLM","HWSs",Samples_info$Samples) #Change PLM for HWS (Hand Washing System Intensive unit)
Samples_info$Samples =sub("HWSsS","HWSs",Samples_info$Samples) #Change HWSsS for HWSs (Hand Washing System Intensive unit)
Samples_info$Samples =sub("HWSsR","HWSr",Samples_info$Samples) #Change HWSsR for HWSr (Hand Washing System Intermediate unit)
Samples_info$Samples =sub("CF","FK",Samples_info$Samples) #Change CF for FK (Family Kitchen)
Samples_info$Samples =sub("SA","LR",Samples_info$Samples) #Change SA for LR (Lactation Room)

ALL.ST <- rbind(na.omit(Strain.ST[1:2]),ST.env)
ALL.ST.f <- data.frame(Samples = ALL.ST$ID, HiSST = ALL.ST$HiSST)

St_inf <- inner_join(Samples_info, ALL.ST.f, by = "Samples")

write.table(St_inf,file="Sp-info_ST_BB+eDNA_maxAbund_10-2022_no47-50.txt", sep="\t",row.names = FALSE, quote = F)


## Combine HiSST profile with ST of each locus ##
HiSST_profile <- full_join(ALL.ST,query.noNA.rm, by = "ID")
HiSST_profile <- data.frame(HiSST = HiSST_profile$HiSST, HiSST_profile[3:5])

write.table(HiSST_profile,file="HiSST-profile_BB+eDNA_maxAbund_10-2022_no47-50.txt", sep="\t",row.names = FALSE, quote = F)




#--------- END ST-Assignation ---------
