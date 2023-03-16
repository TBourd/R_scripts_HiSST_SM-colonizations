#'The script is used to assign a sequence type to Serratia marcescens isolates or dominant eDNA 
#'by using the HiSST scheme. The HiSST scheme is a three-locus sequence-typing (ST) system that includes the bssA, gabR, and dhaM loci. 
#'
#'The script reads in a HiSST database that contains information about each HiSST ST and the corresponding bssA, gabR, and dhaM alleles. 
#'
#'The script then assigns a ST to an isolate or eDNA sample based on the alleles present at the bssA, gabR, and dhaM loci. 
#'There are two options for providing input data: (1) entering ST data manually, or (2) importing ST data from three separate files that contain the bssA, gabR, and dhaM ST data, respectively. 
#'
#'The script includes some additional data cleaning steps, such as removing samples with certain names or replacing incorrect sample names with correct ones.
#'
#'The following packages are used in the script: dplyr.


#' HiSST-Assignation

#' This script assigns a HiSST-profile Sequence Type (ST) to the corresponding S. marcescens isolate(s) or dominant eDNA, based on the Short Sequence Type (SST) of each locus (i.e. bssA, gabR and dhaM) from the HiSST scheme.

#' The file 'v2022-11_HiSST_database.txt' needed for this script is available at https://github.com/TBourd/R_scripts_HiSST_SM-outbreaks/tree/main/Step1_Assign-ST-infos.

library(dplyr)

setwd("~/Step2_Assign-ST_infos/") # set the working directory path of your project

# import HiSST database
HiSST <- read.table("v2022-11_HiSST_database.txt", header = TRUE)

# combine loci-ST profiles of database into one column
combination_ST <- mapply(paste, sep = "-", HiSST$bssA, HiSST$gabR, HiSST$dhaM)
HiSST.cb <- data.frame(HiSST, combination_ST)

# import data
# option 1: without input file
ID <- c("a", "b", "c")
query <- data.frame(ID = ID, bssA = c(2, 4, 3), gabR = c(1, 7, 1), dhaM = c(1, 4, 1))

# option 2: import .txt file
# these tables are obtained from the R script "Script_dada2_HiSST_v04-2022.R" (lines 246-273)
bssA <- read.table("bssA_Samples_and_ASV-ST.txt", header = TRUE)
dhaM <- read.table("dhaM_Samples_and_ASV-ST.txt", header = TRUE)
gabR <- read.table("gabR_Samples_and_ASV-ST.txt", header = TRUE)

bssA <- data.frame(ID = gsub("-bssA.*", "", bssA$Samples), bssA = gsub("ST", "", bssA$ST))
dhaM <- data.frame(ID = gsub("-dhaM.*", "", dhaM$Samples), dhaM = gsub("ST", "", dhaM$ST))
gabR <- data.frame(ID = gsub("-gabR.*", "", gabR$Samples), gabR = gsub("ST", "", gabR$ST))

query <- full_join(bssA, dhaM, by = "ID")
query <- full_join(query, gabR, by = "ID")

# Check for any missing values in the merged dataset
sum(is.na(query))

# Remove rows with missing values
query.noNA <- na.omit(query)

#__________________________________________________________________________________#
#### ____ Optional ____ Edit sample names ####

library(stringr)

# Edit sample names
query.noNA$ID <- gsub(pattern = "(1-BB|ST)", replacement = "BB", x = query.noNA$ID)

# Remove samples to be removed
sp.rm <- c("BB1b", "BB28-Sm1", "BB28-Sm2", "BWD155-1250-Sm4", "BWD162-1313-Sm3", "BB10a", "BB47", "BB48", "BB49", "BB50")
query.noNA.rm <- query.noNA[setdiff(1:nrow(query.noNA), grep(paste(sp.rm, collapse = "|"), query.noNA$ID)),]

# Find sample names that contain characters "00" and ("-1" or "-P")
wrong.names <- query.noNA.rm[str_detect(query.noNA.rm$Samples, "00") & str_detect(query.noNA.rm$Samples, "(-1|-P)"), 1]
wrong.names

# Standardize sample IDs
library(dplyr)

query.noNA.rm$ID <- query.noNA.rm$ID %>%
  sub("BD0021-172","BD172-0021",.) %>%
  sub("BD0035-PLM","BDPLM-0035",.) %>%
  sub("BD0021-173","BD173-0021",.) %>%
  sub("WD0021-PLM","WDPLM-0021",.) %>%
  sub("WD0035-173","WD173-0035",.) %>%
  sub("BD0035-173","BD173-0035",.) %>%
  sub("WD0035-172","WD172-0035",.) %>%
  sub("BD0021-PLM","BDPLM-0021",.) %>%
  sub("WD0035-PLM","WDPLM-0035",.) %>%
  sub("WD0021-173","WD173-0021",.) %>%
  sub("WD0021-172","WD172-0021",.) %>%
  sub("BD0035-172","BD172-0035",.) %>%
  sub("WD0028-180","WD180-0028",.) %>%
  sub("BD0028-180","BD180-0028",.) %>%
  sub("WD0014-180","WD180-0014",.) %>%
  sub("BD0014-180","BD180-0014",.) %>%
  sub("BD0014-155","BD155-0014",.) %>%
  sub("BD0028-152","BD152-0028",.) %>%
  sub("WD0014-152","WD152-0014",.) %>%
  sub("BD0014-152","BD152-0014",.) %>%
  sub("WD0028-152","WD152-0028",.) %>%
  
  sub("BB1a","BB1",.) %>%
  sub("BB10b","BB10",.) %>%
  sub("1PLM","-PLM",.) %>%
  sub("BD1","BD-",.) %>%
  sub("WD1","WD-",.) %>%
  sub("BD-b-2wD","BD1b-2wD",.) %>%
  sub("BDPLM","BD-PLM",.) %>%
  sub("WDPLM","WD-PLM",.) %>%
  
  sub("PLM","HWSs",.) %>%#Change PLM for HWS (Hand Washing System Intensive unit)
  sub("HWSsS","HWSs",.) %>%#Change HWSsS for HWSs (Hand Washing System Intensive unit)
  sub("HWSsR","HWSr",.) %>%#Change HWSsR for HWSr (Hand Washing System Intermediate unit)
  sub("CF","FK",.) %>%#Change CF for FK (Family Kitchen)
  sub("SA","LR",.)#Change SA for LR (Lactation Room)


#### _____ End Import and Edit data _____ ####


# Combine loci-ST profiles of strains in 1 column
combin_query <-setNames(as.data.frame(mapply(paste, sep = "-",query.noNA.rm$bssA,query.noNA.rm$gabR,query.noNA.rm$dhaM)),
                        nm="combination_ST")
row.names(combin_query) <-1:nrow(combin_query)
df_query <- data.frame(ID=query.noNA.rm$ID,combination_ST=combin_query)

# Assignation of HiSST ID for each query Strains with the data base
assign_HiSST_profile=merge(df_query,HiSST.cb, by=c("combination_ST"),all.x = TRUE)

Strain.ST <- data.frame(ID=assign_HiSST_profile$ID,HiSST=assign_HiSST_profile$HiSST,
                        HiSST.profile=assign_HiSST_profile$combination_ST)
#Strain.ST <- with(Strain.ST,  Strain.ST[order(ID) , ]) #Sort data frame by ID 

write.table(Strain.ST,file="HiSST ID of strains_10-2022.txt", sep="\t",row.names = FALSE, quote = F)


#### _____ Assign IDs for each eDNA HiSST profile _____ ####

NA.ST <- Strain.ST[rowSums(is.na(Strain.ST)) > 0, ]

df1 <- aggregate(NA.ST[1], NA.ST[3], FUN = function(X) paste(unique(X), collapse=", ")) # Combining duplicated rows

## Assign ST IDs for each sample ##
library(dplyr)

df2 <- data.frame(df1, HiSST = seq(1:nrow(df1))) # Create ID for each unique eDNA HiSST-profile
df.join <- left_join(NA.ST, df2, by = "HiSST.profile")

ST.env <- data.frame(ID = df.join$ID.x, HiSST = paste0("e",df.join$HiSST.y))


#### _____ Combine clinical strain and eDNA profiles into the file "Samples_info.txt" ___ ####

Samples_info <- read.table("Info_Samples.txt",header = TRUE, sep = "\t")

# Edit sample names
sp.rm <- c("BB47","BB48","BB49","BB50")
Samples_info <- Samples_info %>%
  filter(!grepl(paste(sp.rm,collapse="|"), Samples)) %>%
  mutate(Samples = sub("BDPLM","BD-PLM", Samples),
         Samples = sub("WDPLM","WD-PLM", Samples),
         Samples = sub("PLM","HWSs", Samples),
         Samples = sub("HWSsS","HWSs", Samples),
         Samples = sub("HWSsR","HWSr", Samples),
         Samples = sub("CF","FK", Samples),
         Samples = sub("SA","LR", Samples))


ALL.ST <- rbind(na.omit(Strain.ST[1:2]),ST.env)
ALL.ST.f <- data.frame(Samples = ALL.ST$ID, HiSST = ALL.ST$HiSST)

St_inf <- inner_join(Samples_info, ALL.ST.f, by = "Samples")

write.table(St_inf,file="Sp-info_ST_BB+eDNA_maxAbund_10-2022_no47-50.txt", sep="\t",row.names = FALSE, quote = F)


#--------- END - Optional: Additional step ---------

## Combine HiSST profile with ST of each locus ##
HiSST_profile <- full_join(ALL.ST,query.noNA.rm, by = "ID")
HiSST_profile <- data.frame(HiSST = HiSST_profile$HiSST, HiSST_profile[3:5])

write.table(HiSST_profile,file="HiSST-profile_BB+eDNA_maxAbund_10-2022_no47-50.txt", sep="\t",row.names = FALSE, quote = F)

