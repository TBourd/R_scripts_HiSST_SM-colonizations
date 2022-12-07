library(dplyr)

setwd("~/Step1_Create_dendrogram_and_binary-matrix/")

#### _____ Step 1: import data and compile bssA, gabR and dhaM ____ ####

ST <- read.delim("bssA_ASV_samples.txt") #Change for bssA, then for gabR and for dhaM
#ST <- read.delim("gabR_ASV_samples.txt") 
#ST <- read.delim("dhaM_ASV_samples.txt") 

# Remove unwanted characters within sample names
char.rm <- c("-bssA","-bbsA","-dhaM","-gabR")
ST$Samples =sub(paste(char.rm,collapse="|"),"",ST$Samples)
ST$Samples =sub("1-BB","BB",ST$Samples)
#ST <- ST[-118,] #Remove duplicate data : [ which(duplicated(ST$Samples)) ]

#__ Remove non-expected ASV in clinical isolates = Keep only reads in dominance abundance and descard contamination ASV

# Find strains with more than 2 ASVs per locus
#k <- data.frame(Samples = NA, liste = NA)
#for (i in 1:nrow(BB)) {
#  mylist <- list(which(BB[i,2:ncol(BB)]>0))
#  j <- data.frame(Samples = BB[i,1],liste = sapply(mylist, paste, collapse=", "))
#  k <- rbind(k,j)
#}

BB <- ST[ grep("BB", ST$Samples, invert = FALSE) , ]
BB_max <- BB
for (y in 1:nrow(BB)) {
  for (i in 2:ncol(BB)) {
    ifelse(BB[y,i] < max(BB[y,2:ncol(BB)]) 
           , BB_max[y,i]<-0 , BB_max[y,i] <- BB[y,i])
  }
}
ST <- ST[ grep("BB", ST$Samples, invert = TRUE) , ]
ST <- rbind(ST,BB_max)

### __ Option 1: Remove reads with abundance less than 1% of the total sample abundance, or with less than 8 reads. 
ST_0.01 <- ST
for (y in 1:nrow(ST)) {
  for (i in 2:ncol(ST)) {
    ifelse(ST[y,i] < 0.01*sum(ST[y,2:ncol(ST)]) #| ST[y,i] < 8
           , ST_0.01[y,i]<-0 , ST_0.01[y,i] <- ST[y,i])
  }
}

### __ Option 2: Keep only reads in dominance abundance. (Useful to generate minimum spanning tree by MLST method)

# ! Warning : if a sample has no abundance for a locus, the values will be "1" for all ASV of the locus => Change the value by "0"
#' To verify if any sample has none abundance, 
# run : > which(rowSums(ST[2:ncol(ST)])==0)
# Then : > ST_max[which(rowSums(ST[2:ncol(ST)])==0), 2:ncol(ST_max)] <- 0

ST_max <- ST

for (y in 1:nrow(ST)) {
  for (i in 2:ncol(ST)) {
    ifelse(ST[y,i] < max(ST[y,2:ncol(ST)]) #| ST[y,i] < 10
           , ST_max[y,i]<-0 , ST_max[y,i] <- 1)
  }
}


#___________


## /!\ Change bssA or dhaM or gabR /!\
bssA<- ST_0.01 #or 'ST_max'if option 2 were choosen

bssA_gabR <- full_join(bssA,gabR, by = "Samples")
All_ST <- full_join(bssA_gabR,dhaM, by = "Samples")

All_ST[is.na(All_ST)] <- 0

#All_ST <- na.omit(All_ST) # Remove samples with missing locus (=with NA values)


#' Remove unwanted samples. 
# e.g. remove samples that are not part of the study :
sp.rm <- c("BB1b", "BB28-Sm1", "BB28-Sm2","BWD155-1250-Sm4","BWD162-1313-Sm3","BB10a", "BB47","BB48","BB49","BB50") 
All_ST.rm <- All_ST[ grep(paste(sp.rm,collapse="|"), All_ST$Samples, invert = TRUE) , ] # Delete sample that contain '001-'
#All_ST.rm$Samples =sub("BB10a","BB10",All_ST.rm$Samples)


#__________________________________________________________________________________#


#### _____ Step 2: Optional steps = Edit samples names and group samples ____ ####

####  Option 1 : Edit sample names ####
wrong.names <- All_ST.rm[grepl("00", All_ST.rm$Samples) 
                         & (grepl("-1", All_ST.rm$Samples) | grepl("-P", All_ST.rm$Samples)), 1] # Find sample names that contain characters "00" and ("-1" or "-P")

All_ST.rm$Samples =sub("BD0021-172","BD172-0021",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0035-PLM","BDPLM-0035",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0021-173","BD173-0021",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0021-PLM","WDPLM-0021",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0035-173","WD173-0035",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0035-173","BD173-0035",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0035-172","WD172-0035",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0021-PLM","BDPLM-0021",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0035-PLM","WDPLM-0035",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0021-173","WD173-0021",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0021-172","WD172-0021",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0035-172","BD172-0035",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0028-180","WD180-0028",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0028-180","BD180-0028",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0014-180","WD180-0014",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0014-180","BD180-0014",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0014-155","BD155-0014",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0028-152","BD152-0028",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0014-152","WD152-0014",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD0014-152","BD152-0014",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD0028-152","WD152-0028",All_ST.rm$Samples)

All_ST.rm$Samples =sub("BB1a","BB1",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BB10b","BB10",All_ST.rm$Samples)
All_ST.rm$Samples =sub("1PLM","-PLM",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD1","BD-",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WD1","WD-",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BD-b-2wD","BD1b-2wD",All_ST.rm$Samples)
All_ST.rm$Samples =sub("BDPLM","BD-PLM",All_ST.rm$Samples)
All_ST.rm$Samples =sub("WDPLM","WD-PLM",All_ST.rm$Samples)

All_ST.rm$Samples =sub("PLM","HWSs",All_ST.rm$Samples) #Change PLM for HWS (Hand Washing System Intensive unit)
All_ST.rm$Samples =sub("HWSsS","HWSs",All_ST.rm$Samples) #Change HWSsS for HWSs (Hand Washing System Intensive unit)
All_ST.rm$Samples =sub("HWSsR","HWSr",All_ST.rm$Samples) #Change HWSsR for HWSr (Hand Washing System Intermediate unit)
All_ST.rm$Samples =sub("CF","FK",All_ST.rm$Samples) #Change CF for FK (Family Kitchen)
All_ST.rm$Samples =sub("SA","LR",All_ST.rm$Samples) #Change SA for LR (Lactation Room)

## End edit sample names

#____________________________________________________________________________________________________________#


####  Option 2 : Create colors for each sample group in the dendrogram --- ####

grp.All <- All_ST.rm

### A. Clinical isolates
BB.grp <- data.frame(Samples = All_ST.rm[ grep("BB", All_ST.rm$Samples, invert = FALSE) , 1],group=1)
### B. 2019 sampling
sp.rm <- c("BD1b-2wD","9288","9295")
df.2019.grp <- data.frame(Samples = All_ST.rm[ grep(paste(sp.rm,collapse="|"), All_ST.rm$Samples, invert = FALSE) , 1],group=2)
### C. 2020 sampling
sp.rm <- c("BWD","BB","BD1b-2wD","9288","9295")
df.2020.grp <- data.frame(Samples = All_ST.rm[ grep(paste(sp.rm,collapse="|"), All_ST.rm$Samples, invert = TRUE) , 1],group=3)
### D. 2021 sampling
df.2021.grp <- data.frame(Samples = All_ST.rm[ grep("BWD", All_ST.rm$Samples, invert = FALSE) , 1],group=4)


#____________________________________________________________________________________________________________#


#### _____ Step 3: Compute matrix and create dendrogram____ ####

### remove site with no data ###
All_ST. <- data.frame(grp.All[,2:ncol(grp.All)], row.names = grp.All$Samples)
spe.marg = addmargins(as.matrix(All_ST.))
empty.sites = which(spe.marg["Sum",]==0)
spe <-  All_ST.[,-empty.sites]

#spe.write <- data.frame(Samples = rownames(spe),spe) # For Minimum Spanning tree analysis
#write.table(spe.write, "Sm_bin_ALL2019-2021_MaxAbund_10-2022.txt", sep = "\t", quote = FALSE, row.names = FALSE) # For Minimum Spanning tree analysis

#  Compute a binary matrix for Jaccard distance #
spe.bin <- spe[]
spe.bin[spe.bin> 0] <- 1 
dist.spe = proxy::dist(spe.bin, by_rows = TRUE, method = "Jaccard") # Jaccard distance for binary data

# Or Normalize data to compute Bray-Curtis distance
#spe.hel <- decostand(spe,method = "hellinger") #normalize data
#dist = proxy::dist(spe.hel, by_rows = TRUE, method = "Bray") #Bray-Curtis distance


#### Dendro ####
library(dendextend)

d2 <- hclust(dist.spe ,method = "average")
d1 = as.dendrogram(d2)

## Color labels by groups

# If Step 2 - option 2 was applied, run the following lines:
grp.col <- rbind(BB.grp, df.2019.grp, df.2020.grp, df.2021.grp)
colors_to_use <- data.frame(group=grp.col$group,row.names = grp.col$Samples)

# Apply group colors
lab.d1 <- labels(d1)
colors_to_use.ord <- colors_to_use[lab.d1,,drop=FALSE] #Reorder the row names according to the label order of d1

colors_to_use.ord$group = gsub(1,"orange",colors_to_use.ord$group)
colors_to_use.ord$group = gsub(2,"dark grey",colors_to_use.ord$group)
colors_to_use.ord$group = gsub(3,"navy",colors_to_use.ord$group)
colors_to_use.ord$group = gsub(4,"dark green",colors_to_use.ord$group)
colors_to_use.ord$group = gsub(5,"red4",colors_to_use.ord$group)

labels_colors(d1) <- c(colors_to_use.ord$group,h=0.8) #Color a dendrogram labels according to defined group

## Plot
par(mar=c(0,0,0,0),mgp=c(1,1,0),cex.axis=0.5)
d1 = d1 %>% 
  set("branches_lty", 1) %>%
  set("branches_k_col", h=0.8) %>%
  set("leaves_pch", 20) %>%
  set("leaves_col", "dark cyan") %>%
  set("labels_cex",0.75)


library(circlize)

circos.par(circle.margin=0.2)
circlize_dendrogram(d1)
