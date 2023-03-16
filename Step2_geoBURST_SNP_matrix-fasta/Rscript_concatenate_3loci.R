#'This R script merges DNA sequences from three loci (bssA, gabR, and dhaM) associated with each ST in a dataset of aligned sequences in FASTA format. 
#'
#'Fasta files used in this script are available on https://github.com/TBourd/R_scripts_for_HiSST_scheme/tree/main/Data
#'
#'It creates a single sequence per ST by concatenating the three sequences and calculates a distance matrix based on single-nucleotide polymorphisms between each sequence pair. 
#'The resulting merged sequences was used to the dataset type "Aligned Sequences (FASTA)" employed in Minimum Spanning tree.
#'The resulting distance matrix is used to generate a dendrogram. 
#'The script also includes some additional code for modifying sequence IDs and exporting files. 
#'
#''The 'Biostrings', 'dplyr' and 'dendextend' packages are used.

setwd("~/Rscript_geoBURST_SNP_matrix-fasta/")

# Load required packages
library(Biostrings)
library(dplyr)

# Define a function that modifies sequence IDs in a fasta file
modify_seq_ids <- function(fasta_path) {
  # Read in the fasta file
  fasta_seqs <- read.fasta(fasta_path)
  
  # Modify the sequence IDs
  for (i in seq_along(fasta_seqs)) {
    attr(fasta_seqs[[i]],"name") <- gsub(".*ST", "SST", attr(fasta_seqs[[i]],"Annot"))
  }
  
  # Return the modified fasta sequences
  return(fasta_seqs)
}

# Read in the fasta files as lists of sequences
fasta_seqs_A <- readDNAStringSet("ST_bssA.fasta")
fasta_seqs_B <- readDNAStringSet("ST_gabR.fasta")
fasta_seqs_C <- readDNAStringSet("ST_dhaM.fasta")


# Read in the ST-profiles file as a data frame
st_profiles <- read.delim("v2023-03_eDNA_HiSST_database.txt", sep = "\t")
st_profiles <- st_profiles %>%
  mutate(across(everything(), as.character))


# Create an empty data frame to store the combined sequences
combined_seqs <- data.frame(ST_profile = character(),
                            fasta_file_A = character(),
                            fasta_file_B = character(),
                            fasta_file_C = character(),
                            stringsAsFactors = FALSE)

# Loop over each row of the ST-profiles data frame
for (i in seq_len(nrow(st_profiles))) {
  # Extract the SST numbers for this ST profile
  sst_A <- paste0("ST", st_profiles$bssA[i])
  sst_B <- paste0("ST", st_profiles$gabR[i])
  sst_C <- paste0("ST", st_profiles$dhaM[i])
  
  # Extract the DNA sequences for these SSTs from the fasta files
  dna_seq_A <- as.character(fasta_seqs_A[sst_A])
  dna_seq_B <- as.character(fasta_seqs_B[sst_B])
  dna_seq_C <- as.character(fasta_seqs_C[sst_C])
  
  # Combine the DNA sequences into a single string and add them to the data frame
  combined_seqs <- combined_seqs %>%
    add_row(ST_profile = st_profiles$HiSST[i],
            fasta_file_A = dna_seq_A,
            fasta_file_B = dna_seq_B,
            fasta_file_C = dna_seq_C)
}

# Print the resulting data frame
combined_seqs


# concatenate the fasta files
merged_fasta <- data.frame(HiSST = combined_seqs$ST_profile, 
                           seq = mapply(paste, sep = "", 
                                        combined_seqs$fasta_file_A, 
                                        combined_seqs$fasta_file_B,
                                        combined_seqs$fasta_file_C))


#### Export merged_fasta in fasta files ####

# Remove unwanted STs
#library(readr)

sp_info <- read.delim("Sp-info_ST_BB+eDNA_maxAbund_10-2022_no47-50.txt")
# keep only the rows in snp_matrix where the values in column "HiSST" appear in column "HiSST" of sp_info
merged_fasta_STrm <- merged_fasta[merged_fasta$HiSST %in% sp_info$HiSST, ]

# create a character vector in FASTA format
fasta_lines <- paste0(">", merged_fasta_STrm$HiSST, "\n", merged_fasta_STrm$seq, "\n")
# write the lines to a file
writeLines(fasta_lines, "merged_fasta.fasta")
write.table(merged_fasta_STrm, "merged_fasta_STrm.txt", quote = F, row.names = F, sep = "\t")






#### ______ Optional steps ______ ####



#### Compute a distance matrix based on SNPs ####

# Define a function to calculate SNPs
calculate_snps <- function(sequence1, sequence2) {
  # Convert sequences to character vectors
  sequence1 <- unlist(strsplit(sequence1, ""))
  sequence2 <- unlist(strsplit(sequence2, ""))
  
  # Calculate SNPs
  snps <- sum(sequence1 != sequence2)
  
  return(snps)
}

# Calculate pairwise SNPs
n_sequences <- nrow(merged_fasta)
distance_matrix <- matrix(0, nrow = n_sequences, ncol = n_sequences)
for (i in 1:(n_sequences-1)) {
  for (j in (i+1):n_sequences) {
    distance_matrix[i, j] <- calculate_snps(merged_fasta[[2]][i], merged_fasta[[2]][j])
    distance_matrix[j, i] <- distance_matrix[i, j]
  }
}

distance_df <- as.data.frame(distance_matrix, row.names = merged_fasta$HiSST)
colnames(distance_df) <- merged_fasta$HiSST

write.table(distance_df, "distance_matrix.txt", quote = F, row.names = F, sep = "\t")



##### Create a dendrogram with dendextend ####
library(dendextend)

dend <- as.dendrogram(hclust(as.dist(distance_df)))

par(mar=c(2,2,2,3),mgp=c(1,1,0),cex.axis=0.5)

dend = dend %>% 
  set("labels_col", h=0.8, k = 4) %>%
  set("branches_lwd", 1.5) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color", h=0.8, k = 4) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", "dark cyan") %>%
  set("labels_cex",0.75) 

plot(dend, main = "Phylogenetic Tree", cex = 0.8, edge.width = 1.5, edge.color = "gray")



#### Create a Neighbor-Joining tree ####

# load the necessary packages
library(dendextend)
library(ape)

seqs <- readDNAStringSet("merged_fasta.fasta")

# convert the DNAStringSet object to a DNAbin object
dnabin_sequences <- as.DNAbin(seqs)

# calculate the pairwise distances between sequences
dist_matrix <- dist.dna(dnabin_sequences, model = "K80")

# construct the phylogenetic tree using the distance matrix
tree <- nj(dist_matrix)

# plot the dendrogram with bootstrap values
plot(tree, main = "Phylogenetic Tree", type = "fan", tip.color = "black", cex = 0.8, edge.width = 1.5, edge.color = "gray")
