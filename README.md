# R_scripts_HiSST_SM-colonizations
 Rscripts and data used for _S. marcescens_ colonizations survey, using the HiSST scheme.
 
## - R script "Script_dada2_HiSST_v04-2022.R":
Originate from dada2 processing (https://benjjneb.github.io/dada2/index.html) and was adapted for HiSST analysis. At the end of this script, two text files are generated for further analysis described below.

## - File "Step1_Assign-ST-infos":
R script in this file is used to assign a sequence type to Serratia marcescens isolates or dominant eDNA by using the HiSST scheme. The script then assigns a ST to an isolate or eDNA sample based on the alleles present at the bssA, gabR, and dhaM loci. There are two options for providing input data: (1) entering ST data manually, or (2) importing ST data from three separate files that contain the bssA, gabR, and dhaM ST data, respectively. The script includes some additional data cleaning steps, such as removing samples with certain names or replacing incorrect sample names with correct ones. It creates the file "Sp-info_ST_BB+eDNA_maxAbund_10-2022_no47-50.txt" used in step 2 described below. The following packages are used in the script: dplyr.

## -  File "Step2_geoBURST_SNP_matrix-fasta":
R script in this file merges DNA sequences from three loci (bssA, gabR, and dhaM) associated with each ST in a dataset of aligned sequences in FASTA format. Fasta files used in this script are available on https://github.com/TBourd/R_scripts_for_HiSST_scheme/tree/main/Data
It creates a single sequence per ST by concatenating the three sequences and calculates a distance matrix based on single-nucleotide polymorphisms between each sequence pair. The resulting merged sequences was used to the dataset type "Aligned Sequences (FASTA)" employed in Minimum Spanning tree using the PHYLOViZ platform (see "Fig.2" in Bourdin et al., Appl Environ Microbiol. 2023). 
The script also includes some additional steps to generate a dendrogram and for modifying sequence IDs and exporting files. The 'Biostrings', 'dplyr' and 'dendextend' packages are used.
 
 _______________________________________________________
 ###### For more information on the HiSST method, please see:
 
- Bourdin T, Monnier A, Benoit MÈ, Bédard E, Prévost M, Quach C, Déziel E, Constant P. A High-Throughput Short Sequence Typing Scheme for Serratia marcescens Pure Culture and Environmental DNA. Appl Environ Microbiol. 2021. DOI: https://doi.org/10.1128/AEM.01399-21

- Bourdin T, Benoit MÈ, Monnier A, Bédard E, Prévost M, Charron D, Audy N, Gravel S, Sicard M, Quach C, Déziel E, Constant P. Serratia marcescens colonization in a neonatal intensive care unit has multiple sources, highlighting sink drains as a major reservoir
