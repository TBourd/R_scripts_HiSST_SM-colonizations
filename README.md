# R_scripts_HiSST_SM-outbreaks
 Rscripts and data used for _S. marcescens_ outbreaks survey, using the HiSST scheme.
 
- The first R script in the file entitled "Step1_Create_dendrogram_and_binary-matrix" was used to create a circular dendrogram of the eDNA and the isolated samples from the survey (REF). Lines 164-165 of this script generate a binary matrix used for the Minimum Spanning Tree.

- The second R script in the file entitled "Step2_Assign-ST+infos" allows to identify the HiSST-profile Sequence Type (ST) of the corresponding _S. marcecens_ isolate(s) or dominant eDNA, by using the Short Sequence Type (SST) of each locus (i.e. _bssA_, _gabR_ and _dhaM_) from the HiSST scheme.
 
- The third file "Step3_GeoBurst" gathers the two files used to generate a Minimum Spanning Tree using the PHYLOViZ platform.
 
 For more information on the HiSST method, please see: 
Bourdin T, Monnier A, Benoit MÈ, Bédard E, Prévost M, Quach C, Déziel E, Constant P. A High-Throughput Short Sequence Typing Scheme for Serratia marcescens Pure Culture and Environmental DNA. Appl Environ Microbiol. 2021. doi: 10.1128/AEM.01399-21.
