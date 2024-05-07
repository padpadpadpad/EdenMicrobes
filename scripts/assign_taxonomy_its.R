# assign taxonomy to Eden project sequences
# ITS sequencing

# uses the dada2 and phyloseq pipeline presented here:
# https://f1000research.com/articles/5-1492/v2 
# Bioconductor workflow for Microbiome data analysis: from raw reads to community analysis

# set seed ####
set.seed(42)

# load in packages
librarian::shelf(phyloseq, dada2, DECIPHER, tidyverse)

# download and load in reference databases
# downloaded from here: https://doi.plutof.ut.ee/doi/10.15156/BIO/2938069
ref_db <- 'eden/ref_db/sh_general_release_dynamic_all_18.07.2023.fasta'

# load in sequences
seqs_its <- readRDS('eden/data/eden_fungi_unique_sequence.rds')

# assign taxonomy in batches of 50,000
to_split <- seq(1, nrow(seqs_its), by = 50000)
to_split2 <- c(to_split[2:length(to_split)]-1, nrow(seqs_its))

taxtab = NULL

for(i in 1:length(to_split)){
  seqs_its2 <- seqs_its[to_split[i]:to_split2[i], ]
  taxtab2 <- assignTaxonomy(seqs_its2, refFasta = ref_db, multithread = 4)
  if(!is.null(taxtab)){taxtab <- rbind(taxtab, taxtab2)}
  if(is.null(taxtab)){taxtab <- taxtab2}
}
# save files in case phylogeny does not run
saveRDS(taxtab, 'eden/data/taxtab_fungi.rds')