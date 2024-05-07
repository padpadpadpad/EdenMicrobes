# assign taxonomy to Eden project sequences
# 16s

# uses the dada2 and phyloseq pipeline presented here:
# https://f1000research.com/articles/5-1492/v2 
# Bioconductor workflow for Microbiome data analysis: from raw reads to community analysis

# set seed ####
set.seed(42)

# load in packages
librarian::shelf(phyloseq, dada2, DECIPHER, tidyverse)

# download and load in reference databases
if(!file.exists('eden/ref_db/silva_nr99_v138.1_train_set.fa.gz')){download.file('https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1', 'eden/ref_db/silva_nr99_v138.1_train_set.fa.gz')}
if(!file.exists('eden/ref_db/silva_species_assignment_v138.1.fa.gz')){download.file('https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1', 'eden/ref_db/silva_species_assignment_v138.1.fa.gz')}

silva_db <- 'eden/ref_db/silva_nr99_v138.1_train_set.fa.gz'
silva_spp_db <- 'eden/ref_db/silva_species_assignment_v138.1.fa.gz'

# load in sequences
seqs_16s <- readRDS('eden/data/eden_bact_unique_sequence.rds')

# assign taxonomy in batches of 50,000
to_split <- seq(1, nrow(seqs_16s), by = 50000)
to_split2 <- c(to_split[2:length(to_split)]-1, nrow(seqs_16s))

taxtab = NULL

for(i in 1:length(to_split2)){
  seqs_16s2 <- seqs_16s[to_split[i]:to_split2[i], ]
  taxtab2 <- assignTaxonomy(seqs_16s2, refFasta = silva_db, multithread = 4)
  #if(silva_spp_db != 'NULL'){taxtab2 <- addSpecies(taxtab2, refFasta = silva_spp_db, verbose = TRUE)}
  if(!is.null(taxtab)){taxtab <- rbind(taxtab, taxtab2)}
  if(is.null(taxtab)){taxtab <- taxtab2}
}

# save files in case phylogeny does not run
saveRDS(taxtab, 'eden/data/taxtab_bact.rds')