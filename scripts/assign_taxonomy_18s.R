# assign taxonomy to Eden project sequences
# 18s

# uses the dada2 and phyloseq pipeline presented here:
# https://f1000research.com/articles/5-1492/v2 
# Bioconductor workflow for Microbiome data analysis: from raw reads to community analysis

# set seed ####
set.seed(42)

# load in packages
librarian::shelf(phyloseq, dada2, DECIPHER, tidyverse)

# download and load in reference databases
if(!file.exists('eden/ref_db/silva_132.18s.99_rep_set.dada2.fa.gz')){download.file('https://zenodo.org/records/1447330/files/silva_132.18s.99_rep_set.dada2.fa.gz?download=1', 'eden/ref_db/silva_132.18s.99_rep_set.dada2.fa.gz')}

silva_db <- 'eden/ref_db/silva_132.18s.99_rep_set.dada2.fa.gz'

# load in sequences
seqs_18s <- readRDS('eden/data/eden_euk_unique_sequence.rds')

# assign taxonomy in batches of 50,000
to_split <- seq(1, nrow(seqs_18s), by = 50000)
to_split2 <- c(to_split[2:length(to_split)]-1, nrow(seqs_18s))

taxtab = NULL

for(i in 1:length(to_split)){
  seqs_18s2 <- seqs_18s[to_split[i]:to_split2[i], ]
  taxtab2 <- assignTaxonomy(seqs_18s2, refFasta = silva_db, multithread = 4, tryRC = TRUE)
  if(!is.null(taxtab)){taxtab <- rbind(taxtab, taxtab2)}
  if(is.null(taxtab)){taxtab <- taxtab2}
}
# save files in case phylogeny does not run
saveRDS(taxtab, 'eden/data/taxtab_euk.rds')