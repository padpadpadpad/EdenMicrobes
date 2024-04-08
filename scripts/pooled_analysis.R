# ---------------------------
# Purpose of script: to look at differences in composition and diversity between quadrat corners, quadrats, beds, biomes, and ecosystems 
#
# What this script does:
# 1. combines all levels of replication below quadrat corner level
# 2. rarefies to look at differences in diversity
# 3. keeps non-rarefied data to estimate diversity
# 3. does PCoA to look for differences in composition
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-04-05
#
# Copyright (c) Daniel Padfield, 2024
#
# ---------------------------
#
# Notes: 
# 1. Need to set seed so rarefaction always gives the same answer
# 
# ---------------------------

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)){
  install.packages("librarian")
}
# load packages
librarian::shelf(adw96/breakaway, phyloseq, microbiome, tidyverse)

## ---------------------------

# set seed
set.seed(42)

# load in bacterial phyloseq object
ps_bact <- readRDS('data/sequencing/processed/phyloseq/ps_bact_sub.rds')
ps_bact

summarize_phyloseq(ps_bact)

# grab sample data out
d_samp <- sample_data(ps_bact) %>%
  data.frame()

# aggregate samples by sample
ps_bact2 <- merge_samples(ps_bact, 'sample')
summarize_phyloseq(ps_bact2)

ps_bact
ps_bact2

# look at sequencing numbers
ps_bact2 %>%
  phyloseq::sample_sums() %>%
  sort()
# ok there are 4 samples with under 15000 reads, so lets drop those and rarefy to 16000 reads

# remove samples
ps_bact2 <- prune_samples(sample_sums(ps_bact2) > 15000, ps_bact2)

# check frequency count table
freq_table <- build_frequency_count_tables(otu_table(ps_bact2))
freq_ct <- freq_table[[1]]
