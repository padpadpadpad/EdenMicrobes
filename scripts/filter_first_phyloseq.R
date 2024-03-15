#---------------------------------------------------------------#
# first filter of the phyloseq objects and quick and dirty plot #
#---------------------------------------------------------------#

# for each set of sequencing (bacteria, fungi, eukaryotes):
# read in sample info
# create new sample data
# remove samples with less than 1000 reads
# remove ASVs that are not assigned to a phylum
# remove ASVs with read length < 100

# there is one sample that is there twice
# B802B3r1a and B802B3r1
# Has been sequenced twice as much as all other individual PCRs, B802B3r1 was present in two wells in the PCR libraries
# See how it clusters out

# load in packages
librarian::shelf(phyloseq, tidyverse)

#--------------#
# for fungi ####
#--------------#

# load in fungi phyloseq object
ps_fungi <- readRDS('data/sequencing/processed/phyloseq/ps_fungi.rds')

# load in new sample data
sample_data <- read.csv('data/Eden_PCRs_sample_info.csv') %>%
  filter(!is.na(Ecosystem))

new_sample_data <- sample_data(ps_fungi) %>%
  data.frame()

filter(sample_data,! sample_id %in% new_sample_data$sample_id)

new_sample_data <- sample_data(ps_fungi) %>%
  data.frame() %>%
  rownames_to_column('sample_info') %>%
  # add column for pcr replicate by grabbing the last four characters of sample info
  mutate(pcr = str_sub(sample_info, -1, -1),
         # add column for extraction replicate by grabbing the number between the r and p in sample info
         extraction = str_extract(sample_info, '(?<=r)\\d+(?=p)'),
         extraction = ifelse(is.na(extraction), 1, extraction),
         # extract all text before the r
         sample = str_extract(sample_info, '.*(?=pcr)'),
         sample = str_extract(sample, '.*(?=r)'),
         # add column for quadrat corner
         quadrat_corner = str_sub(sample, -1, -1),
         quadrat = str_sub(sample, -2, -2),
         bed = str_sub(sample, 1, -3)) %>%
  select(-ecosystem) %>%
  left_join(., select(sample_data, sample_id, ecosystem = Ecosystem))
# sample id B802B3r1a needs checking, not sure why a is here

# unique biome assignment
biomes <- tibble(ecosystem = unique(new_sample_data$ecosystem)) %>%
  mutate(biome = ifelse(ecosystem %in% c('South_Africa', 'Australia', 'Vines', 'Citrus', 'Med'), 'Mediterranean', 'Rainforest'))

new_sample_data <- left_join(new_sample_data, biomes) %>%
  select(sample_info, sample_id, biome, ecosystem, bed, quadrat, sample, quadrat_corner, extraction, pcr) %>%
  column_to_rownames('sample_info') 

sample_data(ps_fungi) <- sample_data(new_sample_data)

sample_sums(ps_fungi) %>%
  log10() %>%
  hist()

sample_sums(ps_fungi) %>%
  sort()

# filter out samples with less than 1000 reads
ps_fungi <- prune_samples(sample_sums(ps_fungi) > 1000, ps_fungi)

sample_names(ps_fungi) %>%
  sort()

sample_sums(ps_fungi) %>%
  log10() %>%
  hist()

# remove ASVs where Phylum = NA
ps_fungi_sub <- prune_taxa(taxa_names(ps_fungi)[tax_table(ps_fungi)[, 'Phylum'] != 'NA'], ps_fungi)

# look at distribution of read length
rownames(tax_table(ps_fungi_sub)) %>%
  nchar() %>%
  hist()
rownames(tax_table(ps_fungi_sub)) %>% nchar() %>% summary()

# remove ASVs with read length < 100
ps_fungi_sub <- prune_taxa(taxa_names(ps_fungi_sub)[nchar(rownames(tax_table(ps_fungi_sub))) > 100], ps_fungi_sub)

# saved 
saveRDS(ps_fungi_sub, 'data/sequencing/processed/phyloseq/ps_fungi_sub.rds')

#-------------------------------------#
# load in bacteria phyloseq object ####
#-------------------------------------#

# load in bact data
ps_bact <- readRDS('data/sequencing/processed/phyloseq/ps_bact.rds')

# load in new sample data
sample_data <- read.csv('data/Eden_PCRs_sample_info.csv') %>%
  filter(!is.na(Ecosystem))

# make new sample data
new_sample_data <- sample_data(ps_bact) %>%
  data.frame()

# any samples in our sample data not in the new sample data?
filter(sample_data,! sample_id %in% new_sample_data$sample_id)
# 8 samples are not in our bacteria samples that were in the full sample info

# create a new sample data object
new_sample_data <- sample_data(ps_bact) %>%
  data.frame() %>%
  rownames_to_column('sample_info') %>%
  # add column for pcr replicate by grabbing the last four characters of sample info
  mutate(pcr = str_sub(sample_info, -1, -1),
         # add column for extraction replicate by grabbing the number between the r and p in sample info
         extraction = str_extract(sample_info, '(?<=r)\\d+(?=p)'),
         extraction = ifelse(is.na(extraction), 1, extraction),
         # extract all text before the r
         sample = str_extract(sample_info, '.*(?=pcr)'),
         sample = str_extract(sample, '.*(?=r)'),
         # add column for quadrat corner
         quadrat_corner = str_sub(sample, -1, -1),
         quadrat = str_sub(sample, -2, -2),
         bed = str_sub(sample, 1, -3)) %>%
  select(-ecosystem) %>%
  left_join(., select(sample_data, sample_id, ecosystem = Ecosystem))

# unique biome assignment
biomes <- tibble(ecosystem = unique(new_sample_data$ecosystem)) %>%
  mutate(biome = ifelse(ecosystem %in% c('South_Africa', 'Australia', 'Vines', 'Citrus', 'Med'), 'Mediterranean', 'Rainforest'))

new_sample_data <- left_join(new_sample_data, biomes) %>%
  select(sample_info, sample_id, biome, ecosystem, bed, quadrat, sample, quadrat_corner, extraction, pcr) %>%
  column_to_rownames('sample_info') 

sample_data(ps_bact) <- sample_data(new_sample_data)

sample_sums(ps_bact) %>%
  log10() %>%
  hist()

sample_sums(ps_bact) %>%
  sort()

# filter out samples with less than 1000 reads
ps_bact <- prune_samples(sample_sums(ps_bact) > 1000, ps_bact)

sample_sums(ps_bact) %>%
  log10() %>%
  hist()

# remove ASVs where Phylum = NA
ps_bact_sub <- prune_taxa(taxa_names(ps_bact)[tax_table(ps_bact)[, 'Phylum'] != 'NA'], ps_bact)

# look at distribution of read length
rownames(tax_table(ps_bact_sub)) %>%
  nchar() %>%
  hist()

rownames(tax_table(ps_bact_sub)) %>% nchar() %>% summary()

# remove ASVs with read length < 200
ps_bact_sub <- prune_taxa(taxa_names(ps_bact_sub)[nchar(rownames(tax_table(ps_bact_sub))) > 200], ps_bact_sub)

# saved 
saveRDS(ps_bact_sub, 'data/sequencing/processed/phyloseq/ps_bact_sub.rds')

#--------------------------#
# filter for eukaryotes ####
#--------------------------#

# load in bact data
ps_euk <- readRDS('data/sequencing/processed/phyloseq/ps_euk.rds')

# look at taxonomy table for eukaryotes - compare numbers assigned to each level with that from metabar
d_taxa <- tax_table(ps_euk) %>%
  data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'rank', values_to = 'input') %>%
  group_by(rank) %>%
  # count number of values that are not NA
  summarise(named = sum(!is.na(input)),
            n = n(), .groups = 'drop') %>%
  mutate(prop = named/n)
# so few things were assigned to the phylum level

# look at the taxonomy table from the metabar for eukaryotes

# write custom function to get taxa data
get_taxa_data <- function(file){
  d <- readRDS(file)
  d_taxa <- d$motus
  d_taxa <- select(d_taxa, sequence, contains('silva'))
  rownames(d_taxa) <- NULL
  colnames(d_taxa) <- gsub('_silva', '', colnames(d_taxa))
  return(d_taxa)
}

d_taxa_metabar <- map_df(c('data/sequencing/processed/metabar/eden_euc_rep1_postclean_24.rds',
                          'data/sequencing/processed/metabar/eden_euc_rep2_postclean_24.rds',
                          'data/sequencing/processed/metabar/eden_euc_rep3_postclean_24.rds'),
                        get_taxa_data) %>%
  distinct()  %>%
  select(superkingdom:genus) %>%
  pivot_longer(cols = everything(), names_to = 'rank', values_to = 'input') %>%
  group_by(rank) %>%
  # count number of values that are not NA
  summarise(named = sum(!is.na(input)),
            n = n(), .groups = 'drop') %>%
  mutate(prop = named/n)

# load in new sample data
sample_data <- read.csv('data/Eden_PCRs_sample_info.csv') %>%
  filter(!is.na(Ecosystem))

# make new sample data
new_sample_data <- sample_data(ps_euk) %>%
  data.frame()

# any samples in our sample data not in the new sample data?
filter(sample_data,! sample_id %in% new_sample_data$sample_id)
# 5 samples are not in our bacteria samples that were in the full sample info

# create a new sample data object
new_sample_data <- sample_data(ps_euk) %>%
  data.frame() %>%
  rownames_to_column('sample_info') %>%
  # add column for pcr replicate by grabbing the last four characters of sample info
  mutate(pcr = str_sub(sample_info, -1, -1),
         # add column for extraction replicate by grabbing the number between the r and p in sample info
         extraction = str_extract(sample_info, '(?<=r)\\d+(?=p)'),
         extraction = ifelse(is.na(extraction), 1, as.numeric(extraction)),
         # extract all text before the r
         sample = str_extract(sample_info, '.*(?=pcr)'),
         sample = str_extract(sample, '.*(?=r)'),
         # add column for quadrat corner
         quadrat_corner = str_sub(sample, -1, -1),
         quadrat = str_sub(sample, -2, -2),
         bed = str_sub(sample, 1, -3)) %>%
  select(-ecosystem) %>%
  left_join(., select(sample_data, sample_id, ecosystem = Ecosystem))

# unique biome assignment
biomes <- tibble(ecosystem = unique(new_sample_data$ecosystem)) %>%
  mutate(biome = ifelse(ecosystem %in% c('South_Africa', 'Australia', 'Vines', 'Citrus', 'Med'), 'Mediterranean', 'Rainforest'))

new_sample_data <- left_join(new_sample_data, biomes) %>%
  select(sample_info, sample_id, biome, ecosystem, bed, quadrat, sample, quadrat_corner, extraction, pcr) %>%
  column_to_rownames('sample_info') 

sample_data(ps_euk) <- sample_data(new_sample_data)

sample_sums(ps_euk) %>%
  log10() %>%
  hist()

sample_sums(ps_euk) %>%
  sort()

# filter out samples with less than 1000 reads
ps_euk <- prune_samples(sample_sums(ps_euk) > 1000, ps_euk)

sample_sums(ps_euk) %>%
  log10() %>%
  hist()

# remove ASVs where Phylum = NA
ps_euk_sub <- prune_taxa(taxa_names(ps_euk)[tax_table(ps_euk)[, 'Phylum'] != 'NA'], ps_euk)

sample_sums(ps_euk_sub) %>%
  sort()

# look at distribution of read length
rownames(tax_table(ps_euk_sub)) %>%
  nchar() %>%
  hist()

rownames(tax_table(ps_euk_sub)) %>% nchar() %>% summary()

# remove ASVs with read length < 80
ps_euk_sub <- prune_taxa(taxa_names(ps_euk_sub)[nchar(rownames(tax_table(ps_euk_sub))) > 80], ps_euk_sub)

# saved 
saveRDS(ps_euk_sub, 'data/sequencing/processed/phyloseq/ps_euk_sub.rds')



