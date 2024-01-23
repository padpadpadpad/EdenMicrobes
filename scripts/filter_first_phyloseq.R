# first filter and quick and dirty plot
# try with fungi first

# load in packages
librarian::shelf(phyloseq, tidyverse)

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
