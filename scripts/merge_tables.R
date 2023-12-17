# try and merge the sequence tables across replicates
# each replicate is a different PCR from the same set of samples, that were then run on the different sequencing runs

# want to merge them to be able to ID which PCR the samples came from so that we can use them as a filter later on.
# for example, sum across all PCRs, or say that an OTU must be present in all the PCRs to be included in the analysis

# load packages
librarian::shelf(phyloseq, tidyverse, dada2)

d_example <- readRDS('data/sequencing/example_seqtab.rds')
rownames(d_example)
colnames(d_example)

# load in example dataset
d <- readRDS('data/sequencing/processed/eden_bact_rep1_postclean.rds')
d2 <- readRDS('data/sequencing/processed/eden_bact_rep2_postclean.rds')
d3 <- readRDS('data/sequencing/processed/eden_bact_rep3_postclean.rds')

d_samps <- d$samples
d2_samps <- d2$samples

#------------------------------------#
# try and create a sequence table ####
#------------------------------------#

d

# make function to create seq table
make_seq_table <- function(file){
  d <- readRDS(file)
  d_reads <- d$reads
  colnames(d_reads) <- toupper(d$motus$sequence)
  rownames(d_reads) <- paste(rownames(d_reads), 'pcr', parse_number(file), sep = '')
  return(d_reads)
}

bact1 <- make_seq_table('data/sequencing/processed/eden_bact_rep1_postclean.rds')
bact2 <- make_seq_table('data/sequencing/processed/eden_bact_rep2_postclean.rds')
bact3 <- make_seq_table('data/sequencing/processed/eden_bact_rep3_postclean.rds')

rownames(bact1)
rownames(bact2)
rownames(bact3)

bact_merge <- mergeSequenceTables(tables = list(bact1, bact2, bact3))

euk1 <- make_seq_table('data/sequencing/processed/eden_euc_rep1_postclean.rds')
euk2 <- make_seq_table('data/sequencing/processed/eden_euc_rep2_postclean.rds')
euk3 <- make_seq_table('data/sequencing/processed/eden_euc_rep3_postclean.rds')

euk_merge <- mergeSequenceTables(tables = list(euk1, euk2, euk3))

fungi1 <- make_seq_table('data/sequencing/processed/eden_fungi_rep1_postclean.rds')
fungi2 <- make_seq_table('data/sequencing/processed/eden_fungi_rep2_postclean.rds')
fungi3 <- make_seq_table('data/sequencing/processed/eden_fungi_rep3_postclean.rds')

fungi_merge <- mergeSequenceTables(tables = list(fungi1, fungi2, fungi3))

#---------------------------------------#
# try and create a sample data table ####
#---------------------------------------#

# make a function to extract the sample data from the metabar output

get_sample_data <- function(file){
  d <- readRDS(file)
  d_samples <- d$samples
  # if Ecosystem column is not present then add it
  if(!'Ecosystem' %in% colnames(d_samples)){
    d_samples$Ecosystem <- NA
  }
  d_samples <- dplyr::select(d_samples, sample_id, ecosystem = Ecosystem)
  rownames(d_samples) <- paste(rownames(d_samples), 'pcr', parse_number(file), sep = '')
  return(d_samples)
}

sample_data_bact <- map_df(c('data/sequencing/processed/eden_bact_rep1_postclean.rds',
                             'data/sequencing/processed/eden_bact_rep2_postclean.rds',
                             'data/sequencing/processed/eden_bact_rep3_postclean.rds'),
                           get_sample_data)

sample_data_euk <- map_df(c('data/sequencing/processed/eden_euc_rep1_postclean.rds',
                             'data/sequencing/processed/eden_euc_rep2_postclean.rds',
                             'data/sequencing/processed/eden_euc_rep3_postclean.rds'),
                           get_sample_data)

sample_data_fungi <- map_df(c('data/sequencing/processed/eden_fungi_rep1_postclean.rds',
                              'data/sequencing/processed/eden_fungi_rep2_postclean.rds',
                              'data/sequencing/processed/eden_fungi_rep3_postclean.rds'),
                            get_sample_data)

#-------------------------------------#
# try and create a taxonomic table ####
#-------------------------------------#

# make a function to extract the sample data from the metabar output

get_taxa_data <- function(file){
  d <- readRDS(file)
  d_taxa <- d$motus
  d_taxa <- select(d_taxa, sequence, contains('silva'))
  rownames(d_taxa) <- NULL
  colnames(d_taxa) <- gsub('_silva', '', colnames(d_taxa))
  return(d_taxa)
}

taxa_data_bact <- map_df(c('data/sequencing/processed/eden_bact_rep1_postclean.rds',
                             'data/sequencing/processed/eden_bact_rep2_postclean.rds',
                             'data/sequencing/processed/eden_bact_rep3_postclean.rds'),
                           get_taxa_data) %>%
  distinct() %>%
  mutate(sequence = toupper(sequence))

otus_to_check <- group_by(taxa_data_bact, sequence) %>%
  tally() %>%
  filter(n > 1)

# some exact sequences have different taxonomy
to_check <- filter(taxa_data_bact, sequence %in% otus_to_check$sequence) %>%
  group_by(sequence) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  select(group_id, everything())

# save out these sequences to check
write.csv(to_check, 'data/sequencing/processed/eden_bact_to_check.csv', row.names = FALSE)

# save out list of unique sequences to reassign taxonomy to
select(taxa_data_bact, sequence) %>%
  distinct() %>%
  saveRDS('data/sequencing/processed/eden_bact_unique_sequence.rds')

taxa_data_euk <- map_df(c('data/sequencing/processed/eden_euc_rep1_postclean.rds',
                            'data/sequencing/processed/eden_euc_rep2_postclean.rds',
                            'data/sequencing/processed/eden_euc_rep3_postclean.rds'),
                          get_taxa_data) %>%
  distinct()  %>%
  mutate(sequence = toupper(sequence)) 

otus_to_check <- group_by(taxa_data_euk, sequence) %>%
  tally() %>%
  filter(n > 1)

# some exact sequences have different taxonomy
to_check <- filter(taxa_data_euk, sequence %in% otus_to_check$sequence) %>%
  group_by(sequence) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  select(group_id, everything())

group_by(taxa_data_euk, sequence) %>%
  tally() %>%
  filter(n > 1)

# save out these sequences to check
write.csv(to_check, 'data/sequencing/processed/eden_euk_to_check.csv', row.names = FALSE)

# save out list of unique sequences to reassign taxonomy to
select(taxa_data_euk, sequence) %>%
  distinct() %>%
  saveRDS('data/sequencing/processed/eden_euk_unique_sequence.rds')

get_taxa_data <- function(file){
  d <- readRDS(file)
  d_taxa <- d$motus
  d_taxa <- select(d_taxa, sequence, contains('name'))
  rownames(d_taxa) <- NULL
  colnames(d_taxa) <- gsub('_name', '', colnames(d_taxa))
  return(d_taxa)
}

# do the same with fungi
taxa_data_fungi <- map_df(c('data/sequencing/processed/eden_fungi_rep1_postclean.rds',
                            'data/sequencing/processed/eden_fungi_rep2_postclean.rds',
                            'data/sequencing/processed/eden_fungi_rep3_postclean.rds'),
                          get_taxa_data) %>%
  distinct()  %>%
  mutate(sequence = toupper(sequence))

otus_to_check <- group_by(taxa_data_fungi, sequence) %>%
  tally() %>%
  filter(n > 1)

# some exact sequences have different taxonomy
to_check <- filter(taxa_data_fungi, sequence %in% otus_to_check$sequence) %>%
  group_by(sequence) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  select(group_id, everything())

group_by(taxa_data_fungi, sequence) %>%
  tally() %>%
  filter(n > 1)

# save out these sequences to check
write.csv(to_check, 'data/sequencing/processed/eden_fungi_to_check.csv', row.names = FALSE)

# save out list of unique sequences to reassign taxonomy to
select(taxa_data_fungi, sequence) %>%
  distinct() %>%
  saveRDS('data/sequencing/processed/eden_fungi_unique_sequence.rds')

# some exact sequences have different taxonomy

#-------------------------#
# make phyloseq object ####
#-------------------------#

ps_bact <- phyloseq(otu_table(bact_merge, taxa_are_rows = FALSE),
                    sample_data(sample_data_bact))


