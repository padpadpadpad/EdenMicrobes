# try and merge the sequence tables across replicates
# each replicate is a different PCR from the same set of samples, that were then run on the different sequencing runs

# want to merge them to be able to ID which PCR the samples came from so that we can use them as a filter later on.
# for example, sum across all PCRs, or say that an OTU must be present in all the PCRs to be included in the analysis

# load packages
librarian::shelf(phyloseq, tidyverse, dada2)

#------------------------------------#
# try and create a sequence table ####
#------------------------------------#

# make function to create seq table
make_seq_table <- function(file){
  d <- readRDS(file)
  d_reads <- d$reads
  colnames(d_reads) <- toupper(d$motus$sequence)
  rownames(d_reads) <- paste(rownames(d_reads), 'pcr', parse_number(file), sep = '')
  return(d_reads)
}

bact1 <- make_seq_table('data/sequencing/processed/metabar/eden_bact_rep1_postclean_24.rds')
bact2 <- make_seq_table('data/sequencing/processed/metabar/eden_bact_rep2_postclean_24.rds')
bact3 <- make_seq_table('data/sequencing/processed/metabar/eden_bact_rep3_postclean_24.rds')

rownames(bact1)
rownames(bact2)
rownames(bact3)

bact_merge <- mergeSequenceTables(tables = list(bact1, bact2, bact3))

euk1 <- make_seq_table('data/sequencing/processed/metabar/eden_euc_rep1_postclean_24.rds')
euk2 <- make_seq_table('data/sequencing/processed/metabar/eden_euc_rep2_postclean_24.rds')
euk3 <- make_seq_table('data/sequencing/processed/metabar/eden_euc_rep3_postclean_24.rds')

euk_merge <- mergeSequenceTables(tables = list(euk1, euk2, euk3))

fungi1 <- make_seq_table('data/sequencing/processed/metabar/eden_fungi_rep1_postclean_24.rds')
fungi2 <- make_seq_table('data/sequencing/processed/metabar/eden_fungi_rep2_postclean_24.rds')
fungi3 <- make_seq_table('data/sequencing/processed/metabar/eden_fungi_rep3_postclean_24.rds')

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

sample_data_bact <- map_df(c('data/sequencing/processed/metabar/eden_bact_rep1_postclean_24.rds',
                             'data/sequencing/processed/metabar/eden_bact_rep2_postclean_24.rds',
                             'data/sequencing/processed/metabar/eden_bact_rep3_postclean_24.rds'),
                           get_sample_data)

sample_data_euk <- map_df(c('data/sequencing/processed/metabar/eden_euc_rep1_postclean_24.rds',
                            'data/sequencing/processed/metabar/eden_euc_rep2_postclean_24.rds',
                            'data/sequencing/processed/metabar/eden_euc_rep3_postclean_24.rds'),
                          get_sample_data)

sample_data_fungi <- map_df(c('data/sequencing/processed/metabar/eden_fungi_rep1_postclean_24.rds',
                              'data/sequencing/processed/metabar/eden_fungi_rep2_postclean_24.rds',
                              'data/sequencing/processed/metabar/eden_fungi_rep3_postclean_24.rds'),
                            get_sample_data)

#-------------------------#
# make phyloseq object ####
#-------------------------#

# load in taxonomy for bacteria
tax_bact <- readRDS('data/sequencing/processed/taxtab/taxtab_bact.rds')

ps_bact <- phyloseq(otu_table(bact_merge, taxa_are_rows = FALSE),
                    sample_data(sample_data_bact),
                    tax_table(tax_bact))

ps_bact

# load in taxonomy for eukaryotes
tax_euk <- readRDS('data/sequencing/processed/taxtab/taxtab_euk.rds')

ps_euk <- phyloseq(otu_table(euk_merge, taxa_are_rows = FALSE),
                   sample_data(sample_data_euk),
                   tax_table(tax_euk))

ps_euk

# load in taxonomy for fungi
tax_fungi <- readRDS('data/sequencing/processed/taxtab/taxtab_fungi.rds')

ps_fungi <- phyloseq(otu_table(fungi_merge, taxa_are_rows = FALSE),
                     sample_data(sample_data_fungi),
                     tax_table(tax_fungi))

ps_fungi

saveRDS(ps_euk, 'data/sequencing/processed/phyloseq/ps_euk.rds')
saveRDS(ps_fungi, 'data/sequencing/processed/phyloseq/ps_fungi.rds')
saveRDS(ps_bact, 'data/sequencing/processed/phyloseq/ps_bact.rds')
