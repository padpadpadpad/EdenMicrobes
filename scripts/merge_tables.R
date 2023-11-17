# try and merge the sequence tables across replicates

# load packages
librarian::shelf(phyloseq, tidyverse, dada2)

d_example <- readRDS('data/example_seqtab.rds')
rownames(d_example)
colnames(d_example)

# load in example dataset
d <- readRDS('data/processed/eden_bact_rep1_postclean.rds')
d2 <- readRDS('data/processed/eden_bact_rep2_postclean.rds')
d3 <- readRDS('data/processed/eden_bact_rep3_postclean.rds')


d_samps <- d$samples
d2_samps <- d$samples

# try and create a sequence table
d

# make function to create seq table
make_seq_table <- function(file){
  d <- readRDS(file)
  d_reads <- d$reads
  colnames(d_reads) <- toupper(d$motus$sequence)
  rownames(d_reads) <- paste(rownames(d_reads), 'pcr', parse_number(file), sep = '')
  return(d_reads)
}

bact1 <- make_seq_table('data/processed/eden_bact_rep1_postclean.rds')
bact2 <- make_seq_table('data/processed/eden_bact_rep2_postclean.rds')
bact3 <- make_seq_table('data/processed/eden_bact_rep3_postclean.rds')

rownames(bact1)
rownames(bact2)
rownames(bact3)

bact_merge <- mergeSequenceTables(tables = list(bact1, bact2, bact3))

euk1 <- make_seq_table('data/processed/eden_euc_rep1_postclean.rds')
euk2 <- make_seq_table('data/processed/eden_euc_rep2_postclean.rds')
euk3 <- make_seq_table('data/processed/eden_euc_rep3_postclean.rds')

euk_merge <- mergeSequenceTables(tables = list(euk1, euk2, euk3))
