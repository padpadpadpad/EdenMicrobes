# load in fungal filtered phyloseq object

ps_fungi <- readRDS('data/sequencing/processed/phyloseq/ps_fungi_sub.rds')

ps_fungi

sample_data <- read.csv('data/Eden_PCRs_sample_info.csv') %>%
  select(sample_id, ecosystem = Ecosystem) %>%
  filter(!is.na(ecosystem))

# create a distance matrix for the fungi
dist_bray <- distance(ps_fungi, method = 'bray')

# create sample data for the fungi
d_samp <- sample_data(ps_fungi) %>%
  data.frame()

filter(sample_data,! sample_id %in% d_samp$sample_id)

group_by(d_samp, sample_id) %>%
  tally() %>%
  View()
