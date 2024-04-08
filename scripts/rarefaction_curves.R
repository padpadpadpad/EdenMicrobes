# ---------------------------
# Purpose of script: To run rarefaction curves on the pooled dataset
#
# What this script does:
# 1. Runs a rarefaction curve for the bacteria, fungi, and eukaryotic datasets
# 2. Creates a plot
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-04-08
#
# Copyright (c) Daniel Padfield, 2024
#
# ---------------------------
#
# Notes:
#
# ---------------------------

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)){
  install.packages("librarian")
}
# load packages
librarian::shelf(phyloseq, reshape2, tidyverse, parallel)

## ---------------------------

# custom rarefaction curve
# taken from this GitHub issue: https://github.com/joey711/phyloseq/issues/143
calculate_rarefaction_curves <- function(psdata, measures, depths, parallel=TRUE) {
  
  # set parallel options if required
  if (parallel) {
    paropts  <- list(.packages=c("phyloseq", "reshape2"))
  } else {
    paropts  <- NULL
  }
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- reshape2::melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- plyr::ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive() && ! parallel, 'text', 'none'), .parallel=parallel, .paropts=paropts)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

# set seed
set.seed(42)

# load in bacterial phyloseq object
ps_bact <- readRDS('data/sequencing/processed/phyloseq/ps_bact_sub.rds')
ps_bact

# grab sample data out
d_samp <- sample_data(ps_bact) %>%
  data.frame()

# aggregate samples by sample
ps_bact2 <- merge_samples(ps_bact, 'sample')

# load in sample data and create a new metadata for the merged phyloseq object
d_samp <- sample_data(ps_bact2) %>%
  data.frame() %>%
  select(quadrat_corner) %>%
  mutate(sample = row.names(.),
         quadrat_corner = str_sub(sample, -1, -1),
         quadrat = str_sub(sample, -2, -2),
         bed = str_sub(sample, 1, -3)) %>%
  select(sample, everything())

sample_data <- read.csv('data/Eden_PCRs_sample_info.csv') %>%
  select(sample_id, ecosystem = Ecosystem) %>%
  mutate(sample = str_extract(sample_id, '.*(?=r)')) %>%
  filter(!is.na(ecosystem)) %>%
  select(-sample_id) %>%
  distinct()

d_samp <- left_join(d_samp, sample_data)
row.names(d_samp) <- d_samp$sample

sample_data(ps_bact2) <- d_samp

# Setting up and registering the cluster
cl <- makeCluster(detectCores(all.tests=TRUE))

doParallel::registerDoParallel(cl)

to_do <- seq(1, max(sample_sums(ps_bact2)), length.out = 40)

# calculate rarefaction curves
# did ten at each number of reads in to-do
rarefaction_curves <- calculate_rarefaction_curves(ps_bact2, c('Observed'), rep(to_do, each = 10), parallel=TRUE)

# summarise the rarefaction curves
rarefaction_curves <- group_by(rarefaction_curves, Depth, Sample, Measure) %>%
  summarise(diversity_mean = mean(Alpha_diversity),
            sd = sd(Alpha_diversity),
            .groups = 'drop') %>%
  janitor::clean_names()

rarefaction_curves <- left_join(rarefaction_curves, sample_data)

rarefaction_curves <- rarefaction_curves %>%
  mutate(biome = ifelse(ecosystem %in% c('South_Africa', 'Australia', 'Vines', 'Citrus', 'Med'), 'Mediterranean', 'Rainforest'))

# plot the rarefaction curves
# plot them
ggplot(filter(rarefaction_curves, measure == 'Observed')) +
  geom_line(aes(depth, diversity_mean, group = sample), key_glyph = 'point') +
  geom_ribbon(aes(x = depth, ymin = diversity_mean - sd, ymax = diversity_mean + sd, group = sample), alpha = 0.3) +
  theme_bw(base_size = 14) +
  labs(x = 'Sequencing depth (Number of reads)',
       y = 'Estimated species richness') +
  facet_wrap(~biome + ecosystem, ncol = 5) +
  NULL

ggsave('plots/rarefaction_curves_bact.png', width = 12, height = 6)
