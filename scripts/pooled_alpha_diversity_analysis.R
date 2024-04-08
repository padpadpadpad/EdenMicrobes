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
librarian::shelf(adw96/breakaway, phyloseq, microbiome, tidyverse, furrr, lme4, patchwork)

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

ps_bact
ps_bact2

# look at sequencing numbers
ps_bact2 %>%
  phyloseq::sample_sums() %>%
  sort()

# ok there are 4 samples with under 15000 reads, so lets drop those and rarefy to 16000 reads

#-----------------------------#
# alpha diversity analysis ####
#-----------------------------#

# there are a bunch of ways to do this, some people think rarefying is crap, some people think rarefying is the best way to go. Both of these sets of people are smart. We shall use a couple of different approaches and try and do the analysis

# 1. use breakaway to estimate species richness per sample
d_break <- breakaway(ps_bact2)
plot(d_break)

# make it into a tibble
d_break2 <- summary(d_break) %>% as_tibble()
head(d_break2)

d_break2 <- left_join(rename(d_break2, sample = sample_names), d_samp) %>%
  group_by(ecosystem, quadrat, bed) %>%
  mutate(id = as.character(cur_group_id())) %>%
  ungroup() %>%
  mutate(id2 = paste(id, quadrat_corner, sep = '_'),
         biome = ifelse(ecosystem %in% c('South_Africa', 'Australia', 'Vines', 'Citrus', 'Med'), 'Mediterranean', 'Rainforest'))

# make a plot
ggplot(d_break2_sub, aes(id, log(estimate))) +
  geom_point(aes(group = quadrat_corner), position = position_dodge(width = 0.4)) +
  geom_linerange(aes(ymin = log(lower), ymax = log(upper), group = quadrat_corner), position = position_dodge(width = 0.4)) +
  facet_wrap(~ecosystem, scales = 'free_x') +
  theme_bw()
# a few of these look absolutely bonkers in terms of estimate and uncertainty
# would recommend checking these
# Amazon 1_1, Australia 6_1 and 6_3, Citrus 12_2, Cocoa 13_3, Med 19_1 and 22_1, South Africa 26_2, Vines 28_2, West_Africa 32_2

d_break2_sub <- d_break2 %>%
  filter(!id2 %in% c('1_2', '6_1', '6_3', '12_2', '13_3', '19_1', '22_1', '26_2', '28_2', '32_2')) %>%
  mutate(id_char = as.character(id))

# run statistical analysis 

# want to check whether biome or ecosystems have a significant effect on species richness
# random effect is quadrat within bed, then there are 4 points for each bed
# this is id

# use betta to look at absolute impact of biome
mod1 <- betta_random(formula = estimate ~ 0 + biome | id_char, 
              ses = error,
              data = d_break2_sub)
mod1$table
# it says there is a significant difference between ecosystems but I do not believe it

d_break2_sub <- d_break2_sub %>%
  mutate(weights = 1/(error^2),
         estimate_round = round(estimate, 0))

mod1.2 <- glmer(estimate_round ~ biome + (1|id_char),
                family = 'poisson',
                d_break2_sub,
                weights = 1/error)

summary(mod1.2)
emmeans::emmeans(mod1.2, pairwise ~ biome, type = 'response')

# create empty dataframe
betta_estimates1 <- data.frame(biome = unique(d_break2_sub$biome),
                              estimates = NA,
                              standard_errors = NA,
                              lower = NA,
                              upper = NA, ecosystem = 'all')

# set vector for linear model matrix where everything is switched off
linear_com <- rep(0, times = nrow(betta_estimates1))

for(i in 1:nrow(betta_estimates1)){
  temp_linear_com <- linear_com
  temp_linear_com[i] <- 1
  temp <- betta_lincom(mod1, linear_com = temp_linear_com)
  
  betta_estimates1$estimates[i] <- temp$Estimates
  betta_estimates1$standard_errors[i] <- temp$`Standard Errors`
  betta_estimates1$lower[i] <- temp$`Lower CIs`
  betta_estimates1$upper[i] <- temp$`Upper CIs`
}


# look for variation between ecosystems
mod2 <- betta_random(formula = estimate ~ 0 + ecosystem | id_char, 
                     ses = error,
                     data = d_break2_sub)
mod2$table
mod2$ssq_group

# get confidence interval of each estimate using betta_lincom

# create empty dataframe
betta_estimates2 <- data.frame(ecosystem = unique(d_break2_sub$ecosystem),
                              estimates = NA,
                              standard_errors = NA,
                              lower = NA,
                              upper = NA)

# set vector for linear model matrix where everything is switched off
linear_com <- rep(0, times = nrow(betta_estimates2))

for(i in 1:nrow(betta_estimates2)){
  temp_linear_com <- linear_com
  temp_linear_com[i] <- 1
  temp <- betta_lincom(mod2, linear_com = temp_linear_com)
  
  betta_estimates2$estimates[i] <- temp$Estimates
  betta_estimates2$standard_errors[i] <- temp$`Standard Errors`
  betta_estimates2$lower[i] <- temp$`Lower CIs`
  betta_estimates2$upper[i] <- temp$`Upper CIs`
}

betta_estimates2 <- betta_estimates2 %>% 
  mutate(biome = ifelse(ecosystem %in% c('South_Africa', 'Australia', 'Vines', 'Citrus', 'Med'), 'Mediterranean', 'Rainforest'))

betta_estimates <- bind_rows(betta_estimates1, betta_estimates2) %>%
  mutate(ecosystem2 = gsub('_', ' ', ecosystem))

p1 <- ggplot(betta_estimates, aes(ecosystem2, estimates, ymin = lower, ymax = upper)) +
  geom_pointrange(size = 1.2) +
  theme_bw() +
  facet_wrap(~biome, scales = 'free_x') +
  theme_bw(base_size = 14) +
  labs(x = 'Ecosystem',
       y = 'Richeness',
       title = 'Richness estimates through frequency ratios.') +
  scale_x_discrete(labels = scales::wrap_format(10))

# 2. rarefy the data and calculate richness

# remove samples with less than 15000 reads
ps_bact2 <- prune_samples(sample_sums(ps_bact2) > 15000, ps_bact2)

# setup dataframe to store results
n_boots <- 100

d_rarefy <- tibble(iter = 1:n_boots)

# create function to rarefy data and calculate richness
rarefy_richness <- function(ps, sample_size){
  temp <- rarefy_even_depth(ps, sample.size = sample_size)
  temp <- estimate_richness(temp, split = TRUE) %>%
    rownames_to_column(var = 'sample')
  return(temp)
}

plan(multisession, workers = 4)

d_rarefy <- mutate(d_rarefy, richness = furrr::future_map(iter, ~rarefy_richness(ps_bact2, 15000), .progress = TRUE))

d_rarefy <- unnest(d_rarefy, richness)

d_rarefy_sum <- group_by(d_rarefy, sample) %>%
  summarise(mean_richness = mean(Observed), .groups = 'drop') %>%
  mutate(mean_richness_whole = round(mean_richness, 0))

d_rarefy_sum <- left_join(d_rarefy_sum, d_samp) %>%
  group_by(ecosystem, quadrat, bed) %>%
  mutate(id = as.character(cur_group_id())) %>%
  ungroup() %>%
  mutate(id2 = paste(id, quadrat_corner, sep = '_'),
         biome = ifelse(ecosystem %in% c('South_Africa', 'Australia', 'Vines', 'Citrus', 'Med'), 'Mediterranean', 'Rainforest'))
  
# quick plot
ggplot(d_rarefy_sum, aes(id, mean_richness)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ecosystem, scales = 'free_x')

# right lets look at modelling this using non-linear mixed models
mod3 <- glmer(mean_richness_whole ~ biome + (1|ecosystem/id), family = 'poisson', data = d_rarefy_sum)
mod4 <- glmer(mean_richness_whole ~ 1 + (1|ecosystem/id), family = 'poisson', data = d_rarefy_sum)
summary(mod3)

AIC(mod3, mod4)

emmeans::emmeans(mod3, pairwise ~ biome, type = 'response')

d_preds_biome <- emmeans::emmeans(mod3, pairwise ~ biome, type = 'response')$emmeans %>%
  data.frame() %>%
  mutate(ecosystem = 'all') %>%
  select(biome, ecosystem, pred = rate, lower = asymp.LCL, upper = asymp.UCL)

# get estimates for individual ecosystems by adding the random effect of each level of the random effect onto the fixed effect
ranef(mod3) %>% lattice::dotplot()

d_preds <- select(d_rarefy_sum, biome, ecosystem) %>%
  distinct()

new_predict <- function(x){predict(x, newdata = d_preds, re.form = ~1|ecosystem, type = 'response')}

preds_boot <- bootMer(mod3, new_predict, re.form = ~1|ecosystem, nsim = 1000) %>%
  confint(.)

d_preds <- mutate(d_preds, pred = predict(mod3, newdata = d_preds, re.form = ~1|ecosystem, type = 'response')) %>%
  bind_cols(., preds_boot) %>%
  rename(lower = `2.5 %`, upper = `97.5 %`) %>%
  bind_rows(d_preds_biome) %>%
  mutate(ecosystem2 = gsub('_', ' ', ecosystem))

p2 <- ggplot(d_preds, aes(ecosystem2, pred, ymin = lower, ymax = upper)) +
  geom_pointrange(size = 1.2) +
  theme_bw() +
  facet_wrap(~biome, scales = 'free_x') +
  scale_x_discrete(labels = scales::wrap_format(10)) +
  theme_bw(base_size = 14) +
  labs(x = 'Ecosystem',
       y = 'Richeness',
       title = 'Richness estimates through rarefaction.')

p1 + p2 +
  plot_layout(ncol = 1)

# save plot
ggsave('plots/richness_estimates.png', height = 8, width = 10)
