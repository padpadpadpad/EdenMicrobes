#----------------------------------------------------#
# Try and create a first cluster plot of the dataset #
#----------------------------------------------------#

# load in packages
librarian::shelf(monochromeR, padpadpadpad/MicrobioUoE, phyloseq, tidyverse)

# function to make points lighter but not transparent
alpha_hex <- function(hex_code, alpha = 1){
  # split by ',' and extract number
  temp <- monochromeR::hex_to_rgb(hex_code) %>%
    str_split(., ',') %>%
    .[[1]] %>%
    parse_number()
  
  return(rgba_to_hex(c(temp, alpha)))
}

#-----------------------#
# first for bacteria ####
#-----------------------#

# load in bacteria filtered phyloseq object
ps_bact <- readRDS('data/sequencing/processed/phyloseq/ps_bact_sub.rds')

ps_bact

# read in sample data from the PCR sample info
sample_data <- read.csv('data/Eden_PCRs_sample_info.csv') %>%
  select(sample_id, ecosystem = Ecosystem) %>%
  filter(!is.na(ecosystem))

# check sample data
d_samp <- sample_data(ps_bact) %>%
  data.frame()

filter(sample_data,! sample_id %in% d_samp$sample_id)
# 8 samples are not present but that is ok

group_by(d_samp, sample_id) %>%
  tally() %>%
  View()

# transform counts to relative abundances for ordination
ps_prop <- transform_sample_counts(ps_bact, function(x){x / sum(x)})

# create a distance matrix for the bacteria using the relative abundances

# if 'data/sequencing/processed/phyloseq/bact_dist_bray.rds' exists, read it in, if not, run the code

if(file.exists('data/sequencing/processed/phyloseq/bact_dist_bray.rds')){
  dist_bray <- readRDS('data/sequencing/processed/phyloseq/bact_dist_bray.rds')
} else {
  dist_bray <- distance(ps_prop, method = 'bray')
  saveRDS(dist_bray, 'data/sequencing/processed/phyloseq/bact_dist_bray.rds')
}

# create a pcoa data set - choose only eigenvalues that are > 0
#https://stackoverflow.com/questions/8924488/applying-the-pvclust-r-function-to-a-precomputed-dist-object#27148408
# error message tells how may eigenvalues are > 0
d_pcoa_correct <- ape::pcoa(dist_bray, correction = 'cailliez') 

# correct eigenvalues
correct_eigenvalues <- ape::pcoa(dist_bray, correction = 'cailliez') %>%
  .$values %>%
  janitor::clean_names() %>%
  mutate(axis = 1:n()) %>%
  mutate(axis = paste('axis', axis, sep = '_'))

# how many axes to explain 50% of the data
filter(correct_eigenvalues, cumul_eig < 0.5) %>%
  nrow()
# 38

# plot scree plot for every eigenvalue worth 1%
filter(correct_eigenvalues, relative_eig > 0.01) %>%
  ggplot(., aes(x = parse_number(axis), y = relative_eig)) +
  geom_col() +
  labs(x = 'PCoA', y = 'Proportion of Variance Explained')

# grab out the individual points of the pcoa
d_pcoa <- d_pcoa_correct$vectors %>%
  data.frame() %>%
  janitor::clean_names()

# check rownames of d_pcoa are the same as in d_samp
all(rownames(d_pcoa) == rownames(d_samp))

# add in the sample data
d_samples <- cbind(d_samp, d_pcoa)

# make long format
d_samples_long <- d_samples %>%
  select(-sample_id) %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value')


# calculate extraction centroids based on a each extraction ####

d_ext_centroids <- d_samples %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner, extraction) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

# make long format
d_ext_centroids_long <- d_ext_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_ext') 

# calculate distance of each PCR to its respective centroid
d_distance_ext_centroids <- left_join(d_samples_long, d_ext_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, relative_eig)) %>%
  filter(parse_number(axis) <= 38) %>%
  mutate(dist = relative_eig*(value - value_ext)^2) %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner, extraction, pcr) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_ext_centroids$dist)

# create dataset for plotting lines between raw points and centroids
d_ext_points_2_centroid <- left_join(select(d_samples, biome, ecosystem, bed, quadrat, quadrat_corner, extraction, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2), select(d_ext_centroids, biome, ecosystem, bed, quadrat, quadrat_corner, extraction, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2))

# calculate the quadrat corner centroids from the extraction centroids ####

d_corner_centroids <- d_ext_centroids %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

d_corner_centroids_long <- d_corner_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_corner')

# calculate distance of each PCR to its respective centroid
d_distance_corner_centroids <- left_join(d_ext_centroids_long, d_corner_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, relative_eig)) %>%
  filter(parse_number(axis) <= 38) %>%
  mutate(dist = relative_eig*(value_corner - value_ext)^2) %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner, extraction) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_corner_centroids$dist)

# calculate distance between quadrat corner centroids and extraction centroids
d_corner_points_2_centroid <- left_join(select(d_corner_centroids, biome, ecosystem, bed, quadrat, quadrat_corner, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), select(d_ext_centroids, biome, ecosystem, bed, quadrat, quadrat_corner, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2))

# calculate quadrat centroids ####

d_quad_centroids <- d_corner_centroids %>%
  group_by(biome, ecosystem, bed, quadrat) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

d_quad_centroids_long <- d_quad_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_quad')

# calculate distance of each PCR to its respective centroid
d_distance_quad_centroids <- left_join(d_corner_centroids_long, d_quad_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, relative_eig)) %>%
  filter(parse_number(axis) <= 38) %>%
  mutate(dist = relative_eig*(value_quad - value_corner)^2) %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_quad_centroids$dist)

# calculate distance between quadrat centroids and quadrat corner centroids
d_quad_points_2_centroid <- left_join(select(d_quad_centroids, biome, ecosystem, bed, quadrat, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), select(d_corner_centroids, biome, ecosystem, bed, quadrat, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2))

# calculate bed metrics ####
d_bed_centroids <- d_quad_centroids %>%
  group_by(biome, ecosystem, bed) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

d_bed_centroids_long <- d_bed_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_bed')

# calculate distance of each quad to its respective centroid
d_distance_bed_centroids <- left_join(d_quad_centroids_long, d_bed_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, relative_eig)) %>%
  filter(parse_number(axis) <= 38) %>%
  mutate(dist = relative_eig*(value_bed - value_quad)^2) %>%
  group_by(biome, ecosystem, bed, quadrat) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_bed_centroids$dist)

# calculate distance between bed centroids and quadrat centroids
d_bed_points_2_centroid <- left_join(select(d_bed_centroids, biome, ecosystem, bed, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), select(d_quad_centroids, biome, ecosystem, bed, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2)) %>%
  filter(ecosystem %in% c('Australia', 'South_Africa'))

# calculate ecosystem centroids ####

# which ecosystems have multiple beds?
d_samp %>% select(ecosystem, bed) %>%
  distinct() %>%
  group_by(ecosystem) %>%
  tally()
# australia and south africa

# create centroids from the quadrats when there is only a single bed, and from the bed centroids when there are multiple beds
d_eco_centroids <- bind_rows(filter(d_quad_centroids, !ecosystem %in% c('Australia', 'South_Africa')),
                             filter(d_bed_centroids, ecosystem %in% c('Australia', 'South_Africa'))) %>%
  group_by(biome, ecosystem) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

# make long
d_eco_centroids_long <- d_eco_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_eco')

# calculate distance of each quadrat centroid (or bed centroid) from its respective centroid
d_distance_eco_centroids <- left_join(bind_rows(filter(d_quad_centroids_long, !ecosystem %in% c('Australia', 'South_Africa')),
                                                filter(d_bed_centroids_long, ecosystem %in% c('Australia', 'South_Africa')) %>% rename(value_quad = value_bed)), d_eco_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, relative_eig)) %>%
  filter(parse_number(axis) <= 38) %>%
  mutate(dist = relative_eig*(value_eco - value_quad)^2) %>%
  group_by(biome, ecosystem, bed, quadrat) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_eco_centroids$dist)

# calculate distance between ecosystem centroids and quadrat centroids
d_eco_points_2_centroid <- left_join(select(d_eco_centroids, biome, ecosystem, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), bind_rows(filter(d_quad_centroids, !ecosystem %in% c('Australia', 'South_Africa')) %>%
                                                                                                                                             select(biome, ecosystem, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2),
                                                                                                                                           filter(d_bed_centroids, ecosystem %in% c('Australia', 'South_Africa')) %>%
                                                                                                                                             select(biome, ecosystem, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2)))

# calculate biome centroids ####

# create centroids from the ecosystems
d_biome_centroids <- d_eco_centroids %>%
  group_by(biome) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

# make long
d_biome_centroids_long <- d_biome_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_biome')

# calculate distance of each ecosystem centroid from its respective centroid
d_distance_biome_centroids <- left_join(d_eco_centroids_long, d_biome_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, relative_eig)) %>%
  filter(parse_number(axis) <= 38) %>%
  mutate(dist = relative_eig*(value_biome - value_eco)^2) %>%
  group_by(biome, ecosystem) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_biome_centroids$dist)

# set colours for plot
colours <- RColorBrewer::brewer.pal(11, 'Spectral')

# make plot - pretty gnarly but think THINK it works
ggplot() +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_ext_points_2_centroid, alpha = 0.1) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_corner_points_2_centroid, alpha = 0.3) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_quad_points_2_centroid, alpha = 0.3) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_bed_points_2_centroid, alpha = 0.3) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_eco_points_2_centroid, alpha = 0.7) +
  geom_point(aes(axis_1, axis_2), d_biome_centroids, size = 10) +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_eco_centroids, size = 8) + 
  scale_colour_manual(values = colours) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), filter(d_bed_centroids, ecosystem %in% c('Australia', 'South Africa')), size = 6, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.8)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_quad_centroids, size = 5, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.6)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_corner_centroids, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.4)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_ext_centroids, size = 1.5, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.2)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), size = 0.5, d_samples, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.2)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = 'Axis 1 (12.6%)',
       y = 'Axis 2 (5%)') +
  #facet_wrap(~ecosystem, scales = 'free') +
  NULL

ggsave('plots/pcoa_all_bacteria.png', width = 12, height = 8)

# combine the distances into one dataframe
d_distance_all <- bind_rows(
  select(d_distance_ext_centroids, biome, ecosystem, dist) %>% mutate(type = 'PCR',
                                                                     order = 1),
  select(d_distance_corner_centroids, biome, ecosystem, dist) %>% mutate(type = 'Extraction',
                                                                     order = 2),
  select(d_distance_quad_centroids, biome, ecosystem, dist) %>% mutate(type = 'Quadrat corner',
                                                                     order = 4),
  select(d_distance_bed_centroids, biome, ecosystem, dist) %>% mutate(type = 'Quadrat',
                                                                      order = 5),
  select(d_distance_eco_centroids, biome, ecosystem, dist) %>% mutate(type = 'Bed',
                                                                     order = 6),
  select(d_distance_biome_centroids, biome, ecosystem, dist) %>% mutate(type = 'Ecosystem',
                                                                     order = 7)
)

ggplot(d_distance_all, aes(forcats::fct_reorder(type, order), dist)) +
  geom_pretty_boxplot(aes(fill = type, col = type)) +
  geom_point(aes(col = type), shape = 21, fill = 'white', position = position_jitter(width = 0.3)) +
  theme_bw(base_size = 12) +
  labs(x = 'Level of replication',
       y = 'Distance between replicates at that level') +
  scale_color_brewer('Level of replication', palette = 'Set1', type = 'qual') +
  scale_fill_brewer('Level of replication', palette = 'Set1', type = 'qual') +
  facet_wrap(~biome) +
  NULL

ggsave('plots/distance_all_bacteria.png', width = 8, height = 4)

#---------------------#
# second for fungi ####
#---------------------#

# load in bacteria filtered phyloseq object
ps_fungi <- readRDS('data/sequencing/processed/phyloseq/ps_fungi_sub.rds')

ps_fungi

# read in sample data from the PCR sample info
sample_data <- read.csv('data/Eden_PCRs_sample_info.csv') %>%
  select(sample_id, ecosystem = Ecosystem) %>%
  filter(!is.na(ecosystem))

# check sample data
d_samp <- sample_data(ps_fungi) %>%
  data.frame()

filter(sample_data,! sample_id %in% d_samp$sample_id)
# 77 samples are not present but that is ok

group_by(d_samp, sample_id) %>%
  tally() %>%
  View()

# transform counts to relative abundances for ordination
ps_prop <- transform_sample_counts(ps_fungi, function(x){x / sum(x)})

# create a distance matrix for the bacteria using the relative abundances

# if 'data/sequencing/processed/phyloseq/bact_dist_bray.rds' exists, read it in, if not, run the code

if(file.exists('data/sequencing/processed/phyloseq/fungi_dist_bray.rds')){
  dist_bray <- readRDS('data/sequencing/processed/phyloseq/fungi_dist_bray.rds')
} else {
  dist_bray <- distance(ps_prop, method = 'bray')
  saveRDS(dist_bray, 'data/sequencing/processed/phyloseq/fungi_dist_bray.rds')
}

# create a pcoa data set - choose only eigenvalues that are > 0
#https://stackoverflow.com/questions/8924488/applying-the-pvclust-r-function-to-a-precomputed-dist-object#27148408
# error message tells how may eigenvalues are > 0
d_pcoa_correct <- ape::pcoa(dist_bray, correction = 'cailliez') 

# correct eigenvalues
correct_eigenvalues <- ape::pcoa(dist_bray, correction = 'cailliez') %>%
  .$values %>%
  janitor::clean_names() %>%
  mutate(axis = 1:n()) %>%
  mutate(axis = paste('axis', axis, sep = '_'))

# how many axes to explain 50% of the data
num_axes <- filter(correct_eigenvalues, cum_corr_eig < 0.5) %>%
  nrow()
# 22

# plot scree plot for every eigenvalue worth 1%
filter(correct_eigenvalues, rel_corr_eig > 0.01) %>%
  ggplot(., aes(x = parse_number(axis), y = rel_corr_eig)) +
  geom_col() +
  labs(x = 'PCoA', y = 'Proportion of Variance Explained')

# grab out the individual points of the pcoa
d_pcoa <- d_pcoa_correct$vectors %>%
  data.frame() %>%
  janitor::clean_names()

# check rownames of d_pcoa are the same as in d_samp
all(rownames(d_pcoa) == rownames(d_samp))

# add in the sample data
d_samples <- cbind(d_samp, d_pcoa)

# make long format
d_samples_long <- d_samples %>%
  select(-sample_id) %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value')

# calculate extraction centroids based on a each extraction ####

d_ext_centroids <- d_samples %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner, extraction) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

# make long format
d_ext_centroids_long <- d_ext_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_ext') 

# calculate distance of each PCR to its respective centroid
d_distance_ext_centroids <- left_join(d_samples_long, d_ext_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value - value_ext)^2) %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner, extraction, pcr) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_ext_centroids$dist)

# create dataset for plotting lines between raw points and centroids
d_ext_points_2_centroid <- left_join(select(d_samples, biome, ecosystem, bed, quadrat, quadrat_corner, extraction, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2), select(d_ext_centroids, biome, ecosystem, bed, quadrat, quadrat_corner, extraction, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2))

# calculate the quadrat corner centroids from the extraction centroids ####

d_corner_centroids <- d_ext_centroids %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

d_corner_centroids_long <- d_corner_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_corner')

# calculate distance of each PCR to its respective centroid
d_distance_corner_centroids <- left_join(d_ext_centroids_long, d_corner_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_corner - value_ext)^2) %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner, extraction) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_corner_centroids$dist)

# calculate distance between quadrat corner centroids and extraction centroids
d_corner_points_2_centroid <- left_join(select(d_corner_centroids, biome, ecosystem, bed, quadrat, quadrat_corner, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), select(d_ext_centroids, biome, ecosystem, bed, quadrat, quadrat_corner, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2))

# calculate quadrat centroids ####

d_quad_centroids <- d_corner_centroids %>%
  group_by(biome, ecosystem, bed, quadrat) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

d_quad_centroids_long <- d_quad_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_quad')

# calculate distance of each PCR to its respective centroid
d_distance_quad_centroids <- left_join(d_corner_centroids_long, d_quad_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_quad - value_corner)^2) %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_quad_centroids$dist)

# calculate distance between quadrat centroids and quadrat corner centroids
d_quad_points_2_centroid <- left_join(select(d_quad_centroids, biome, ecosystem, bed, quadrat, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), select(d_corner_centroids, biome, ecosystem, bed, quadrat, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2))

# calculate bed metrics ####
d_bed_centroids <- d_quad_centroids %>%
  group_by(biome, ecosystem, bed) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

d_bed_centroids_long <- d_bed_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_bed')

# calculate distance of each quad to its respective centroid
d_distance_bed_centroids <- left_join(d_quad_centroids_long, d_bed_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_bed - value_quad)^2) %>%
  group_by(biome, ecosystem, bed, quadrat) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_bed_centroids$dist)

# calculate distance between bed centroids and quadrat centroids
d_bed_points_2_centroid <- left_join(select(d_bed_centroids, biome, ecosystem, bed, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), select(d_quad_centroids, biome, ecosystem, bed, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2)) %>%
  filter(ecosystem %in% c('Australia', 'South_Africa'))

# calculate ecosystem centroids ####

# which ecosystems have multiple beds?
d_samp %>% select(ecosystem, bed) %>%
  distinct() %>%
  group_by(ecosystem) %>%
  tally()
# australia and south africa

# create centroids from the quadrats when there is only a single bed, and from the bed centroids when there are multiple beds
d_eco_centroids <- bind_rows(filter(d_quad_centroids, !ecosystem %in% c('Australia', 'South_Africa')),
                             filter(d_bed_centroids, ecosystem %in% c('Australia', 'South_Africa'))) %>%
  group_by(biome, ecosystem) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

# make long
d_eco_centroids_long <- d_eco_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_eco')

# calculate distance of each quadrat centroid (or bed centroid) from its respective centroid
d_distance_eco_centroids <- left_join(bind_rows(filter(d_quad_centroids_long, !ecosystem %in% c('Australia', 'South_Africa')),
                                                filter(d_bed_centroids_long, ecosystem %in% c('Australia', 'South_Africa')) %>% rename(value_quad = value_bed)), d_eco_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_eco - value_quad)^2) %>%
  group_by(biome, ecosystem, bed, quadrat) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_eco_centroids$dist)

# calculate distance between ecosystem centroids and quadrat centroids
d_eco_points_2_centroid <- left_join(select(d_eco_centroids, biome, ecosystem, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), bind_rows(filter(d_quad_centroids, !ecosystem %in% c('Australia', 'South_Africa')) %>%
                                                                                                                                              select(biome, ecosystem, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2),
                                                                                                                                            filter(d_bed_centroids, ecosystem %in% c('Australia', 'South_Africa')) %>%
                                                                                                                                              select(biome, ecosystem, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2)))

# calculate biome centroids ####

# create centroids from the ecosystems
d_biome_centroids <- d_eco_centroids %>%
  group_by(biome) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

# make long
d_biome_centroids_long <- d_biome_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_biome')

# calculate distance of each ecosystem centroid from its respective centroid
d_distance_biome_centroids <- left_join(d_eco_centroids_long, d_biome_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_biome - value_eco)^2) %>%
  group_by(biome, ecosystem) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_biome_centroids$dist)

# set colours for plot
colours <- RColorBrewer::brewer.pal(11, 'Spectral')

# make plot - pretty gnarly but think THINK it works
ggplot() +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_ext_points_2_centroid, alpha = 0.1) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_corner_points_2_centroid, alpha = 0.3) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_quad_points_2_centroid, alpha = 0.3) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_bed_points_2_centroid, alpha = 0.3) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_eco_points_2_centroid, alpha = 0.7) +
  geom_point(aes(axis_1, axis_2), d_biome_centroids, size = 10) +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_eco_centroids, size = 8) + 
  scale_colour_manual(values = colours) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), filter(d_bed_centroids, ecosystem %in% c('Australia', 'South Africa')), size = 6, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.8)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_quad_centroids, size = 5, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.6)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_corner_centroids, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.4)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_ext_centroids, size = 1.5, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.2)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), size = 0.5, d_samples, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.2)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = 'Axis 1 (5.9%)',
       y = 'Axis 2 (4.3%)') +
  #facet_wrap(~ecosystem, scales = 'free') +
  NULL

ggsave('plots/pcoa_all_fungi.png', width = 12, height = 8)

# combine the distances into one dataframe
d_distance_all <- bind_rows(
  select(d_distance_ext_centroids, biome, ecosystem, dist) %>% mutate(type = 'PCR',
                                                                      order = 1),
  select(d_distance_corner_centroids, biome, ecosystem, dist) %>% mutate(type = 'Extraction',
                                                                         order = 2),
  select(d_distance_quad_centroids, biome, ecosystem, dist) %>% mutate(type = 'Quadrat corner',
                                                                       order = 4),
  select(d_distance_bed_centroids, biome, ecosystem, dist) %>% mutate(type = 'Quadrat',
                                                                      order = 5),
  select(d_distance_eco_centroids, biome, ecosystem, dist) %>% mutate(type = 'Bed',
                                                                      order = 6),
  select(d_distance_biome_centroids, biome, ecosystem, dist) %>% mutate(type = 'Ecosystem',
                                                                        order = 7)
)

ggplot(d_distance_all, aes(forcats::fct_reorder(type, order), dist)) +
  geom_pretty_boxplot(aes(fill = type, col = type)) +
  geom_point(aes(col = type), shape = 21, fill = 'white', position = position_jitter(width = 0.3)) +
  theme_bw(base_size = 12) +
  labs(x = 'Level of replication',
       y = 'Distance between replicates at that level') +
  scale_color_brewer('Level of replication', palette = 'Set1', type = 'qual') +
  scale_fill_brewer('Level of replication', palette = 'Set1', type = 'qual') +
  facet_wrap(~biome)

ggsave('plots/distance_all_fungi.png', width = 8, height = 4)

#-----------------------------#
# third for all eukaryotes ####
#-----------------------------#

# load in bacteria filtered phyloseq object
ps_euk <- readRDS('data/sequencing/processed/phyloseq/ps_euk_sub.rds')

ps_euk

# read in sample data from the PCR sample info
sample_data <- read.csv('data/Eden_PCRs_sample_info.csv') %>%
  select(sample_id, ecosystem = Ecosystem) %>%
  filter(!is.na(ecosystem))

# check sample data
d_samp <- sample_data(ps_euk) %>%
  data.frame()

filter(sample_data,! sample_id %in% d_samp$sample_id)
# 6 samples are not present but that is ok

group_by(d_samp, sample_id) %>%
  tally() %>%
  View()

# transform counts to relative abundances for ordination
ps_prop <- transform_sample_counts(ps_euk, function(x){x / sum(x)})

# create a distance matrix for the bacteria using the relative abundances

# if 'data/sequencing/processed/phyloseq/bact_dist_bray.rds' exists, read it in, if not, run the code

if(file.exists('data/sequencing/processed/phyloseq/euk_dist_bray.rds')){
  dist_bray <- readRDS('data/sequencing/processed/phyloseq/euk_dist_bray.rds')
} else {
  dist_bray <- distance(ps_prop, method = 'bray')
  saveRDS(dist_bray, 'data/sequencing/processed/phyloseq/euk_dist_bray.rds')
}

# create a pcoa data set - choose only eigenvalues that are > 0
#https://stackoverflow.com/questions/8924488/applying-the-pvclust-r-function-to-a-precomputed-dist-object#27148408
# error message tells how may eigenvalues are > 0
d_pcoa_correct <- ape::pcoa(dist_bray, correction = 'cailliez') 

# correct eigenvalues
correct_eigenvalues <- ape::pcoa(dist_bray, correction = 'cailliez') %>%
  .$values %>%
  janitor::clean_names() %>%
  mutate(axis = 1:n()) %>%
  mutate(axis = paste('axis', axis, sep = '_'))

# how many axes to explain 50% of the data
num_axes <- filter(correct_eigenvalues, cum_corr_eig < 0.5) %>%
  nrow()
# 22

# plot scree plot for every eigenvalue worth 1%
filter(correct_eigenvalues, rel_corr_eig > 0.01) %>%
  ggplot(., aes(x = parse_number(axis), y = rel_corr_eig)) +
  geom_col() +
  labs(x = 'PCoA', y = 'Proportion of Variance Explained')

# grab out the individual points of the pcoa
d_pcoa <- d_pcoa_correct$vectors %>%
  data.frame() %>%
  janitor::clean_names()

# check rownames of d_pcoa are the same as in d_samp
all(rownames(d_pcoa) == rownames(d_samp))

# add in the sample data
d_samples <- cbind(d_samp, d_pcoa)

# make long format
d_samples_long <- d_samples %>%
  select(-sample_id) %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value')

# calculate extraction centroids based on a each extraction ####

d_ext_centroids <- d_samples %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner, extraction) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

# make long format
d_ext_centroids_long <- d_ext_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_ext') 

# calculate distance of each PCR to its respective centroid
d_distance_ext_centroids <- left_join(d_samples_long, d_ext_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value - value_ext)^2) %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner, extraction, pcr) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_ext_centroids$dist)

# create dataset for plotting lines between raw points and centroids
d_ext_points_2_centroid <- left_join(select(d_samples, biome, ecosystem, bed, quadrat, quadrat_corner, extraction, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2), select(d_ext_centroids, biome, ecosystem, bed, quadrat, quadrat_corner, extraction, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2))

# calculate the quadrat corner centroids from the extraction centroids ####

d_corner_centroids <- d_ext_centroids %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

d_corner_centroids_long <- d_corner_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_corner')

# calculate distance of each PCR to its respective centroid
d_distance_corner_centroids <- left_join(d_ext_centroids_long, d_corner_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_corner - value_ext)^2) %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner, extraction) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_corner_centroids$dist)

# calculate distance between quadrat corner centroids and extraction centroids
d_corner_points_2_centroid <- left_join(select(d_corner_centroids, biome, ecosystem, bed, quadrat, quadrat_corner, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), select(d_ext_centroids, biome, ecosystem, bed, quadrat, quadrat_corner, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2))

# calculate quadrat centroids ####

d_quad_centroids <- d_corner_centroids %>%
  group_by(biome, ecosystem, bed, quadrat) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

d_quad_centroids_long <- d_quad_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_quad')

# calculate distance of each PCR to its respective centroid
d_distance_quad_centroids <- left_join(d_corner_centroids_long, d_quad_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_quad - value_corner)^2) %>%
  group_by(biome, ecosystem, bed, quadrat, quadrat_corner) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_quad_centroids$dist)

# calculate distance between quadrat centroids and quadrat corner centroids
d_quad_points_2_centroid <- left_join(select(d_quad_centroids, biome, ecosystem, bed, quadrat, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), select(d_corner_centroids, biome, ecosystem, bed, quadrat, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2))

# calculate bed metrics ####
d_bed_centroids <- d_quad_centroids %>%
  group_by(biome, ecosystem, bed) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

d_bed_centroids_long <- d_bed_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_bed')

# calculate distance of each quad to its respective centroid
d_distance_bed_centroids <- left_join(d_quad_centroids_long, d_bed_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_bed - value_quad)^2) %>%
  group_by(biome, ecosystem, bed, quadrat) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_bed_centroids$dist)

# calculate distance between bed centroids and quadrat centroids
d_bed_points_2_centroid <- left_join(select(d_bed_centroids, biome, ecosystem, bed, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), select(d_quad_centroids, biome, ecosystem, bed, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2)) %>%
  filter(ecosystem %in% c('Australia', 'South_Africa'))

# calculate ecosystem centroids ####

# which ecosystems have multiple beds?
d_samp %>% select(ecosystem, bed) %>%
  distinct() %>%
  group_by(ecosystem) %>%
  tally()
# australia and south africa

# create centroids from the quadrats when there is only a single bed, and from the bed centroids when there are multiple beds
d_eco_centroids <- bind_rows(filter(d_quad_centroids, !ecosystem %in% c('Australia', 'South_Africa')),
                             filter(d_bed_centroids, ecosystem %in% c('Australia', 'South_Africa'))) %>%
  group_by(biome, ecosystem) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

# make long
d_eco_centroids_long <- d_eco_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_eco')

# calculate distance of each quadrat centroid (or bed centroid) from its respective centroid
d_distance_eco_centroids <- left_join(bind_rows(filter(d_quad_centroids_long, !ecosystem %in% c('Australia', 'South_Africa')),
                                                filter(d_bed_centroids_long, ecosystem %in% c('Australia', 'South_Africa')) %>% rename(value_quad = value_bed)), d_eco_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_eco - value_quad)^2) %>%
  group_by(biome, ecosystem, bed, quadrat) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_eco_centroids$dist)

# calculate distance between ecosystem centroids and quadrat centroids
d_eco_points_2_centroid <- left_join(select(d_eco_centroids, biome, ecosystem, centroid_PCoA1 = axis_1, centroid_PCoA2 = axis_2), bind_rows(filter(d_quad_centroids, !ecosystem %in% c('Australia', 'South_Africa')) %>%
                                                                                                                                              select(biome, ecosystem, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2),
                                                                                                                                            filter(d_bed_centroids, ecosystem %in% c('Australia', 'South_Africa')) %>%
                                                                                                                                              select(biome, ecosystem, sample_PCoA1 = axis_1, sample_PCoA2 = axis_2)))

# calculate biome centroids ####

# create centroids from the ecosystems
d_biome_centroids <- d_eco_centroids %>%
  group_by(biome) %>%
  summarise(across(starts_with('axis'), mean), .groups = 'drop')

# make long
d_biome_centroids_long <- d_biome_centroids %>%
  pivot_longer(cols = starts_with('axis'), names_to = 'axis', values_to = 'value_biome')

# calculate distance of each ecosystem centroid from its respective centroid
d_distance_biome_centroids <- left_join(d_eco_centroids_long, d_biome_centroids_long) %>%
  left_join(., select(correct_eigenvalues, axis, rel_corr_eig)) %>%
  filter(parse_number(axis) <= num_axes) %>%
  mutate(dist = rel_corr_eig*(value_biome - value_eco)^2) %>%
  group_by(biome, ecosystem) %>%
  # calculate distance
  summarise(dist = abs(sqrt(sum(dist))),
            .groups = 'drop')

hist(d_distance_biome_centroids$dist)

# set colours for plot
colours <- RColorBrewer::brewer.pal(11, 'Spectral')

# make plot - pretty gnarly but think THINK it works
ggplot() +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_ext_points_2_centroid, alpha = 0.1) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_corner_points_2_centroid, alpha = 0.3) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_quad_points_2_centroid, alpha = 0.3) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_bed_points_2_centroid, alpha = 0.3) +
  geom_segment(aes(x = centroid_PCoA1, y = centroid_PCoA2, yend = sample_PCoA2, xend = sample_PCoA1, col = ecosystem), show.legend = FALSE, d_eco_points_2_centroid, alpha = 0.7) +
  geom_point(aes(axis_1, axis_2), d_biome_centroids, size = 10) +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_eco_centroids, size = 8) + 
  scale_colour_manual(values = colours) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), filter(d_bed_centroids, ecosystem %in% c('Australia', 'South Africa')), size = 6, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.8)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_quad_centroids, size = 5, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.6)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_corner_centroids, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.4)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), d_ext_centroids, size = 1.5, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.2)) +
  ggnewscale::new_scale_colour() +
  geom_point(aes(axis_1, axis_2, col = ecosystem), size = 0.5, d_samples, show.legend = FALSE) +
  scale_colour_manual(values = purrr::map_vec(colours, alpha_hex, alpha = 0.2)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = 'Axis 1 (6.8%)',
       y = 'Axis 2 (4%)') +
  #facet_wrap(~ecosystem, scales = 'free') +
  NULL

ggsave('plots/pcoa_all_euk.png', width = 12, height = 8)

# combine the distances into one dataframe
d_distance_all <- bind_rows(
  select(d_distance_ext_centroids, biome, ecosystem, dist) %>% mutate(type = 'PCR',
                                                                      order = 1),
  select(d_distance_corner_centroids, biome, ecosystem, dist) %>% mutate(type = 'Extraction',
                                                                         order = 2),
  select(d_distance_quad_centroids, biome, ecosystem, dist) %>% mutate(type = 'Quadrat corner',
                                                                       order = 4),
  select(d_distance_bed_centroids, biome, ecosystem, dist) %>% mutate(type = 'Quadrat',
                                                                      order = 5),
  select(d_distance_eco_centroids, biome, ecosystem, dist) %>% mutate(type = 'Bed',
                                                                      order = 6),
  select(d_distance_biome_centroids, biome, ecosystem, dist) %>% mutate(type = 'Ecosystem',
                                                                        order = 7)
)

ggplot(d_distance_all, aes(forcats::fct_reorder(type, order), dist)) +
  geom_pretty_boxplot(aes(fill = type, col = type)) +
  geom_point(aes(col = type), shape = 21, fill = 'white', position = position_jitter(width = 0.3)) +
  theme_bw(base_size = 12) +
  labs(x = 'Level of replication',
       y = 'Distance between replicates at that level') +
  scale_color_brewer('Level of replication', palette = 'Set1', type = 'qual') +
  scale_fill_brewer('Level of replication', palette = 'Set1', type = 'qual')

ggsave('plots/distance_all_euk.png', width = 8, height = 4)

