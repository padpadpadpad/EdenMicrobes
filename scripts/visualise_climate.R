# Script to analyse environmental variables of Eden Project biomes

#--------------------------#
# what this script does ####
#--------------------------#

# load in Eden GPS data as it has all site info
d_site <- read.csv('data/Eden_GPS_data.csv') %>%
  janitor::clean_names() %>%
  select(ecosystem, diversity, biome) %>%
  distinct() %>%
  mutate(ecosystem = gsub(' ', '_', ecosystem))

# load in packages ####
librarian::shelf(tidyverse, patchwork)

# load in data ####

d <- read.csv('data/climate/Eden_microclimate_nov23.csv') %>%
  janitor::clean_names() %>%
  rename(ecosystem = habitat) %>%
  left_join(., d_site) %>%
  # remove any random time in the date column
  # remove seconds from time
  mutate(date = gsub(' .*', '', date),
         time = ifelse(str_count(time, ':') == 2,
                       gsub(':\\d{2}$', '', time),
                       time),
         time_date = dmy_hm(paste(date, time, sep = ' ')),
         time2 = hm(time),
         hour = hour(time2))

# ok should be able to plot these now
ggplot(d, aes(time_date, temperature_c)) +
  geom_point(aes(groups = datalogger)) +
  facet_wrap(~ecosystem, scales = 'free_x')
  
# average hourly time periods for each datalogger in each rep
d_hour <- group_by(d, biome, ecosystem, diversity, datalogger, rep, date, hour) %>%
  summarise(ave_temp = mean(temperature_c),
            sd_temp = sd(temperature_c),
            ave_moisture = mean(volumetricwatercontent),
            sd_moisture = sd(volumetricwatercontent),
            .groups = 'drop') %>%
  mutate(time_date = dmy_h(paste(date, hour, sep = ' ')))

# ok should be able to plot these now
ggplot(d_hour, aes(time_date, ave_temp)) +
  geom_point(aes(groups = datalogger)) +
  facet_wrap(~ecosystem + biome, scales = 'free_x')

# ok should be able to plot these now
ggplot(d_hour, aes(time_date, ave_moisture)) +
  geom_point(aes(groups = datalogger)) +
  facet_wrap(~ecosystem + biome, scales = 'free_x')

# create summary stats for each ecosystem per day
# max temp
# min temp
# mean temperature
# median temperature
# do the same for moisture
d_daily <- group_by(d, ecosystem, biome, diversity, date, datalogger) %>%
  summarise(max_temp = max(temperature_c, na.rm = TRUE),
            min_temp = min(temperature_c, na.rm = TRUE),
            mean_temp = mean(temperature_c, na.rm = TRUE),
            median_temp = median(temperature_c, na.rm = TRUE),
            max_moisture = max(volumetricwatercontent, na.rm = TRUE),
            min_moisture = min(volumetricwatercontent, na.rm = TRUE),
            mean_moisture = mean(volumetricwatercontent, na.rm = TRUE),
            median_moisture = median(volumetricwatercontent, na.rm = TRUE),
            .groups = 'drop')

# keep only the temperature columns
d_daily_temp <- d_daily %>%
  select(-contains('moisture'))

# pivot longer
d_daily_temp <- pivot_longer(d_daily_temp, cols = contains('temp'), names_to = 'metric', values_to = 'value') %>%
  mutate(metric = gsub('_temp', '', metric),
         ecosystem2 = gsub('_', ' ', ecosystem),
         biome2 = ifelse(biome == 'Humid', 'Rainforest', 'Mediterranean'))
  

# of the daily values, create averages for each ecosystem
# also include start date and end date and number of days measured and middle date between those times
d_summary_temp <- mutate(d_daily_temp, date2 = dmy(date)) %>%
  group_by(ecosystem, biome, ecosystem2, biome2, diversity, metric) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            start_date = min(date2),
            end_date = max(date2),
            n_days = n(),
            middle_date = min(date2) + (max(date2) - min(date2))/2,
            .groups = 'drop')

# look at average date of measurements
select(d_summary_temp, biome2, middle_date) %>%
  distinct() %>%
  group_by(biome2) %>%
  summarise(middle_date = mean(middle_date))

# make big plot of temperatures
p1 <- ggplot(d_summary_temp) +
  geom_point(aes(ecosystem2, value, colour = metric), shape = 21, fill = 'white', position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), d_daily_temp) +
  geom_linerange(aes(ecosystem2, ymin = mean - sd, ymax = mean + sd, group = metric), position = position_dodge(width = 0.6)) +
  geom_point(aes(ecosystem2, mean, fill = metric), shape = 21, position = position_dodge(width = 0.6), size = 3) +
  facet_wrap(~biome2, scales = 'free_x') +
  theme_bw(base_size = 16) +
  labs(x = 'habitat', y = 'Temperature (°C)') +
  scale_y_continuous(breaks = seq(5, 35, by = 5), limits = c(5,35)) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  scale_color_brewer('Metric', type = 'seq', palette = 'Reds', direction = -1) +
  scale_fill_brewer('Metric', type = 'seq', palette = 'Reds', direction = -1)

# keep only the moisture columns
d_daily_moist <- d_daily %>%
  select(-contains('temp'))

# pivot longer
d_daily_moist <- pivot_longer(d_daily_moist, cols = contains('moisture'), names_to = 'metric', values_to = 'value') %>%
  mutate(metric = gsub('_moisture', '', metric),
         ecosystem2 = gsub('_', ' ', ecosystem),
         biome2 = ifelse(biome == 'Humid', 'Rainforest', 'Mediterranean'))

# of the daily values, create averages for each ecosystem
# also include start date and end date and number of days measured and middle date between those times
d_summary_moist <- mutate(d_daily_moist, date2 = dmy(date)) %>%
  group_by(ecosystem, biome, ecosystem2, biome2, diversity, metric) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            start_date = min(date2),
            end_date = max(date2),
            n_days = n(),
            middle_date = min(date2) + (max(date2) - min(date2))/2,
            .groups = 'drop')

p2 <- ggplot(d_summary_moist) +
  geom_point(aes(ecosystem2, value, colour = metric), shape = 21, fill = 'white', position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), d_daily_moist) +
  geom_linerange(aes(ecosystem2, ymin = mean - sd, ymax = mean + sd, group = metric), position = position_dodge(width = 0.6)) +
  geom_point(aes(ecosystem2, mean, fill = metric), shape = 21, position = position_dodge(width = 0.6), size = 3) +
  facet_wrap(~biome2, scales = 'free_x') +
  theme_bw(base_size = 16) +
  labs(x = 'Habitat', y = 'Moisture (%)') +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  scale_color_brewer('Metric', type = 'seq', palette = 'Blues', direction = -1) +
  scale_fill_brewer('Metric', type = 'seq', palette = 'Blues', direction = -1)

p2

p1 + theme(axis.title.x = element_blank()) + p2 + plot_layout(ncol = 1)

ggsave('plots/temperature_moisture.png', width = 12, height = 8)
