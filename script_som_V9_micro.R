################################################################################################################################
# Environmental Associations and Distribution of Mixoplankton within Protist Communities Across Global Oceanographic Gradients #
# MetaPR2 data + Mixoplankton Database: analysis for V9-micro dataset                                                          #
# Code by Suzana G Leles, October 28 2025                                                                                         #
################################################################################################################################

rm(list = ls())

# Load required libraries
library(maps)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(SOMbrero)
library(dplyr)
library(patchwork)
library(tidyverse)
library(vegan)

# General plotting settings ####
mytheme <-   theme(panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.background = element_rect(colour = "black", fill = "white", 
                                                   size = 1.5),
                   axis.text = element_text(size = 15, colour ="black"),
                   axis.text.x = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),
                   axis.text.y = element_text(margin = unit(c(0.5,0.5,0.5,0.5), "cm")),                   
                   axis.title = element_text(size = 17, colour ="black"),
                   axis.title.y = element_text(margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                   axis.title.x = element_text(margin = unit(c(-0.2,-0.2,-0.2,-0.2), "cm")),
                   axis.ticks.length = unit(-0.25, "cm"),
                   plot.margin = unit(c(1, 1.5, 0.5, 0.5), "lines"),
                   legend.position = c(0.15, 0.75),
                   legend.title = element_blank(),
                   legend.text = element_text(size=15, colour ="black"),
                   legend.key = element_rect(fill = "white"),
                   plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                   strip.background = element_rect(color = "black", fill = "white"),
                   strip.text = element_text(size = 15, color = "black")
)

# Import data from MetaPR2
data_meta <- read.delim("metapr2_ASVs_selected_abundance_Eukaryota_2023-08-30.tsv")
# remove depth_level as this is not the correct column (we will import it from Sample Metadata)
data_meta <- data_meta %>% select(c(-depth_level, -temperature, -salinity, -depth))
nrow(data_meta)
ncol(data_meta)

# check number of unique ASVs
nlevels(factor(data_meta$asv_code))

# Import data from Mixo Database
data_asv <- read.csv("ASVs_and_Metadata_v5.csv", fileEncoding = "Latin1")
colnames(data_asv)
nrow(data_asv)
ncol(data_asv)

# check number of unique ASVs
nlevels(factor(data_asv$asv_code))

# Add a column for a variable that combines size class and mixo type
data_asv <- data_asv %>%
  unite("type_size", mixoplankton_functional_type, size.class, sep = "_", remove = FALSE)
nrow(data_asv)
ncol(data_asv)

# Get environmental data
data_env <- read.csv("metadata_v5.csv")
# checking total number of stations
nlevels(factor(data_env$file_code))
# select the depth levels we are going to work with (only sunlit ocean)
data_env <- data_env %>% filter(depth_level %in% c("euphotic", "surface"))
# checking total number of stations
nlevels(factor(data_env$file_code))
# select focused variables
data_env <- data_env %>% select(file_code, temperature, salinity, NH4, NO3, PO4, Si, Chla, season, depth_level)

# set all negative nitrate values to zero
data_env <- data_env %>%
  mutate(NO3 = if_else(NO3 < 0, 0, NO3))

# Find common column names
intersect(colnames(data_meta), colnames(data_asv))
intersect(colnames(data_env), colnames(data_meta))
intersect(colnames(data_env), colnames(data_asv))

# Merge the MDB ASV file with the MetaPR2
d <- merge(data_meta, data_asv)
# checking we are merging all ASVs present in both datasets
nlevels(factor(d$asv_code))
# checking we are merging with all stations
nlevels(factor(d$file_code))

# Merge now with the environmental data
d <- merge(d, data_env, by = "file_code")
# checking we are merging all ASVs present in both datasets
nlevels(factor(d$asv_code))
# checking we are merging with all stations present in both datasets
nlevels(factor(d$file_code))
# check depth_level
nlevels(factor(d$depth_level))

# Exclude seaweeds (categorized as "excluded" within the column trophic_strategy)
levels(factor(d$trophic_strategy))
d1 <- d[d$trophic_strategy != "excluded", ]
# checking total number of unique ASVs
nlevels(factor(d1$asv_code)) 

# Adjust the relative abundance n_reads_pct
colnames(d1)
nlevels(factor(d1$file_code)) # checking unique station ID

# Sum reads by station ID (which is unique)
summarized_data <- d1 %>%
  group_by(file_code) %>%
  summarize(sum_reads_corrected = sum(n_reads))

# Merge the summarized data back to the original dataframe
d2 <- d1 %>% left_join(summarized_data, by = "file_code")
nlevels(factor(d2$file_code)) # checking number of unique station ID remains the same
nlevels(factor(d1$asv_code)) 

# Now we calculate the corrected relative abundance (without seaweeds)
d2$n_reads_pct_corrected <- (d2$n_reads)/d2$sum_reads_corrected

# Checking the relative abundance sums up to 100% - it does!
su <- d2 %>%
  group_by(file_code) %>%
  summarize(total_n_reads_pct_corr = sum(n_reads_pct_corrected, na.rm = TRUE))

max(su$total_n_reads_pct_corr)
min(su$total_n_reads_pct_corr)

# export merged dataset for euphotic and surface (V4 and V9) by functional type
# d2d <- d2 %>% filter(plankton_functional_type %in% c("phytoplankton"))
# write_csv(d2d, "v5_merged_data_surface_euphotic_diatoms.csv")
# d2p <- d2 %>% filter(plankton_functional_type %in% c("protozooplankton"))
# write_csv(d2p, "v5_merged_data_surface_euphotic_protozooplankton.csv")
# d2par <- d2 %>% filter(plankton_functional_type %in% c("parasite"))
# write_csv(d2par, "v5_merged_data_surface_euphotic_parasites.csv")
# d2u <- d2 %>% filter(plankton_functional_type %in% c("not assessed"))
# write_csv(d2u, "v5_merged_data_surface_euphotic_unassessed.csv")
# d2cm <- d2 %>% filter(plankton_functional_type %in% c("CM"))
# write_csv(d2cm, "v5_merged_data_surface_euphotic_CM.csv")
# d2ncm <- d2 %>% filter(plankton_functional_type %in% c("NCM"))
# write_csv(d2ncm, "v5_merged_data_surface_euphotic_NCM.csv")

# Filtering by gene region
d2 <- d2 %>% filter(d2$gene_region == "V9")
nlevels(factor(d2$file_code))
nlevels(factor(d2$asv_code))

levels(factor(d2$genus))

# Exporting ASV entries and their respective species and genus
unique_combinations_asv_species_genus <- d2 %>%
  distinct(asv_code, species, genus)

write_csv(unique_combinations_asv_species_genus, "all_v5_asv_species_genus_V9region.csv")

# Exporting merged data for V9 and euphotic (without seaweeds, i.e. "excluded")
write_csv(d2, "all_v5_all_merged_data_V9_euphotic.csv")

# Working with a subset of the variables
data_eco <- data.frame(asv = d2$asv_code, 
                       relab = d2$n_reads_pct_corrected, 
                       station = d2$file_code, 
                       type = d2$mixoplankton_functional_type,
                       size = d2$size.class,
                       type_size = d2$type_size,
                       lat = d2$latitude,
                       long = d2$longitude,
                       season = d2$season,
                       date = d2$date,
                       depth = d2$depth_level,
                       sizefraction = d2$fraction_name)

nlevels(factor(data_eco$station))
nlevels(factor(data_eco$asv))

#####################################################
# Run Self-Organizing Map for community (ASVs) data # 
#####################################################

# Filter by filter size (sampling strategy)
data_eco_wide <- data_eco %>% 
  filter(sizefraction == "micro") %>%
  select(station, lat, long, asv, relab, depth, season)

nlevels(factor(data_eco_wide$station))
nlevels(factor(data_eco_wide$asv))

# Transform dataframe from long to wide format
data_eco_wide2 <- data_eco_wide %>%
  pivot_wider(names_from = asv, values_from = relab, values_fill = list(relab = 0))

data_eco_mat <- as.matrix(data_eco_wide2[,7:ncol(data_eco_wide2)])

samples <- data_eco_wide2$station

# Scale the dataset to perform the SOM
data_eco_mat <- apply(data_eco_mat, 2, scale)

rownames(data_eco_mat) <- samples

# Run the SOM
set.seed(31441)
eco_som <- trainSOM(data_eco_mat, scaling = "none", dimension = c(6, 6), nb.save = 10, maxit = 2000)

# Perform the hierarchical clustering
plot(eco_som, what = "energy")
plot(superClass(eco_som))
eco.clust <- superClass(eco_som, k = 5)
plot(eco.clust, plot.var = FALSE)

# extract clustering for mapping
ids <- eco.clust$som$clustering
length(ids)
clusters <- eco.clust$cluster
length(clusters)
som_ids <- clusters[ids]
length(som_ids)

# Add the cluster ids to our full dataset
data_eco_wide2$comclust <- som_ids
data_eco_wide2 %>% count(comclust)
# retrieve community SOM cluster from SOM analysis for RDA
comclust <- data_eco_wide2$comclust
som_ids <- data_eco_wide2$comclust

# Plot map and sampling stations based on the clustering

# Get world map data
world <- map_data("world")

p <- ggplot() +
  geom_map(data = world, map = world, aes(x = long, y = lat, map_id = region),
           colour = "darkgray", fill = "darkgray") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white", size = 0.5, colour = 'darkgray'),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.width = unit(.5, "cm"),
        legend.key.height = unit(.95, "cm"),
        strip.background = element_rect(colour="darkgray", fill="white"),
        strip.text = element_text(size = 19),
        plot.title = element_text(hjust = 0, face = "bold", size = 22)
  )

# Load the required package
library(RColorBrewer)

col_com <- c("#762A83", "#C2A5CF", "#FFEE99", "#ACD39E", "#1B7837")

pc1 <- p + geom_point(data = data_eco_wide2, aes(x = long, y = lat, fill = factor(comclust)), size = 4, shape = 21) +
  scale_fill_manual(values = col_com, name = "Community Clusters", guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(legend.key = element_blank()) +
  theme(legend.title = element_text(size = 15, colour = "black"),
        legend.position = "top")

pc1

ggsave("fig_SOM_map_V9_micro.png", pc1, device = "png", width = 6, height = 4, dpi = 300) # 8.66

pbyc <- p + geom_point(data = data_eco_wide2, aes(x = long, y = lat, fill = factor(comclust)), size = 4, shape = 21) +
  scale_fill_manual(values = col_com, name = "Community Clusters", guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  facet_wrap(~ comclust) +
  theme(legend.key = element_blank()) +
  theme(legend.title = element_text(size = 15, colour = "black"),
        legend.position = "top")

pbyc

ggsave("fig_SOM_map_by_cluster_V9_micro.png", pbyc, device = "png", width = 8, height = 5, dpi = 300) # 8.66

station_counts <- data_eco_wide2 %>%
  group_by(lat, long) %>%
  summarise(count = n_distinct(station))

station_counts
max(station_counts$count)

# Characterize the environment across the clusters
env_data <- data.frame(station = data_env$file_code,
                       temperature = data_env$temperature, 
                       salinity = data_env$salinity,
                       NO3 = data_env$NO3,
                       depth = data_env$depth,
                       season = data_env$season)

data_eco_wide_env <- data_eco_wide2 %>%
  left_join(env_data, by = "station")

max(data_eco_wide_env$salinity, na.rm = "TRUE")
min(data_eco_wide_env$salinity, na.rm = "TRUE")

pc2 <- ggplot(data_eco_wide_env, aes(x = salinity, y = temperature, fill = factor(comclust), size = NO3)) +
  geom_point(shape = 21) +
  xlim(30, 38) +
  #ylim(20, 32) +
  ylab(expression("Temperature ( " *  degree * "C)")) +
  xlab("Salinity") +
  theme_bw() +
  scale_fill_manual(values = col_com, name = "Communities", guide = "none") +
  guides(size = guide_legend(title = expression(NO[3]))) +
  mytheme +
  #scale_size_continuous(breaks = c(0.1, 1, 10, 30), range = c(1, 10),name = expression(NO[3])) +
  theme(legend.title = element_text(size=15, colour ="black"),
        legend.position = c(0.17,0.68),
        legend.key = element_blank())
pc2

pc1 + pc2

som_map_V9_micro <- data_eco_wide_env %>% select(long,lat,comclust,temperature,NO3,salinity) %>%
  mutate(dataset = "V9-micro")

write.csv(som_map_V9_micro, "som_map_V9_micro.csv")

pt <- ggplot(data_eco_wide_env, aes(x = factor(comclust), y = temperature, fill = factor(comclust))) +
  geom_violin() +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "black", color = "black") +
  ylab(expression("Temperature ( " *  degree * "C)")) +
  xlab("Community clusters") +
  scale_fill_manual(values = col_com, name = "Communities", guide = "none") +
  mytheme
pt
ps <- ggplot(data_eco_wide_env, aes(x = factor(comclust), y = salinity, fill = factor(comclust))) +
  geom_violin() +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "black", color = "black") +
  ylab("Salinity") +
  xlab("Community clusters") +
  scale_fill_manual(values = col_com, name = "Communities", guide = "none") +
  mytheme
ps
pn <- ggplot(data_eco_wide_env, aes(x = factor(comclust), y = NO3, fill = factor(comclust))) +
  geom_violin() +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "black", color = "black") +
  ylab(expression("Nitrate (" *  mu * "M)")) +
  xlab("Community clusters") +
  scale_fill_manual(values = col_com, name = "Communities", guide = "none") +
  mytheme
pn
ggplot(data_eco_wide_env, aes(x = factor(comclust), y = abs(lat), fill = factor(comclust))) +
  geom_violin() +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "black", color = "black") +
  ylab("Latitude") +
  xlab("Community clusters") +
  scale_fill_manual(values = col_com, name = "Communities", guide = "none") +
  mytheme

pt/ps/pn

ggsave("fig_V9_micro_env_V9_micro.png", pt+ps+pn, device = "png", width = 13, height = 4, dpi = 300) # 8.66

pnhist <- ggplot(data_eco_wide_env, aes(x = NO3, fill = factor(comclust))) +
  geom_histogram() +
  facet_wrap(~ comclust) +
  scale_fill_manual(values = col_com, name = "Communities", guide = "none") +
  mytheme

ggsave("fig_V9_micro_nitrate_hist_V9_micro.png", pnhist, device = "png", width = 10, height = 6, dpi = 300) # 8.66

pthist <- ggplot(data_eco_wide_env, aes(x = temperature, fill = factor(comclust))) +
  geom_histogram() +
  facet_wrap(~ comclust) +
  scale_fill_manual(values = col_com, name = "Communities", guide = "none") +
  mytheme

ggsave("fig_V9_micro_temp_hist_V9_micro.png", pthist, device = "png", width = 10, height = 6, dpi = 300) # 8.66

pshist <- ggplot(data_eco_wide_env, aes(x = salinity, fill = factor(comclust))) +
  geom_histogram() +
  facet_wrap(~ comclust) +
  scale_fill_manual(values = col_com, name = "Communities", guide = "none") +
  mytheme

ggsave("fig_V9_micro_sal_hist_V9_micro.png", pshist, device = "png", width = 10, height = 6, dpi = 300) # 8.66

#############################################################################################
# Statistical tests to evaluate if environmental conditions differ among community clusters #
#############################################################################################

# TEMPERATURE
ggplot(data_eco_wide_env, aes(x = temperature)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(x = "Temperature", y = "Count") 

anova_model <- aov(temperature ~ factor(comclust), data = data_eco_wide_env)

qqnorm(residuals(anova_model))
qqline(residuals(anova_model))
shapiro.test(residuals(anova_model))

# Given that assumptions of parametric tests were violated, let's do non-parametric tests
kruskal.test(temperature ~ factor(comclust), data = data_eco_wide_env)
# Now compare each cluster pair
library(FSA)
dunnTest(temperature ~ factor(comclust), data = data_eco_wide_env, method = "bonferroni")

# SALINITY
ggplot(data_eco_wide_env, aes(x = salinity)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(x = "Salinity", y = "Count") 

anova_model <- aov(salinity ~ factor(comclust), data = data_eco_wide_env)

qqnorm(residuals(anova_model))
qqline(residuals(anova_model))
shapiro.test(residuals(anova_model))

kruskal.test(salinity ~ factor(comclust), data = data_eco_wide_env)
dunnTest(salinity ~ factor(comclust), data = data_eco_wide_env, method = "bonferroni")

max(data_eco_wide_env$salinity, na.rm = "TRUE")

# NITRATE
ggplot(data_eco_wide_env, aes(x = log(NO3+1))) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(x = "Nitrate", y = "Count") 

anova_model <- aov(NO3 ~ factor(comclust), data = data_eco_wide_env)

qqnorm(residuals(anova_model))
qqline(residuals(anova_model))
shapiro.test(residuals(anova_model))

kruskal.test(NO3 ~ factor(comclust), data = data_eco_wide_env)
dunnTest(NO3 ~ factor(comclust), data = data_eco_wide_env, method = "bonferroni")

#################################################################
# Quantify the contribution of rare ASVs to patterns in the SOM #
#################################################################

type.som <- eco_som

# first transform back to long format
data_eco_long <- data_eco_wide2 %>%
  pivot_longer(cols = 6:(ncol(data_eco_wide2)-1),
               names_to = "asv",
               values_to = "relab")

# now add plankton functional type information
asv_type <- data.frame(asv = data_asv$asv_code, type = data_asv$mixoplankton_functional_type)

# add functional type information
data_eco_long2 <- data_eco_long %>%
  left_join(asv_type, by = "asv")

# group by ASV and count occurrences
asv_occurrence <- data_eco_long2 %>%
  group_by(asv) %>%
  summarize(occurrence = sum(relab > 0))

# total number of stations 
length(unique(data_eco_long2$station)) * 0.1 
# change the value in occurrence to set the percent threshold of occurrence across stations
asv_classification <- asv_occurrence %>%
  mutate(classification = ifelse(occurrence < 70, "rare", "not rare"))

# Merge the classification back into the original long dataframe
data_eco_long2 <- data_eco_long2 %>%
  left_join(asv_classification, by = "asv")

classification_counts <- data_eco_long2 %>%
  distinct(asv, classification) %>%  # Ensure each ASV is counted only once
  count(classification)
#classification_counts

# now let's summarize the data
data_eco_sum <- data_eco_long2 %>% group_by(comclust, station, classification, type) %>%
  summarise(relab = sum(relab)) %>% ungroup %>% group_by(comclust,classification, type) %>%
  summarise(mean_relab = mean(relab),
            sd_relab = sd(relab))

data_eco_sum$type <- factor(data_eco_sum$type, 
                            levels = c("CM", "eSNCM", "pSNCM", "GNCM", "not assessed", "parasite", "phytoplankton", "protozooplankton"),
                            labels = c("CM", "eSNCM", "pSNCM", "GNCM", "other phytoplankton", "parasites", "diatoms", "protozooplankton"))

# Generate 8 colors from the BuRd color palette
colors <- brewer.pal(8, "RdBu")

# plot rare VERSUS non-rare
pc4 <- ggplot(data_eco_sum, aes(x = factor(comclust), y = mean_relab, fill= classification)) +
  geom_bar(stat = "identity",position = "fill") + 
  facet_wrap(~type, ncol = 1) + 
  scale_fill_manual(values = c("gray", "gray30")) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Percetage of ASVs") +
  xlab("Community cluster id") +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = scales::percent_format()) +
  mytheme +
  theme(legend.position = "right")
pc4

ggsave("fig_rare_not_rare_V9_micro.png", pc4, device = "png", width = 5, height = 11, dpi = 300) # 8.66

###########################################
# Plot heatmaps with most abundant genera #
###########################################

# Plot the relative abundance by functional type and by community cluster (BY GENUS)

# make it long format
data_eco_wide3 <- data_eco_wide2 %>%
  pivot_longer(cols = 6:(ncol(data_eco_wide2)-1),
               names_to = "asv",
               values_to = "relab")

n_distinct(data_eco_wide3$asv)

# bring genus information
asv_genus <- data.frame(asv = data_asv$asv_code, genus = data_asv$genus)

n_distinct(asv_genus$asv)

# add genus information
data_eco_wide3 <- data_eco_wide3 %>%
  left_join(asv_genus, by = "asv")

# add type information
data_eco_wide3 <- data_eco_wide3 %>%
  left_join(asv_type, by = "asv")

any(is.na(data_eco_wide3$genus))

# count number of ASVs by genus and by functional type
d_asv_by_genus_by_type <- data_eco_wide3 %>%
  group_by(genus, type) %>%
  summarise(unique_asvs = n_distinct(asv), .groups = "drop")

sum(d_asv_by_genus_by_type$unique_asvs)

# summarize (summing all ASVs by genus)
data_eco_wide4 <- data_eco_wide3 %>%
  group_by(genus, station, comclust, type) %>%
  summarize(sum = sum(relab, na.rm = TRUE), .groups = "drop")

# get the average across stations by community cluster
data_eco_wide5 <- data_eco_wide4 %>%
  group_by(genus, comclust, type) %>%
  summarize(avg = mean(sum, na.rm = TRUE), .groups = "drop")

# checking zeros and non-zeros
data_eco_wide5 %>%
  summarize(
    zeros = sum(avg == 0, na.rm = TRUE),
    non_zeros = sum(avg != 0, na.rm = TRUE)
  )

# calculate the cumulative relative abundance in the dataframe
df_70 <- data_eco_wide5 %>%
  arrange(comclust, desc(avg)) %>%
  group_by(comclust) %>%
  mutate(cum_abundance = cumsum(avg)) %>%
  filter(cum_abundance <= 0.7) 

write.csv(df_70, "v5_sum_by_genus_averaged_by_station_V9micro.csv")

# check how many genus entries we have by comclust
df_70 %>%
  group_by(comclust) %>%
  summarize(genus_count = n_distinct(genus))

levels(factor(df_70$genus))

# Prepare for heatmap by reshaping data
df_heatmap <- data_eco_wide5 %>%
  filter(genus %in% df_70$genus & type %in% df_70$type)

any(is.na(df_heatmap))

# remove any text that appears after "_", including it
df_heatmap <- df_heatmap %>% mutate(genus = str_remove(genus, "_.*"))

df_heatmap <- df_heatmap %>%
  mutate(genus = str_replace_all(genus, "Phalachroma", "Phalacroma"))

# Now we have entries for all genus by comclust! 
genus_count <- df_heatmap %>%
  group_by(comclust) %>%
  summarize(genus_count = n_distinct(genus))

genus_count

# extract the total number of genus selected for the analysis
total_number_of_genus <- genus_count %>%
  summarise(avg_genus = mean(genus_count, na.rm = TRUE)) %>%
  pull(avg_genus)

total_number_of_genus

# order "genus" by "type"
df_heatmap$type <- factor(df_heatmap$type, levels = c("CM", "eSNCM", "pSNCM", "not assessed", "phytoplankton", "protozooplankton", "parasite")[7:1])
df_heatmap$genus <- factor(df_heatmap$genus, levels = unique(df_heatmap$genus)[total_number_of_genus:1])

df_type_color <- df_heatmap %>% group_by(genus, type) %>% summarise(type = type[1]) %>%
  group_by(type) %>% arrange(genus, .by_group = TRUE)

colors_type <- c("#B2182B", "#D6604D", "#F4A582", "lightgray", "#92C5DE", "#4393C3", "#2166AC")
colors_type2 <- c("#2166AC", "#4393C3", "#92C5DE", "lightgray","#F4A582","#D6604D","#B2182B")

# plot heatmap

library(ggtext)

# make sure genera name appear italicized
df_type_color <- df_type_color %>%
  mutate(y_label = if_else(str_detect(genus, "-"), paste0(genus), paste0("*", genus, "*")))

fig4 <- ggplot() +
  geom_tile(data = df_heatmap, aes(y = interaction(genus, type), x = factor(comclust), fill = avg)) +
  geom_point(data = df_type_color, aes(x = 0.425,y = interaction(genus, type),color = type), pch = 15, size = 4) +
  scale_y_discrete(labels = df_type_color$y_label) +
  scale_color_manual(values = colors_type2, name = "Functional Type",
                     labels = c("parasites", "protozooplankton", "diatoms","other phytoplankton", "pSNCMs", "eSNCMs", "CMs")) +
  xlab("Community clusters") +
  ylab("") +
  scale_fill_gradientn(
    colours = c("white", brewer.pal(n = 9, name = "YlOrBr")),
    limits = c(0, 0.142),
    breaks = c(0, 0.07, 0.14),
    name = expression("Relative abundance"),
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = 12, barheight = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5, 5.5), clip = "off") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.y = element_markdown(hjust = 1, margin = margin(r = 10)),
        axis.ticks.y = element_blank(), 
        legend.box = "vertical") +
  guides(color = "none")
fig4

ggsave("fig_heatmap_genus_V9micro.png", fig4, device = "png", width = 4.5, height = 7.5, dpi = 300) 

# Now do the same analysis but grouping by major oceanic biomes

# summarize (summing all ASVs by genus)
data_eco_wide4 <- data_eco_wide3 %>%
  group_by(genus, station, comclust, type) %>%
  summarize(sum = sum(relab, na.rm = TRUE), .groups = "drop")

# combine community clusters based on major oceanic biomes (analysed and identified in the R script "script_soms_comparison")
rules <- tibble(
  comclust = c(1, 2, 3, 4, 5),
  biome = c("temperate/subtropical", "subtropical/tropical", "temperate/subtropical", "subtropical/tropical", "subpolar/temperate")
)

# merge into the dataframe
data_eco_wide4 <- data_eco_wide4 %>%
  left_join(rules, by = c("comclust"))

# get the average across stations by biome
data_eco_wide5 <- data_eco_wide4 %>%
  group_by(genus, biome, type) %>%
  summarize(avg = mean(sum, na.rm = TRUE), .groups = "drop")

# checking zeros and non-zeros
data_eco_wide5 %>%
  summarize(
    zeros = sum(avg == 0, na.rm = TRUE),
    non_zeros = sum(avg != 0, na.rm = TRUE)
  )

# calculate the cumulative relative abundance in the dataframe
df_70 <- data_eco_wide5 %>%
  arrange(biome, desc(avg)) %>%
  group_by(biome) %>%
  mutate(cum_abundance = cumsum(avg)) %>%
  filter(cum_abundance <= 0.7) 

write.csv(df_70, "v5_sum_by_genus_averaged_by_station&biome_V9micro.csv")

# check how many genus entries we have by comclust
df_70 %>%
  group_by(biome) %>%
  summarize(genus_count = n_distinct(genus))

levels(factor(df_70$genus))

# Prepare for heatmap by reshaping data
df_heatmap <- data_eco_wide5 %>%
  filter(genus %in% df_70$genus & type %in% df_70$type)

any(is.na(df_heatmap))

# remove any text that appears after "_", including it
df_heatmap <- df_heatmap %>% mutate(genus = str_remove(genus, "_.*"))

df_heatmap <- df_heatmap %>%
  mutate(genus = str_replace_all(genus, "Phalachroma", "Phalacroma"))

# Now we have entries for all genus by biome! 
genus_count <- df_heatmap %>%
  group_by(biome) %>%
  summarize(genus_count = n_distinct(genus))

genus_count

# extract the total number of genus selected for the analysis
total_number_of_genus <- genus_count %>%
  summarise(avg_genus = mean(genus_count, na.rm = TRUE)) %>%
  pull(avg_genus)

total_number_of_genus

# order "genus" by "type"
df_heatmap$type <- factor(df_heatmap$type, levels = c("CM", "eSNCM", "pSNCM", "not assessed", "phytoplankton", "protozooplankton", "parasite")[7:1])
df_heatmap$genus <- factor(df_heatmap$genus, levels = unique(df_heatmap$genus)[total_number_of_genus:1])

df_type_color <- df_heatmap %>% group_by(genus, type) %>% summarise(type = type[1]) %>%
  group_by(type) %>% arrange(genus, .by_group = TRUE)

colors_type <- c("#B2182B", "#D6604D", "#F4A582", "lightgray", "#92C5DE", "#4393C3", "#2166AC")
colors_type2 <- c("#2166AC", "#4393C3", "#92C5DE", "lightgray","#F4A582","#D6604D","#B2182B")

# plot heatmap

library(ggtext)

# make sure genera name appear italicized
df_type_color <- df_type_color %>%
  mutate(y_label = if_else(str_detect(genus, "-"), paste0(genus), paste0("*", genus, "*")))

df_heatmap$biome <- factor(df_heatmap$biome, levels = c("subpolar/temperate", "temperate/subtropical", "subtropical/tropical"),
                           labels = c("Sp-Te", "Te-St", "St-Tr"))

fig4b <- ggplot() +
  geom_tile(data = df_heatmap, aes(y = interaction(genus, type), x = factor(biome), fill = avg)) +
  geom_point(data = df_type_color, aes(x = 0.425,y = interaction(genus, type),color = type), pch = 15, size = 4) +
  scale_y_discrete(labels = df_type_color$y_label) +
  scale_color_manual(values = colors_type2, name = "Trophic Category",
                     labels = c("parasites", "protozooplankton", "diatoms","other phytoplankton", "pSNCMs", "eSNCMs", "CMs")) +
  xlab("") +
  ylab("") +
  scale_fill_gradientn(
    colours = c("white", brewer.pal(n = 9, name = "YlOrBr")),
    limits = c(0, 0.142),
    breaks = c(0, 0.07, 0.14),
    name = expression("Relative abundance"),
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = 6, barheight = 0.5)) +
  coord_cartesian(expand = FALSE, xlim = c(0.6, 3.5), clip = "off") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.y = element_markdown(hjust = 1, margin = margin(r = 10)),
        axis.ticks.y = element_blank(), 
        legend.box = "vertical",
        axis.text.x = element_text(angle = 90)) +
  guides(color = "none")
fig4b

ggsave("fig_heatmap_genus_V9micro_by_biome.png", fig4b, device = "png", width = 2.8, height = 6.3, dpi = 300)

# Now plot the relative abundance for the different functional types

data_eco_wide4 <- data_eco_wide3 %>%
  group_by(type, station, comclust) %>%
  summarize(avg = sum(relab, na.rm = TRUE))

data_eco_wide5 <- data_eco_wide4 %>%
  group_by(type, comclust) %>%
  summarize(avg = mean(avg, na.rm = TRUE))

# Generate 8 colors from the BuRd color palette
colors <- brewer.pal(8, "RdBu")

data_eco_wide5$type <- factor(data_eco_wide5$type, 
                              levels = c("CM", "eSNCM", "pSNCM", "GNCM", "not assessed", "phytoplankton", "protozooplankton", "parasite"),
                              labels = c("CM", "eSNCM", "pSNCM", "GNCM", "other phytoplankton", "diatoms", "protozooplankton", "parasite"))

colors <- c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "lightgray", "#92C5DE", "#4393C3", "#2166AC")

pc5 <- ggplot(data_eco_wide5, aes(x = factor(comclust), y = avg, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors, name = "Trophic category") +
  ylab("Relative abundance") +
  xlab("Community clusters") +
  mytheme +
  #scale_y_continuous(expand = c(0, 0)) +  
  #scale_x_discrete(expand = c(0, 0)) +    
  theme(legend.position = "right", legend.key = element_blank(), legend.title = element_text(size = 16))
pc5

ggsave("fig_rel_ab_V9_micro.png", pc5, device = "png", width = 6.5, height = 4, dpi = 300)

#####################################################################################################################
# Plot the relationships between environmental drivers and relative abundance of trophic categories across clusters #
#####################################################################################################################

# first transform back to long format
data_eco_long <- data_eco_wide2 %>%
  pivot_longer(cols = 6:(ncol(data_eco_wide2)-1),
               names_to = "asv",
               values_to = "relab")

# now add plankton functional type information
asv_type <- data.frame(asv = data_asv$asv_code, type = data_asv$mixoplankton_functional_type)

# add functional type information
data_eco_long2 <- data_eco_long %>%
  left_join(asv_type, by = "asv")

# add environmental variables
data_eco_long2 <- data_eco_long2 %>%
  left_join(env_data, by = "station")

# group by com clust
data_eco_long3 <- data_eco_long2 %>%
  group_by(type, comclust, station, temperature, salinity, NO3) %>%
  summarize(sum = sum(relab, na.rm = TRUE), .groups = "drop")

data_eco_long3$type <- factor(data_eco_long3$type, 
                              levels = c("CM", "eSNCM", "pSNCM", "GNCM", "not assessed", "phytoplankton", "protozooplankton", "parasite"),
                              labels = c("CM", "eSNCM", "pSNCM", "GNCM", "other phytoplankton", "diatoms", "protozooplankton", "parasites"))

# Test if linear models are better than GLMs or GAMs

library(mgcv)
library(broom)

# Nest data by type
model_data <- data_eco_long3 %>%
  group_by(type) %>%
  nest()

# Fit lm and gam models
model_fits <- model_data %>%
  mutate(lm_fit = map(data, ~ lm(sum ~ temperature, data = .x)),
         glm_fit = map(data, ~ glm(sum ~ temperature, data = .x, family = binomial())),
         gam_fit = map(data, ~ gam(sum ~ s(temperature), family = betar(link = "logit"), data = .x)))

# Extract AIC values
model_comparison <- model_fits %>%
  mutate(
    lm_aic = map_dbl(lm_fit, AIC),
    gam_aic = map_dbl(gam_fit, AIC),
    glm_aic = map_dbl(glm_fit, AIC),
    best_model = case_when(
      glm_aic < gam_aic & glm_aic < lm_aic ~ "glm",
      gam_aic < lm_aic ~ "gam",
      TRUE ~ "lm"
    )
  ) %>%
  select(type, lm_aic, gam_aic, glm_aic, best_model)

# Extract GAM stats
gam_stats <- model_fits %>%
  mutate(gam_summary = map(gam_fit, summary),
         dev_expl = map_dbl(gam_summary, "dev.expl"),
         gam_tidy = map(gam_fit, ~ tidy(.x, parametric = FALSE))) %>%
  unnest(gam_tidy) %>%
  select(type, edf, statistic, p.value, dev_expl) %>%
  rename(gam_edf = edf, gam_F = statistic, gam_p = p.value)

# Extract LM stats
lm_stats <- model_fits %>%
  mutate(lm_glance = map(lm_fit, glance)) %>%
  unnest(lm_glance) %>%
  select(type, adj.r.squared, statistic, p.value, AIC) %>%
  rename(lm_adj_r2 = adj.r.squared,
         lm_F = statistic,
         lm_p = p.value)

# Extract GLM stats
glm_stats <- model_fits %>%
  mutate(glm_summary = map(glm_fit, summary),
         glm_tidy = map(glm_fit, ~ tidy(.x, conf.int = TRUE, conf.level = 0.95)),
         glm_glance = map(glm_fit, glance)) %>%
  unnest(glm_tidy) %>%
  select(type, term, estimate, std.error, statistic, p.value) %>%
  rename(glm_estimate = estimate,
         glm_std_error = std.error,
         glm_z = statistic,
         glm_p = p.value)

# Combine everything
model_summary <- model_comparison %>%
  left_join(gam_stats, by = "type") %>%
  select(type, best_model,gam_edf, gam_F, gam_p, dev_expl, lm_aic, glm_aic, gam_aic)

write_csv(model_summary, "model_comparison_summary_micro-V9_temp.csv")

# Generate newdata for prediction based on beta-regression GAM
pred_grid <- model_data %>%
  mutate(newdata = map(data, ~ {
    rng <- range(.x$temperature, na.rm = TRUE)
    tibble(temperature = seq(rng[1], rng[2], length.out = 200))
  })) %>%
  left_join(model_fits %>% select(type, gam_fit), by = "type") %>%
  mutate(
    gam_pred = map2(gam_fit, newdata, ~ {
      pr <- predict(.x, newdata = .y, type = "response", se.fit = TRUE)
      .y %>%
        mutate(
          fit   = pr$fit,
          lower = pr$fit - 1.96 * pr$se.fit,
          upper = pr$fit + 1.96 * pr$se.fit
        )
    })
  ) %>%
  select(type, gam_pred) %>%
  unnest(cols = gam_pred)

data_eco_long3 <- data_eco_long3 %>%
  filter(!(type == "GNCM" & sum > 0.002))

pgamt <- ggplot(data_eco_long3, aes(x = temperature, y = sum)) +
  geom_point(size = 2, alpha = 0.2) +
  #scale_color_manual(values = col_com, name = "Communities") +
  geom_ribbon(data = pred_grid,
              aes(x = temperature, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "grey80", alpha = 0.4) +
  geom_line(data = pred_grid,
            aes(x = temperature, y = fit),
            inherit.aes = FALSE,
            color = "black", size = 1.5) +
  facet_wrap(~ type, scales = "free_y", nrow = 1) +
  xlab(expression("Temperature (" * degree * "C)")) +
  ylab("Relative abundance") +
  mytheme +
  theme(legend.position = "none")

pgamt

ggsave("fig_gam_temp_micro_V9.png", pgamt, device = "png", width = 26, height = 3.5, dpi = 300)

# Plot now only for cluster 5
model_data <- data_eco_long3 %>% filter(comclust == 5, type %in% c("eSNCM", "protozooplankton", "diatoms")) %>%
  group_by(type) %>%
  nest()

model_fits <- model_data %>%
  mutate(gam_fit = map(data, ~ gam(sum ~ s(temperature), family = betar(link = "logit"), data = .x)))

gam_p <- model_fits %>%  # your data with gam_fit
  mutate(gam_tidy = map(gam_fit, ~ tidy(.x, parametric = FALSE))) %>%
  unnest(gam_tidy) %>%
  select(type, term, statistic, p.value) %>%
  rename(gam_F = statistic, gam_p = p.value)

# Generate newdata for prediction based on beta-regression GAM
pred_grid <- model_data %>%
  mutate(newdata = map(data, ~ {
    rng <- range(.x$temperature, na.rm = TRUE)
    tibble(temperature = seq(rng[1], rng[2], length.out = 200))
  })) %>%
  left_join(model_fits %>% select(type, gam_fit), by = "type") %>%
  mutate(
    gam_pred = map2(gam_fit, newdata, ~ {
      pr <- predict(.x, newdata = .y, type = "response", se.fit = TRUE)
      .y %>%
        mutate(
          fit   = pr$fit,
          lower = pr$fit - 1.96 * pr$se.fit,
          upper = pr$fit + 1.96 * pr$se.fit
        )
    })
  ) %>%
  select(type, gam_pred) %>%
  unnest(cols = gam_pred)

pc5_temp <- ggplot(data_eco_long3 %>% filter(comclust == 5, type %in% c("eSNCM", "protozooplankton", "diatoms")), aes(x = temperature, y = sum)) +
  geom_point(aes(col = factor(comclust)), size = 2, alpha = 0.6, color = "#C2A5CF") +
  scale_color_manual(values = col_com, name = "Communities") +
  geom_ribbon(data = pred_grid,
              aes(x = temperature, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "grey80", alpha = 0.4) +
  geom_line(data = pred_grid,
            aes(x = temperature, y = fit),
            inherit.aes = FALSE,
            color = "black") +
  facet_wrap(~ type, scales = "free_y") +
  xlab(expression("Temperature (" * degree * "C)")) +
  ylab("Relative abundance") +
  mytheme +
  theme(legend.position = "none")
pc5_temp

ggsave("fig_cluster5_micro_V9.png", pc5_temp, device = "png", width = 12, height = 3.4, dpi = 300)

# Do the same now for salinity

# removing outliers with low salinity values
data_eco_long3 <- data_eco_long3 %>%
  filter(between(salinity, 32, 38))

# Nest data by type
model_data <- data_eco_long3 %>%
  group_by(type) %>%
  nest()

# Fit lm and gam models
model_fits <- model_data %>%
  mutate(lm_fit = map(data, ~ lm(sum ~ salinity, data = .x)),
         glm_fit = map(data, ~ glm(sum ~ salinity, data = .x, family = binomial())),
         gam_fit = map(data, ~ gam(sum ~ s(salinity), family = betar(link = "logit"), data = .x)))

# Extract AIC values
model_comparison <- model_fits %>%
  mutate(
    lm_aic = map_dbl(lm_fit, AIC),
    gam_aic = map_dbl(gam_fit, AIC),
    glm_aic = map_dbl(glm_fit, AIC),
    best_model = case_when(
      glm_aic < gam_aic & glm_aic < lm_aic ~ "glm",
      gam_aic < lm_aic ~ "gam",
      TRUE ~ "lm"
    )
  ) %>%
  select(type, lm_aic, gam_aic, glm_aic, best_model)

# Extract GAM stats
gam_stats <- model_fits %>%
  mutate(gam_summary = map(gam_fit, summary),
         dev_expl = map_dbl(gam_summary, "dev.expl"),
         gam_tidy = map(gam_fit, ~ tidy(.x, parametric = FALSE))) %>%
  unnest(gam_tidy) %>%
  select(type, edf, statistic, p.value, dev_expl) %>%
  rename(gam_edf = edf, gam_F = statistic, gam_p = p.value)

# Extract LM stats
lm_stats <- model_fits %>%
  mutate(lm_glance = map(lm_fit, glance)) %>%
  unnest(lm_glance) %>%
  select(type, adj.r.squared, statistic, p.value, AIC) %>%
  rename(lm_adj_r2 = adj.r.squared,
         lm_F = statistic,
         lm_p = p.value)

# Extract GLM stats
glm_stats <- model_fits %>%
  mutate(glm_summary = map(glm_fit, summary),
         glm_tidy = map(glm_fit, ~ tidy(.x, conf.int = TRUE, conf.level = 0.95)),
         glm_glance = map(glm_fit, glance)) %>%
  unnest(glm_tidy) %>%
  select(type, term, estimate, std.error, statistic, p.value) %>%
  rename(glm_estimate = estimate,
         glm_std_error = std.error,
         glm_z = statistic,
         glm_p = p.value)

# Combine everything
model_summary <- model_comparison %>%
  left_join(gam_stats, by = "type") %>%
  select(type, best_model,gam_edf, gam_F, gam_p, dev_expl, lm_aic, glm_aic, gam_aic)

write_csv(model_summary, "model_comparison_summary_micro-V9_sal.csv")

# Generate newdata for prediction based on beta-regression GAM
pred_grid <- model_data %>%
  mutate(newdata = map(data, ~ {
    rng <- range(.x$salinity, na.rm = TRUE)
    tibble(salinity = seq(rng[1], rng[2], length.out = 200))
  })) %>%
  left_join(model_fits %>% select(type, gam_fit), by = "type") %>%
  mutate(
    gam_pred = map2(gam_fit, newdata, ~ {
      pr <- predict(.x, newdata = .y, type = "response", se.fit = TRUE)
      .y %>%
        mutate(
          fit   = pr$fit,
          lower = pr$fit - 1.96 * pr$se.fit,
          upper = pr$fit + 1.96 * pr$se.fit
        )
    })
  ) %>%
  select(type, gam_pred) %>%
  unnest(cols = gam_pred)

pgams <- ggplot(data_eco_long3, aes(x = salinity, y = sum)) +
  geom_point(size = 2, alpha = 0.2) +
  scale_color_manual(values = col_com, name = "Communities") +
  geom_ribbon(data = pred_grid,
              aes(x = salinity, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "grey80", alpha = 0.4) +
  geom_line(data = pred_grid,
            aes(x = salinity, y = fit),
            inherit.aes = FALSE,
            color = "black", size = 1.5) +
  facet_wrap(~ type, scales = "free_y", nrow = 1) +
  xlab("Salinity") +
  ylab("Relative abundance") +
  mytheme +
  theme(legend.position = "none")

pgams

ggsave("fig_gam_sal_micro_V9.png", pgams, device = "png", width = 26, height = 3.5, dpi = 300)

# Plot now only for cluster 5
model_data <- data_eco_long3 %>% filter(comclust == 5, type %in% c("eSNCM", "protozooplankton", "diatoms")) %>%
  group_by(type) %>%
  nest()

model_fits <- model_data %>%
  mutate(gam_fit = map(data, ~ gam(sum ~ s(salinity), family = betar(link = "logit"), data = .x)))

gam_p <- model_fits %>%  # your data with gam_fit
  mutate(gam_tidy = map(gam_fit, ~ tidy(.x, parametric = FALSE))) %>%
  unnest(gam_tidy) %>%
  select(type, term, statistic, p.value) %>%
  rename(gam_F = statistic, gam_p = p.value)

# Generate newdata for prediction based on beta-regression GAM
pred_grid <- model_data %>%
  mutate(newdata = map(data, ~ {
    rng <- range(.x$salinity, na.rm = TRUE)
    tibble(salinity = seq(rng[1], rng[2], length.out = 200))
  })) %>%
  left_join(model_fits %>% select(type, gam_fit), by = "type") %>%
  mutate(
    gam_pred = map2(gam_fit, newdata, ~ {
      pr <- predict(.x, newdata = .y, type = "response", se.fit = TRUE)
      .y %>%
        mutate(
          fit   = pr$fit,
          lower = pr$fit - 1.96 * pr$se.fit,
          upper = pr$fit + 1.96 * pr$se.fit
        )
    })
  ) %>%
  select(type, gam_pred) %>%
  unnest(cols = gam_pred)

pc5_sal <- ggplot(data_eco_long3 %>% filter(comclust == 5, type %in% c("eSNCM", "protozooplankton", "diatoms")), aes(x = salinity, y = sum)) +
  geom_point(aes(col = factor(comclust)), size = 2, alpha = 0.6, color = "#C2A5CF") +
  scale_color_manual(values = col_com, name = "Communities") +
  geom_ribbon(data = pred_grid,
              aes(x = salinity, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "grey80", alpha = 0.4) +
  geom_line(data = pred_grid,
            aes(x = salinity, y = fit),
            inherit.aes = FALSE,
            color = "black") +
  facet_wrap(~ type, scales = "free_y") +
  xlab("Salinity") +
  ylab("Relative abundance") +
  mytheme +
  theme(legend.position = "none")
pc5_sal

ggsave("fig_cluster5_sal_micro_V9.png", pc5_sal, device = "png", width = 12, height = 3.4, dpi = 300)

# Do the same for nitrate

# Nest data by type
model_data <- data_eco_long3 %>%
  group_by(type) %>%
  nest()

# Fit lm and gam models
model_fits <- model_data %>%
  mutate(lm_fit = map(data, ~ lm(sum ~ NO3, data = .x)),
         glm_fit = map(data, ~ glm(sum ~ NO3, data = .x, family = binomial())),
         gam_fit = map(data, ~ gam(sum ~ s(NO3), family = betar(link = "logit"), data = .x)))

# Extract AIC values
model_comparison <- model_fits %>%
  mutate(
    lm_aic = map_dbl(lm_fit, AIC),
    gam_aic = map_dbl(gam_fit, AIC),
    glm_aic = map_dbl(glm_fit, AIC),
    best_model = case_when(
      glm_aic < gam_aic & glm_aic < lm_aic ~ "glm",
      gam_aic < lm_aic ~ "gam",
      TRUE ~ "lm"
    )
  ) %>%
  select(type, lm_aic, gam_aic, glm_aic, best_model)

# Extract GAM stats
gam_stats <- model_fits %>%
  mutate(gam_summary = map(gam_fit, summary),
         dev_expl = map_dbl(gam_summary, "dev.expl"),
         gam_tidy = map(gam_fit, ~ tidy(.x, parametric = FALSE))) %>%
  unnest(gam_tidy) %>%
  select(type, edf, statistic, p.value, dev_expl) %>%
  rename(gam_edf = edf, gam_F = statistic, gam_p = p.value)

# Extract LM stats
lm_stats <- model_fits %>%
  mutate(lm_glance = map(lm_fit, glance)) %>%
  unnest(lm_glance) %>%
  select(type, adj.r.squared, statistic, p.value, AIC) %>%
  rename(lm_adj_r2 = adj.r.squared,
         lm_F = statistic,
         lm_p = p.value)

# Extract GLM stats
glm_stats <- model_fits %>%
  mutate(glm_summary = map(glm_fit, summary),
         glm_tidy = map(glm_fit, ~ tidy(.x, conf.int = TRUE, conf.level = 0.95)),
         glm_glance = map(glm_fit, glance)) %>%
  unnest(glm_tidy) %>%
  select(type, term, estimate, std.error, statistic, p.value) %>%
  rename(glm_estimate = estimate,
         glm_std_error = std.error,
         glm_z = statistic,
         glm_p = p.value)

# Combine everything
model_summary <- model_comparison %>%
  left_join(gam_stats, by = "type") %>%
  select(type, best_model,gam_edf, gam_F, gam_p, dev_expl, lm_aic, glm_aic, gam_aic)

write_csv(model_summary, "model_comparison_summary_micro-V9_nitrate.csv")

# Generate newdata for prediction based on beta-regression GAM
pred_grid <- model_data %>%
  mutate(newdata = map(data, ~ {
    rng <- range(.x$NO3, na.rm = TRUE)
    tibble(NO3 = seq(rng[1], rng[2], length.out = 200))
  })) %>%
  left_join(model_fits %>% select(type, gam_fit), by = "type") %>%
  mutate(
    gam_pred = map2(gam_fit, newdata, ~ {
      pr <- predict(.x, newdata = .y, type = "response", se.fit = TRUE)
      .y %>%
        mutate(
          fit   = pr$fit,
          lower = pr$fit - 1.96 * pr$se.fit,
          upper = pr$fit + 1.96 * pr$se.fit
        )
    })
  ) %>%
  select(type, gam_pred) %>%
  unnest(cols = gam_pred)

pgams <- ggplot(data_eco_long3, aes(x = NO3, y = sum)) +
  geom_point(size = 2, alpha = 0.2) +
  scale_color_manual(values = col_com, name = "Communities") +
  geom_ribbon(data = pred_grid,
              aes(x = NO3, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "grey80", alpha = 0.4) +
  geom_line(data = pred_grid,
            aes(x = NO3, y = fit),
            inherit.aes = FALSE,
            color = "black", size = 1.5) +
  facet_wrap(~ type, scales = "free_y", nrow = 1) +
  xlab(expression("Nitrate (" *  mu * "M)")) +
  ylab("Relative abundance") +
  mytheme +
  theme(legend.position = "none")

pgams

ggsave("fig_gam_nitrate_micro_V9.png", pgams, device = "png", width = 26, height = 3.5, dpi = 300)

# Plot now only for cluster 5
model_data <- data_eco_long3 %>% filter(comclust == 5, type %in% c("eSNCM", "protozooplankton", "diatoms")) %>%
  group_by(type) %>%
  nest()

model_fits <- model_data %>%
  mutate(gam_fit = map(data, ~ gam(sum ~ s(NO3), family = betar(link = "logit"), data = .x)))

gam_p <- model_fits %>%  # your data with gam_fit
  mutate(gam_tidy = map(gam_fit, ~ tidy(.x, parametric = FALSE))) %>%
  unnest(gam_tidy) %>%
  select(type, term, statistic, p.value) %>%
  rename(gam_F = statistic, gam_p = p.value)

# Generate newdata for prediction based on beta-regression GAM
pred_grid <- model_data %>%
  mutate(newdata = map(data, ~ {
    rng <- range(.x$NO3, na.rm = TRUE)
    tibble(NO3 = seq(rng[1], rng[2], length.out = 200))
  })) %>%
  left_join(model_fits %>% select(type, gam_fit), by = "type") %>%
  mutate(
    gam_pred = map2(gam_fit, newdata, ~ {
      pr <- predict(.x, newdata = .y, type = "response", se.fit = TRUE)
      .y %>%
        mutate(
          fit   = pr$fit,
          lower = pr$fit - 1.96 * pr$se.fit,
          upper = pr$fit + 1.96 * pr$se.fit
        )
    })
  ) %>%
  select(type, gam_pred) %>%
  unnest(cols = gam_pred)

pc5_no3 <- ggplot(data_eco_long3 %>% filter(comclust == 5, type %in% c("eSNCM", "protozooplankton", "diatoms")), aes(x = NO3, y = sum)) +
  geom_point(aes(col = factor(comclust)), size = 2, alpha = 0.6, color = "#C2A5CF") +
  scale_color_manual(values = col_com, name = "Communities") +
  geom_ribbon(data = pred_grid,
              aes(x = NO3, ymin = lower, ymax = upper),
              inherit.aes = FALSE,
              fill = "grey80", alpha = 0.4) +
  geom_line(data = pred_grid,
            aes(x = NO3, y = fit),
            inherit.aes = FALSE,
            color = "black") +
  facet_wrap(~ type, scales = "free_y") +
  xlab(expression("Nitrate (" *  mu * "M)")) +
  ylab("Relative abundance") +
  mytheme +
  theme(legend.position = "none")
pc5_no3

ggsave("fig_cluster5_nitrate_micro_V9.png", pc5_no3, device = "png", width = 12, height = 3.4, dpi = 300)


##################################################################
# plot the scaled relative abundance of trophic types by cluster #
##################################################################

md <- merge(data_eco, 
            data_eco_wide2[, c("lat", "long", "depth", "season", "comclust")], 
            by = c("lat", "long", "depth", "season"), 
            all.x = TRUE)

data_eco_wide_test <- md %>% filter(sizefraction == "micro") %>%
  group_by(type, station, comclust) %>%
  summarize(avg = sum(relab, na.rm = TRUE))

nrow(data_eco_wide_test)
sum(is.na(data_eco_wide_test$comclust))

data_eco_wide6 <- data_eco_wide5 %>%
  group_by(type) %>%
  mutate(mean_relab_clust = mean(avg, na.rm = TRUE)) %>%  
  ungroup() %>%
  mutate(relab_adjusted = avg - mean_relab_clust)

data_eco_wide6 <- data_eco_wide6 %>% mutate(dataset = "V9-micro")
write.csv(data_eco_wide6, "som_map_V9_micro_fig3.csv")

pc6 <- ggplot(data_eco_wide6, aes(x = factor(comclust), y = relab_adjusted, group = type)) +
  geom_line(size = 1) + 
  geom_point(size = 4, shape = 21, aes(fill = factor(comclust))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = col_com, name = "Functional Type") +
  facet_wrap(~type, ncol = 2, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12)), breaks = function(x) pretty(x, n = 3),) +
  mytheme +
  ylab("Scaled relative abundance") +
  xlab("Community clusters") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 12),
        strip.background = element_blank())
pc6

(pc1 + pc2) / pc6

##############################
# Running dbRDA for ASV data # dbRDA with Hellinger + Euclidean  
##############################

# Working with a subset of the variables
data_eco <- data.frame(asv = d2$asv_code, 
                       relab = d2$n_reads_pct_corrected, 
                       station = d2$file_code, 
                       type = d2$mixoplankton_functional_type,
                       size = d2$size.class,
                       type_size = d2$type_size,
                       lat = d2$latitude,
                       long = d2$longitude,
                       season = d2$season,
                       date = d2$date,
                       depth = d2$depth_level,
                       #envclust = d2$envclust,
                       sizefraction = d2$fraction_name,
                       temp = d2$temperature, 
                       salinity = d2$salinity,
                       NO3 = d2$NO3, 
                       PO4 = d2$PO4, 
                       Si = d2$Si, 
                       Chla = d2$Chla,
                       NH4 = d2$NH4)

# Summarize to get sum relab by functional type by station
tot_relab <- data_eco %>%
  group_by(station, type) %>%
  summarize(totrelab = sum(relab, na.rm = TRUE)) %>%
  pivot_wider(names_from = type, values_from = totrelab, values_fill = list(totrelab = 0))

# Hellinger transformed data for clustering and multivariate analyses
data_eco2 <- data_eco %>%
  mutate(hellinger = sqrt(relab))

# Filter only fraction_name = micro
data_eco_wide <- data_eco2 %>% 
  filter(sizefraction == "micro") %>%
  select(station, lat, long, temp, salinity, NO3, PO4, NH4, Chla, asv, hellinger)

min(data_eco_wide$NO3, na.rm = "TRUE")

# Transform dataframe from long to wide format
data_eco_wide2 <- data_eco_wide %>%
  pivot_wider(names_from = asv, values_from = hellinger, values_fill = list(hellinger = 0))  # Fill NA with 0 if needed

nrow(data_eco_wide2)
ncol(data_eco_wide2)

# som_ids has SOM community cluster and samples contain stations
som_samples <- data.frame(station = samples, comclust = som_ids)

# merge comclust to data for RDA
data_eco_wide2 <- data_eco_wide2 %>%
  left_join(som_samples, by = "station")

# Importing environmental and community data
data_rda <- data_eco_wide2

# Select only stations that have data for temp, salinity, and NO3 (n = 1767)
data_rda <- data_rda[!(is.na(data_rda$temp) & is.na(data_rda$salinity) & is.na(data_rda$NO3)), ]
data_rda <- data_rda[complete.cases(data_rda$temp, data_rda$salinity, data_rda$NO3), ]

rda_env2 <- data_rda %>%
  left_join(tot_relab, by = "station")

# re-ordering columns
rda_env2 <- rda_env2 %>%
  select(names(tot_relab), everything())

# Preparing abiotic and biotic datasets
rda_env <- rda_env2[,c(2:9,12:14)]
nrow(rda_env)
rda_com <- rda_env2[,18:(ncol(rda_env2)-1)]
nrow(rda_com)

# Standardization of the environmental data
env <- decostand(rda_env, "standardize")
env <- as.matrix(env)

com <- rda_com

# Run dbRDA
dbRDA <- capscale(com ~ env, distance = "euclidean", add = TRUE)
alias(dbRDA)
cor(env)
vif.cca(dbRDA)

# Extract eigenvalues
eigenvalues <- dbRDA$CCA$eig

# Compute % variance explained for each axis
variance_explained <- eigenvalues / sum(eigenvalues) * 100
variance_explained

plot(dbRDA) # plots ASVs in red and stations in black
plot(dbRDA, display = "sites", type = "points")

dbRDAsum <- summary(dbRDA)

# Extract site scores (rows = sites)
df1 <- as.data.frame(scores(dbRDA, display = "sites"))
# Extract species scores (rows = ASVs)
df2 <- as.data.frame(scores(dbRDA, display = "species"))
# Extract biplot scores (rows = environmental variables)
df3 <- as.data.frame(scores(dbRDA, display = "bp"))

#df1 <- data.frame(dbRDAsum$sites[,1:2])
#df2 <- data.frame(dbRDAsum$species[,1:2]) 
#df3 <- data.frame(dbRDAsum$biplot[,1:2])

df2 <- df2 %>% mutate(asv = rownames(df2))
df1$station <- data_rda$station
df1$comclust <- data_rda$comclust

write.csv(rda_env2, "rda_V9_micro.csv", row.names = FALSE)

df1_unique <- data_eco2 %>%
  group_by(asv) %>%
  summarize(type = first(type), .groups = 'drop')

df2_with_type <- df2 %>%
  left_join(df1_unique, by = "asv")

rownames(df3) <- sub("^env", "", rownames(df3))

# if running V9-micro
rownames(df3) <- c("CM", "GNCM", "eSNCM", "other phytoplankton", "parasite", "diatoms", "protozooplankton", "temperature", "salinity", "nitrate")

pc7 <- ggplot(df1, aes(x = CAP1, y = CAP2, fill = factor(comclust))) +
  geom_point(size = 2.5, shape = 21) +
  scale_fill_manual(values = col_com, name = "Communities") +
  geom_segment(data = df3, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
               color = 'black', arrow = arrow(length = unit(0.01,"npc")), size = 0.5, inherit.aes = FALSE) +
  geom_label_repel(data=df3, aes(x=CAP1, y=CAP2, label=rownames(df3)), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE, 
                   box.padding = 0.1, point.padding = 0.2, max.overlaps = Inf) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlab('db-RDA1 (27.5%)') +
  ylab('db-RDA2 (18.0%)') +
  mytheme +
  xlim(-2,3) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none", 
        legend.title = element_text(),
        legend.key = element_blank())

pc7

df1 <- df1 %>% mutate(dataset = "V9-micro")
df3 <- df3 %>% mutate(dataset = "V9-micro")
write.csv(df1, "rda_V9_micro1.csv")
write.csv(df3, "rda_V9_micro2.csv")

fig3a <- pc1 + pc2 + plot_layout(widths = c(1.0, 0.6))
fig3b <- pc6 + pc7 + plot_layout(widths = c(1.0, 1.0))
fig3 <- fig3a / fig3b + plot_layout(heights = c(0.7, 1.0)) + plot_annotation(tag_levels = 'A')
fig3 <- fig3 & theme(plot.tag = element_text(size = 24))

ggsave("fig3_V9_micro.png", fig3, device = "png", width = 10, height = 10, dpi = 600) 
