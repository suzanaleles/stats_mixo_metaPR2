#########################################################################################################
# Mapping Mixoplankton: Environmental Drivers and Global Distribution within Marine Protist Communities #
# MetaPR2 data + Mixoplankton Database: analysis for V9-micro, V9-total, and V4-total datasets          #
# Code by Suzana G Leles, July 25 2025                                                                  #
#########################################################################################################

rm(list = ls())

# Load required libraries
library(maps)
library(ggplot2)
library(SOMbrero)
library(dplyr)
library(patchwork)
library(tidyverse)
library(SOMbrero)
library(vegan)
library(ggrepel)

# Plotting settings ####
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

# Importing data used for RDA from V9-micro dataset
v9micro <- read.csv("rda_V9_micro.csv")
ncol(v9micro) - 11 # 7967 ASVs
v9micro$dataset <- rep("V9-micro", nrow(v9micro))

# Importing data used for RDA from V9-total dataset
v9total <- read.csv("rda_V9_total.csv")
ncol(v9total) - 11 # 11874 ASVs
v9total$dataset <- rep("V9-total", nrow(v9total))

# Importing data used for RDA from V4-total dataset
v4total <- read.csv("rda_V4_total.csv")
ncol(v4total) - 11 # 18669 ASVs
v4total$dataset <- rep("V4-total", nrow(v4total))

# Merging datasets
length(intersect(names(v9micro), names(v9total))) - 11 # 6259 ASVs in commmon
length(intersect(names(v9total), names(v4total))) - 11 # 7 ASVs in common
length(intersect(names(v9micro), names(v4total))) - 11 # 7 ASVs in common

data_all <- bind_rows(v9micro, v9total, v4total)

###############################################################################################
# Characterizing samples based on major oceanic biomes across datasets and community clusters #
###############################################################################################

# First, run PCA considering only abiotic variables. We want to color samples and clusters by major biomes
pca_abiotic <- data_all[,12:14]

# Standardization of the environmental data
abiotic_scaled <- decostand(pca_abiotic, "standardize")

# Run PCA
pca_res <- prcomp(abiotic_scaled, center = FALSE, scale. = FALSE)
summary(pca_res)
variance_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
variance_explained
pca_scores <- as.data.frame(pca_res$x[, 1:2])
pca_scores$comclust <- data_all$comclust
pca_scores$dataset <- data_all$dataset

load_raw <- as.data.frame(pca_res$rotation[, 1:2]*2)
rownames(load_raw) <- c("temperature", "salinity", "nitrate")
load_raw$Variable <- rownames(load_raw)
datasets <- unique(pca_scores$dataset)
loadings <- expand.grid(Variable = load_raw$Variable, dataset = datasets)
loadings <- merge(loadings, load_raw, by = "Variable")

col_com <- c("#762A83", "#C2A5CF", "#FFEE99", "#ACD39E", "#1B7837")

pcabiotic <- ggplot(data = pca_scores, aes(x = PC1, y = PC2, fill = factor(comclust))) +
  facet_wrap(~ dataset) +
  geom_point(color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = col_com) +
  xlab('PC1 (67.0%)') +
  ylab('PC2 (22.2%)')  +
  xlim(-2,4) +
  ylim(-2,2)+
  geom_segment(data = loadings, aes(x = 0, xend = PC1, y = 0, yend = PC2), 
               color = 'black', arrow = arrow(length = unit(0.05,"npc")), size = 1.5, inherit.aes = FALSE) +
  geom_label_repel(data = loadings, aes(x = PC1, y = PC2, label = Variable), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE, 
                   box.padding = 0.1, point.padding = 0.2, max.overlaps = Inf) + 
  mytheme +
  theme(legend.position = "top",
        legend.key = element_blank())

pcabiotic

pcabiotic2 <- ggplot(data = pca_scores, aes(x = PC1, y = PC2, fill = factor(dataset))) +
  geom_point(color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = col_com) +
  xlab('PC1 (67.0%)') +
  ylab('PC2 (22.2%)')  +
  xlim(-2,4) +
  ylim(-2,2)+
  geom_segment(data = load_raw, aes(x = 0, xend = PC1, y = 0, yend = PC2), 
               color = 'black', arrow = arrow(length = unit(0.05,"npc")), size = 1.5, inherit.aes = FALSE) +
  geom_label_repel(data = load_raw, aes(x = PC1, y = PC2, label = Variable), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE, 
                   box.padding = 0.1, point.padding = 0.2, max.overlaps = Inf) + 
  mytheme +
  theme(legend.position = "top",
        legend.key = element_blank())

pcabiotic2

# create a copy of our PCA results
pca_scores2 <- pca_scores

# define the oceanic regions for each community cluster from each dataset (based on the above PCA)
rules <- tibble(
  dataset  = c("V4-total", "V4-total", "V4-total", "V4-total", "V4-total",
               "V9-total", "V9-total", "V9-total", "V9-total", "V9-total",
               "V9-micro", "V9-micro", "V9-micro", "V9-micro", "V9-micro"),
  comclust = c(1, 2, 3, 4, 5,
               1, 2, 3, 4, 5,
               1, 2, 3, 4, 5),
  biome    = c("subpolar/temperate", "polar", "temperate/subtropical", "subtropical/tropical", "polar",
               "temperate/subtropical", "temperate/subtropical", "subpolar/temperate", "subtropical/tropical", "subtropical/tropical",
               "temperate/subtropical", "subtropical/tropical", "temperate/subtropical", "subtropical/tropical", "subpolar/temperate"))

# merge into the dataframe
pca_scores3 <- pca_scores2 %>%
  left_join(rules, by = c("dataset", "comclust"))

pca_scores3$biome <- factor(pca_scores3$biome, 
                            levels = c("polar", "subpolar/temperate", "temperate/subtropical", "subtropical/tropical"))


pca_scores3$dataset <- factor(pca_scores3$dataset, levels = c("V9-micro","V9-total","V4-total"))

# now plot PCA coloring by oceanic region to confirm our choices are aligning 

pcbiomes <- ggplot(data = pca_scores3, aes(x = PC1, y = PC2, fill = factor(biome))) +
  facet_wrap(~ dataset, ncol = 1) +
  geom_point(color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = col_com) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  xlab('PC1 (67.0%)') +
  ylab('PC2 (22.2%)')  +
  xlim(-2,4) +
  ylim(-2,2)+
  geom_segment(data = loadings, aes(x = 0, xend = PC1, y = 0, yend = PC2), 
               color = 'black', arrow = arrow(length = unit(0.04,"npc")), size = 1.0, inherit.aes = FALSE) +
  geom_label_repel(data = loadings, aes(x = PC1, y = PC2, label = Variable), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE, 
                   box.padding = 0.1, point.padding = 0.2, max.overlaps = Inf) + 
  mytheme +
  theme(legend.key = element_blank(), strip.text = element_text(size = 15, color = "black"),
        legend.position = "top", legend.title = element_blank(), strip.background = element_blank())

pcbiomes

pca_biomes <- ggplot(data = pca_scores3, aes(x = PC1, y = PC2, fill = factor(biome), shape = dataset)) +
  geom_point(size = 4) +
  scale_shape_manual(values = c(22,24,21)) +
  scale_fill_manual(values = col_com,
                    guide = guide_legend(override.aes = list(shape = 21))) +
  xlab('PC1 (67.0%)') +
  ylab('PC2 (22.2%)')  +
  xlim(-2,4) +
  ylim(-2,2)+
  geom_segment(data = load_raw, aes(x = 0, xend = PC1, y = 0, yend = PC2), 
               color = 'black', arrow = arrow(length = unit(0.05,"npc")), size = 1.5, inherit.aes = FALSE) +
  geom_label_repel(data = load_raw, aes(x = PC1, y = PC2, label = Variable), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE, 
                   box.padding = 0.1, point.padding = 0.2, max.overlaps = Inf) + 
  mytheme +
  theme(legend.position = "right",
        legend.key = element_blank())

pca_biomes

# Importing data for all datasets with SOM clustering to make global maps and T-S diagrams to evaluate if these support our biome classification
v9micro2 <- read.csv("som_map_V9_micro.csv")
v9total2 <- read.csv("som_map_V9_total.csv")
v4total2 <- read.csv("som_map_V4_total.csv")

data_all2 <- bind_rows(v9micro2, v9total2, v4total2)

data_all2 <- data_all2 %>%
  left_join(rules, by = c("dataset", "comclust"))

data_all2$biome <- factor(data_all2$biome, 
                         levels = c("polar", "subpolar/temperate", "temperate/subtropical", "subtropical/tropical"))

data_all2$dataset <- factor(data_all2$dataset, levels = c("V9-micro","V9-total","V4-total"))

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

pmaps <- p + geom_point(data = data_all2, aes(x = long, y = lat, fill = factor(biome)), size = 4, shape = 21) +
  facet_wrap(~dataset, ncol = 1) +
  scale_fill_manual(values = col_com) +
  theme(legend.key = element_blank(), strip.text = element_text(size = 15, color = "black"),
        legend.position = "none", legend.title = element_blank(), strip.background = element_blank())

pmaps

pcbiomes + pmaps

pts <- ggplot(data_all2, aes(x = salinity, y = temperature, fill = factor(biome), size = NO3)) +
  geom_point(shape = 21) +
  facet_wrap(~dataset, ncol = 1) +
  xlim(32.5, 38) +
  #ylim(20, 32) +
  ylab(expression("Temperature ( " *  degree * "C)")) +
  xlab("Salinity") +
  theme_bw() +
  scale_fill_manual(values = col_com, guide = "none") +
  guides(size = guide_legend(title = expression(NO[3]))) +
  mytheme +
  theme(legend.title = element_text(size=15, colour ="black"),
        legend.position = c(0.87,0.45), legend.key = element_blank(),
        strip.background = element_blank())

pts

combined <- pcbiomes + pmaps + pts + plot_layout(widths = c(1, 1.7, 1)) &   
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20)) 

combined <- pmaps + pcbiomes + pts + plot_layout(widths = c(1.7, 1, 1)) &   
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20)) 

ggsave("fig2.png", combined, device = "png", width = 14, height = 10, dpi = 600) 

# Now plot violin plots for environmental variables by dataset

# Give unique ids to community clusters to order by biome
data_all2 <- data_all2 %>%
  group_by(dataset) %>%
  mutate(comclust_lab = case_match(dataset,
                                   "V9-micro" ~ recode_factor(comclust,`5` = "C1", `1` = "C2", `3` = "C3", `2` = "C4", `4` = "C5"),
                                   "V9-total" ~ recode_factor(comclust,`3` = "C6", `1` = "C7", `2` = "C8", `4` = "C9", `5` = "C10"),
                                   "V4-total" ~ recode_factor(comclust,`2` = "C11", `5` = "C12", `1` = "C13", `3` = "C14", `4` = "C15"),
                                   .default = factor(comclust))) %>%
  ungroup()

pt <- ggplot(data_all2, aes(x = factor(comclust_lab), y = temperature, fill = factor(biome))) +
  geom_violin() +
  facet_wrap(~ dataset, scales = "free_x") +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "black", color = "black") +
  ylab(expression("Temperature ( " *  degree * "C)")) +
  xlab("") +
  scale_fill_manual(values = col_com) +
  mytheme +
  theme(legend.position = "top",
        legend.key = element_blank())

ps <- ggplot(data_all2, aes(x = factor(comclust_lab), y = salinity, fill = factor(biome))) +
  geom_violin() +
  facet_wrap(~ dataset, scales = "free_x") +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "black", color = "black") +
  ylab("Salinity") +
  xlab("") +
  scale_fill_manual(values = col_com) +
  mytheme +
  theme(legend.position = "none",
        legend.key = element_blank())

data_all2 <- data_all2 %>% mutate(NO3_log = log1p(NO3))

pn <- ggplot(data_all2, aes(x = factor(comclust_lab), y = NO3, fill = factor(biome))) +
  geom_violin() +
  facet_wrap(~ dataset, scales = "free_x") +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "black", color = "black") +
  ylab(expression("Nitrate (" *  mu * "M)")) +
  xlab("Community clusters") +
  scale_fill_manual(values = col_com) +
  mytheme +
  theme(legend.position = "none",
        legend.key = element_blank())

pt/ps/pn

combined <- pt/ps/pn + plot_layout(widths = c(1, 1, 1)) &   
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20)) 

ggsave("figS2.png", combined, device = "png", width = 9, height = 10, dpi = 600) 


############################################################################
# Evaluating changes in community composition across clusters and datasets #
############################################################################

v9micro3 <- read.csv("som_map_V9_micro_fig3.csv")
v9total3 <- read.csv("som_map_V9_total_fig3.csv")
v4total3 <- read.csv("som_map_V4_total_fig3.csv")

data_all3 <- bind_rows(v9micro3, v9total3, v4total3)

data_all3 <- data_all3 %>%
  left_join(rules, by = c("dataset", "comclust"))

data_all3$biome <- factor(data_all3$biome, 
                          levels = c("polar", "subpolar/temperate", "temperate/subtropical", "subtropical/tropical"))

data_all3$dataset <- factor(data_all3$dataset, levels = c("V9-micro","V9-total","V4-total"))

data_all3$type <- factor(data_all3$type, levels = c("CM","eSNCM","pSNCM", "GNCM", "unassessed", "diatoms", "protozooplankton", "parasite"))

data_all3 <- data_all3 %>%
  group_by(dataset) %>%
  mutate(comclust_lab = case_match(dataset,
                                   "V9-micro" ~ recode_factor(comclust,`5` = "C1", `1` = "C2", `3` = "C3", `2` = "C4", `4` = "C5"),
                                   "V9-total" ~ recode_factor(comclust,`3` = "C6", `1` = "C7", `2` = "C8", `4` = "C9", `5` = "C10"),
                                   "V4-total" ~ recode_factor(comclust,`2` = "C11", `5` = "C12", `1` = "C13", `3` = "C14", `4` = "C15"),
                                   .default = factor(comclust))) %>%
  ungroup()

pcom1 <- ggplot(data_all3 %>% filter(dataset == "V9-micro"), aes(x = comclust_lab, y = relab_adjusted, group = type)) +
  geom_line(size = 1) + 
  geom_point(size = 4, shape = 21, aes(fill = factor(biome))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c(col_com[2], col_com[3], col_com[4])) +
  #facet_wrap(~ type, ncol = 1) +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12)), breaks = function(x) pretty(x, n = 3),) +
  mytheme +
  ylab("Scaled relative abundance") +
  xlab("") +
  labs(title = "V9-micro") +
  theme(plot.title = element_text(size = 15, family = "sans", face = "plain"),
        legend.position = "none",
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        legend.key = element_blank())

pcom2 <- ggplot(data_all3 %>% filter(dataset == "V9-total"), aes(x = comclust_lab, y = relab_adjusted, group = type)) +
  geom_line(size = 1) + 
  geom_point(size = 4, shape = 21, aes(fill = factor(biome))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c(col_com[2], col_com[3], col_com[4])) +
  #facet_wrap(~ type, ncol = 1) +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12)), breaks = function(x) pretty(x, n = 3),) +
  mytheme +
  ylab("") +
  xlab("Community clusters") +
  labs(title = "V9-total") +
  theme(plot.title = element_text(size = 15, family = "sans", face = "plain"),
        legend.position = "none",
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        legend.key = element_blank())

pcom3 <- ggplot(data_all3 %>% filter(dataset == "V4-total"), aes(x = comclust_lab, y = relab_adjusted, group = type)) +
  geom_line(size = 1) + 
  geom_point(size = 4, shape = 21, aes(fill = factor(biome))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c(col_com[1], col_com[2], col_com[3], col_com[4])) +
  #facet_wrap(~ type, ncol = 1) +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12)), breaks = function(x) pretty(x, n = 4),) +
  mytheme +
  ylab("") +
  xlab("") +
  labs(title = "V4-total") +
  theme(
    plot.title = element_text(size = 15, family = "sans", face = "plain"),
    legend.position = "bottom",                     
    legend.direction = "horizontal",
    legend.justification.top = "left",           
    legend.location = "plot",                    
    legend.margin = margin(l = -400, t = 0, r = 0, b = 0), 
    legend.key = element_blank(),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 12),
    strip.background = element_blank()
  )

combined <- pcom1 + pcom2 + pcom3 + plot_layout(widths = c(1, 1, 1)) &   
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20)) 

ggsave("fig3.png", combined, device = "png", width = 9, height = 10, dpi = 600) 

###########################################################################################################
# Running RDA to help see differences across community clusters based on environment and trophic category #
###########################################################################################################

rda_V9_micro1 <- read.csv("rda_V9_micro1.csv")
rda_V9_micro2 <- read.csv("rda_V9_micro2.csv")

rda_V9_total1 <- read.csv("rda_V9_total1.csv")
rda_V9_total2 <- read.csv("rda_V9_total2.csv")

rda_V4_total1 <- read.csv("rda_V4_total1.csv")
rda_V4_total2 <- read.csv("rda_V4_total2.csv")

data_all4a <- bind_rows(rda_V9_micro1, rda_V9_total1, rda_V4_total1)
data_all4b <- bind_rows(rda_V9_micro2, rda_V9_total2, rda_V4_total2)

data_all4a <- data_all4a %>%
  left_join(rules, by = c("dataset", "comclust"))

data_all4a$biome <- factor(data_all4a$biome, 
                          levels = c("polar", "subpolar/temperate", "temperate/subtropical", "subtropical/tropical"))

prda1 <- ggplot(data_all4a %>% filter(dataset == "V9-micro"), aes(x = CAP1, y = CAP2, fill = factor(biome))) +
  geom_point(size = 2.5, shape = 21) +
  scale_fill_manual(values = c(col_com[2], col_com[3], col_com[4])) +
  geom_segment(data = rda_V9_micro2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
               color = 'black', arrow = arrow(length = unit(0.01,"npc")), size = 0.8, inherit.aes = FALSE) +
  geom_label_repel(data = rda_V9_micro2, aes(x=CAP1, y=CAP2, label = X), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE, 
                   box.padding = 0.6, point.padding = 0.2, max.overlaps = Inf) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlab('db-RDA1 (27.5%)') +
  ylab('db-RDA2 (18.0%)') +
  mytheme +
  xlim(-2,3) +
  labs(title = "V9-micro") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none", 
        legend.title = element_text(),
        legend.key = element_blank(),
        plot.title = element_text(size = 17, family = "sans", face = "plain"))

prda2 <- ggplot(data_all4a %>% filter(dataset == "V9-total"), aes(x = CAP1, y = CAP2, fill = factor(biome))) +
  geom_point(size = 2.5, shape = 21) +
  scale_fill_manual(values = c(col_com[2], col_com[3], col_com[4])) +
  geom_segment(data = rda_V9_total2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
               color = 'black', arrow = arrow(length = unit(0.01,"npc")), size = 0.8, inherit.aes = FALSE) +
  geom_label_repel(data = rda_V9_total2, aes(x=CAP1, y=CAP2, label = X), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE, 
                   box.padding = 0.4, point.padding = 0.2, max.overlaps = Inf) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlab('db-RDA1 (28.5%)') +
  ylab('db-RDA2 (19.5%)')  +
  mytheme +
  xlim(-2,3) +
  labs(title = "V9-total") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none", 
        legend.title = element_text(),
        legend.key = element_blank(),
        plot.title = element_text(size = 17, family = "sans", face = "plain"))

prda3 <- ggplot(data_all4a %>% filter(dataset == "V4-total"), aes(x = CAP1, y = CAP2, fill = factor(biome))) +
  geom_point(size = 2.5, shape = 21) +
  scale_fill_manual(values = c(col_com[1], col_com[2], col_com[3], col_com[4])) +
  geom_segment(data = rda_V4_total2, aes(x = 0, xend = CAP1, y = 0, yend = CAP2), 
               color = 'black', arrow = arrow(length = unit(0.01,"npc")), size = 0.8, inherit.aes = FALSE) +
  geom_label_repel(data = rda_V4_total2, aes(x=CAP1, y=CAP2, label = X), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE, 
                   box.padding = 0.4, point.padding = 0.2, max.overlaps = Inf) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlab('db-RDA1 (37.9%)') +
  ylab('db-RDA2 (29.3%)')  +
  mytheme +
  xlim(-2,3) +
  labs(title = "V4-total") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(size = 17, family = "sans", face = "plain"),
        legend.position = "top",                     
        legend.direction = "horizontal",
        legend.justification.top = "left",           
        legend.location = "plot",                    
        legend.margin = margin(l = -400, t = 0, r = 0, b = 0),)

combined <- prda1 + prda2 + prda3 + plot_layout(widths = c(1, 1, 1)) &   
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 22)) 

ggsave("fig3.png", combined, device = "png", width = 16, height = 6, dpi = 600) 

############################################################################################
# Assessing differences in community composition based on functional types across datasets #
############################################################################################

# Run PCA considering only biotic variables
pca_biotic <- data_all[,2:9]

# Standardization of the environmental data
biotic_scaled <- decostand(pca_biotic, "standardize")

# Run PCA
pca_res <- prcomp(biotic_scaled, center = FALSE, scale. = FALSE)
summary(pca_res)
variance_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
variance_explained
pca_scores <- as.data.frame(pca_res$x[, 1:2])
pca_scores$comclust <- data_all$comclust
pca_scores$dataset <- data_all$dataset

load_raw <- as.data.frame(pca_res$rotation[, 1:2]*2)
rownames(load_raw) <- c("CM", "GNCM", "eSNCM", "unassessed", "parasite", "diatoms", "protozooplankton", "pSNCM")
load_raw$Variable <- rownames(load_raw)
datasets <- unique(pca_scores$dataset)
loadings <- expand.grid(Variable = load_raw$Variable, dataset = datasets)
loadings <- merge(loadings, load_raw, by = "Variable")

pca_scores$dataset <- factor(pca_scores$dataset, levels = c("V9-micro", "V9-total", "V4-total"))

pcbiotic1 <- ggplot(data = pca_scores, aes(x = PC1, y = PC2, fill = factor(dataset))) +
  geom_point(color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = c("gray10", "gray50", "gray90")) +
  xlab('PC1 (27.5%)') +
  ylab('PC2 (16.5%)')  +
  xlim(-2.5,2.5) +
  ylim(-2.5,3)+
  geom_segment(data = load_raw, aes(x = 0, xend = PC1, y = 0, yend = PC2),
               color = 'deepskyblue3', arrow = arrow(length = unit(0.05,"npc")), size = 1.5, inherit.aes = FALSE) +
  geom_label_repel(data = load_raw, aes(x = PC1, y = PC2, label = Variable), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE,
                   box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf) +
  mytheme +
  theme(legend.position = "top",
        legend.key = element_blank())

pcbiotic1

# Since V4 is so different, run PCA only for V9s
datav9 <- bind_rows(v9micro, v9total)

# Run biotic PCA
pca_biotic <- datav9[,2:9]

# Standardization of the environmental data
biotic_scaled <- decostand(pca_biotic, "standardize")

# Run PCA
pca_res <- prcomp(biotic_scaled, center = FALSE, scale. = FALSE)
summary(pca_res)
variance_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
variance_explained
pca_scores <- as.data.frame(pca_res$x[, 1:2])
pca_scores$comclust <- datav9$comclust
pca_scores$dataset <- datav9$dataset

load_raw <- as.data.frame(pca_res$rotation[, 1:2]*2)
rownames(load_raw) <- c("CM", "GNCM", "eSNCM", "unassessed", "parasite", "diatoms", "protozooplankton", "pSNCM")
load_raw$Variable <- rownames(load_raw)
datasets <- unique(pca_scores$dataset)
loadings <- expand.grid(Variable = load_raw$Variable, dataset = datasets)
oadings <- merge(loadings, load_raw, by = "Variable")

pcbiotic2 <- ggplot(data = pca_scores, aes(x = PC1, y = PC2, fill = factor(dataset))) +
  geom_point(color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = c("gray10", "gray50")) +
  xlab('PC1 (28.2%)') +
  ylab('PC2 (18.4%)')  +
  xlim(-2.5,2.5) +
  ylim(-2.5,3)+
  geom_segment(data = load_raw, aes(x = 0, xend = PC1, y = 0, yend = PC2),
               color = 'deepskyblue3', arrow = arrow(length = unit(0.05,"npc")), size = 1.5, inherit.aes = FALSE) +
  geom_label_repel(data = load_raw, aes(x = PC1, y = PC2, label = Variable), fill = "white",
                   color = 'black', size = 3.5, inherit.aes = FALSE,
                   box.padding = 0.6, point.padding = 0.2, max.overlaps = Inf) +
  mytheme +
  theme(legend.position = "top",
        legend.key = element_blank())

pcbiotic2

combined <- pcbiotic1 + pcbiotic2 + plot_layout(widths = c(1, 1)) &
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20))

ggsave("figSM_PCAs.png", combined, device = "png", width = 10, height = 5.5, dpi = 600)

