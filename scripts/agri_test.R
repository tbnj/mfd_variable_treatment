library(tidyverse)
library(vegan)
library(ggspatial)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(patchwork)

setwd("/mfd_variable_treatment")

### Import data
## Dist
dist <- data.table::fread("data/2025-02-19_MFD_BC_distance_subset.csv") %>%
  column_to_rownames(var = "V1")

## Samples
samples <- colnames(dist)

### Metadata from all agriculture samples
## P04_6 does have agricultural samples, but are from lowland and peatland soils of which the 
## grasslands are either in rotation or with permanent crop. If included they result in formation
## of a separate cluster

## Frozen samples
meta.all.frozen <- readxl::read_excel("data/2025-04-14_mfd_db.xlsx") %>%
filter(fieldsample_barcode %in% samples,
       mfd_hab1 == "Fields",
       !is.na(mfd_hab2),
       !project_id %in% c("P04_3", "P04_5", "P04_6")) %>%
  mutate(treatment = "Frozen")

## Dried samples
meta.all.dried <- readxl::read_excel("data/2025-04-14_mfd_db.xlsx") %>%
  filter(fieldsample_barcode %in% samples,
         mfd_hab1 == "Fields",
         !is.na(mfd_hab2),
         project_id %in% c("P04_3", "P04_5")) %>%
  mutate(treatment = "Dried")

### Metadata from ariculture grid samples
## Frozen samples
meta.grid.frozen <- data.table::fread("data/2025-04-23_MFD_samples_grid_10km.tsv", na.strings = "") %>%
  filter(mfd_hab1 == "Fields",
         !is.na(mfd_hab2),
         !project_id %in% c("P04_3", "P04_5", "P04_6")) %>%
  mutate(treatment = "Frozen") %>%
  rename(latitude = lat,
         longitude = long) %>%
  mutate(across(sampling_date, ~as.character(.)))

## Dried samples
meta.grid.dried <- data.table::fread("data/2025-04-23_MFD_samples_grid_10km.tsv", na.strings = "") %>%
  filter(mfd_hab1 == "Fields",
         !is.na(mfd_hab2),
         project_id %in% c("P04_3", "P04_5")) %>%
  mutate(treatment = "Dried") %>%
  rename(latitude = lat,
         longitude = long) %>%
  mutate(across(sampling_date, ~as.character(.)))

## Non-agriculture frozen samples
meta.grid.natural <- data.table::fread("data/2025-04-23_MFD_samples_grid_10km.tsv", na.strings = "") %>%
  filter(mfd_sampletype == "Soil",
         mfd_areatype == "Natural",
         !is.na(mfd_hab2)) %>%
  mutate(treatment = "Frozen") %>%
  rename(latitude = lat,
         longitude = long) %>%
  mutate(across(sampling_date, ~as.character(.)))

## Crop summary for selection of samples
summary.crop <- meta.grid.dried %>%
  group_by(mfd_hab2) %>%
  summarise(dried_grid = n()) %>%
  left_join(meta.grid.frozen %>%
              group_by(mfd_hab2) %>%
              summarise(frozen_grid = n())) %>%
  left_join(meta.all.dried %>%
              group_by(mfd_hab2) %>%
              summarise(dried_all = n()) %>%
              left_join(meta.all.frozen %>%
                          group_by(mfd_hab2) %>%
                          summarise(frozen_all = n())))

## Find natural groups with at least 30 samples
summary.natural <- meta.grid.natural %>%
  group_by(mfd_hab1) %>%
  summarise(frozen_grid = n()) %>%
  filter(frozen_grid >= 30)

## Pull number of crop type samples
filter.crop <- summary.crop %>%
  filter(!is.na(frozen_grid)) %>%
  group_by(mfd_hab2) %>%
  mutate(n = min(frozen_grid, dried_grid)) %>%
  pull(n)

## Name crop vector
names(filter.crop) <- summary.crop %>%
  filter(!is.na(frozen_grid)) %>%
  pull(mfd_hab2)

## Pull number of natural samples
filter.natural <- summary.natural %>%
  filter(!is.na(frozen_grid)) %>%
  group_by(mfd_hab1) %>%
  mutate(n = min(frozen_grid)) %>%
  pull(n)

## Name naturla vector
names(filter.natural) <- summary.natural %>%
  filter(!is.na(frozen_grid)) %>%
  pull(mfd_hab1)


## Make functions for filtering
filter.grid.agri <- function(data, filter) {
  out <- lst()
  
  for (i in 1:length(filter)) {
    tmp <- data %>%
      filter(mfd_hab2 %in% names(filter[i])) %>%
      slice_sample(n = filter[i])
    
    out[[i]] <- tmp
  }
  
  res <- bind_rows(out)
  
  return(res)
}

filter.grid.natural <- function(data, filter) {
  out <- lst()
  
  for (i in 1:length(filter)) {
    tmp <- data %>%
      filter(mfd_hab1 %in% names(filter[i])) %>%
      slice_sample(n = 30)
    
    out[[i]] <- tmp
  }
  
  res <- bind_rows(out)
  
  return(res)
}

## Run filtering
set.seed(123)
meta.dried.select <- filter.grid.agri(meta.grid.dried, filter.crop)
meta.frozen.select <- filter.grid.agri(meta.grid.frozen, filter.crop)
meta.frozen.natural <- filter.grid.natural(meta.grid.natural, filter.natural)

## Combine metadata
meta.comb <- meta.dried.select %>%
  rbind(meta.frozen.select) %>%
  rbind(meta.frozen.natural)


### Create DK map with frozen and dried agriculture samples
meta.agri <- meta.dried.select %>%
  rbind(meta.frozen.select)


## Import world map
world <- ne_countries(scale = "large", returnclass = "sf")

## DK map with points
map <- ggplot(data = world) + 
  geom_sf(fill = "antiquewhite") + 
  geom_point(data = meta.agri, aes(x = longitude, y = latitude, fill = mfd_hab2, shape = treatment), 
             size = 4, alpha = 1) +
  theme_bw(base_size = 16) +
  scale_shape_manual(values = c(24, 21)) +
  # scale_fill_manual(values = color_vector3) +
  # geom_text(data= world_points, aes(x=X, y=Y, label=name), color = "darkblue", fontface = "bold", check_overlap = FALSE) + 
  # ggrepel::geom_text_repel(data = points, aes(x = longitude, y = latitude, label = site_name), nudge_y = 0.3, nudge_x = 0.3) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13.5), ylim = c(54.5, 58), expand = FALSE) + 
  guides(fill = "none",
         shape = "none") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme_minimal(base_size = 19) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = "bottom")

map


### Test effect
## Pull samples from combined metdata
samples.comb <- meta.comb %>%
  pull(fieldsample_barcode)

## Pull samples from agriculture
samples.agri <- meta.agri %>%
  pull(fieldsample_barcode)

## Subset dist object
dist.comb <- dist %>%
  select(all_of(sort(samples.comb))) %>%
  filter(rownames(.) %in% sort(samples.comb))

dist.agri <- dist %>%
  select(all_of(sort(samples.agri))) %>%
  filter(rownames(.) %in% sort(samples.agri))

rm(dist)
gc()

### Perform the pcoa for agriculture samples
PCOA.agri <- ape::pcoa(dist.agri)

## plot the eigenvalues and interpret 
barplot(PCOA.agri$values$Relative_eig[1:10])

## Get percentage of variance explained by the first 3 principal coordinates
sum(as.vector(PCOA.agri$value$Relative_eig)[1])
sum(as.vector(PCOA.agri$value$Relative_eig)[2])
sum(as.vector(PCOA.agri$value$Relative_eig)[3])

## PCO1 and PCO2 explains 30% of the variation in the communities
sum(as.vector(PCOA.agri$value$Relative_eig)[1:2])

## ANOSIM 
set.seed(123)
anosim.agri <- anosim(dist.agri, meta.agri$mfd_hab2, permutations = 9999)

anosim.agri

### PERMANOVA agriculture treatment
## Test dispersion
betadisp.agri.treatment <- betadisper(as.dist(dist.agri), meta.agri$treatment)
anova(betadisp.agri.treatment)

plot(betadisp.agri.treatment)

betadisp.agri.habitat <- betadisper(as.dist(dist.agri), meta.agri$mfd_hab2)
anova(betadisp.agri.habitat)

plot(betadisp.agri.habitat)

## Perform PERMANOVA
set.seed(123)
permanova.agri <- vegan::adonis2(dist.agri ~ treatment, by = "terms",
                                 data = meta.agri, permutations = 9999)

df.permanova.agri <- permanova.agri %>%
  data.frame() %>%
  rownames_to_column("Variable") %>%
  mutate(across(Variable, ~str_replace(., "mfd_hab2", "crop"))) %>%
  mutate(across(R2, ~ as.character(round(. * 100, digits = 3)))) %>%
  mutate(pval = ifelse(Pr..F. < 0.001, "<0.001", round(Pr..F., 3))) %>%
  select(Variable, R2, pval)

df.permanova.agri

## Extract the scores for plotting with ggplot
PCOA.agri.scores <- PCOA.agri$vectors %>%
  as.data.frame() %>%
  select(1:6) %>%
  rename_with(., ~str_replace(., "Axis.", "PCO")) %>%
  cbind(fieldsample_barcode = colnames(dist.agri)) %>%
  relocate(fieldsample_barcode, .before = "PCO1") %>%
  left_join(meta.agri) %>%
  group_by(project_id) %>%
  mutate(PCO1.ctr = mean(PCO1), PCO2.ctr = mean(PCO2)) %>%
  ungroup() %>%
  mutate(label = str_c("__ANOSIM:__ ", "*R* = ", round(anosim.agri$statistic, 3), ", *p* = ", anosim.agri$signif,
                       "<br>", "__PERMANOVA:__ ", "*R* = ", df.permanova.agri$R2[1], ", *p* = ", df.permanova.agri$pval[1]))

## Create object for labelling
label.agri <- PCOA.agri.scores %>%
  select(PCO1.ctr, PCO2.ctr, project_id) %>%
  distinct()

### Visualise results
## Plot of all samples in the subset - colored by crop type (MFDO2)
p.pcoa.1v2.agri <- ggplot(data = PCOA.agri.scores, aes(x = PCO1, y = PCO2)) +
  geom_segment(aes(x = PCO1.ctr, y = PCO2.ctr, xend = PCO1, yend = PCO2), 
               lwd = 1, alpha = 0.5) +
  geom_point(aes(shape = treatment, fill = mfd_hab2), size = 4, alpha = 1) +
  scale_shape_manual(values = c(24, 21)) +
  ggrepel::geom_label_repel(data = label.agri, aes(x = PCO1.ctr, y = PCO2.ctr, label = project_id),
                            max.overlaps = Inf) +
  stat_ellipse(aes(group = treatment, linetype = treatment)) +
  theme_minimal(base_size = 19) +
  # scale_fill_manual(values = mfdo1.palette) +
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 19),
        plot.margin = margin(0,0,0,0),
        strip.text = ggtext::element_markdown(hjust = 0)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 22, size = 5, alpha = 1))) +
  labs(fill = "Crop",
       shape = "Treatment",
       linetype = "Treatment") +
  xlab("PCO1 - 12.0%") +
  ylab("PCO2 - 5.7%") +
  facet_wrap(vars(label), ncol = 1)

## Render plot
p.pcoa.1v2.agri



### Perform the pcoa for natural vs agriculture
PCOA.comb <- ape::pcoa(dist.comb)

## plot the eigenvalues and interpret 
barplot(PCOA.comb$values$Relative_eig[1:10])

## Get percentage of variance explained by the first 3 principal coordinates
sum(as.vector(PCOA.comb$value$Relative_eig)[1])
sum(as.vector(PCOA.comb$value$Relative_eig)[2])
sum(as.vector(PCOA.comb$value$Relative_eig)[3])

## PCO1 and PCO2 explains 30% of the variation in the communities
sum(as.vector(PCOA.comb$value$Relative_eig)[1:2])

## ANOSIM 
set.seed(123)
anosim.comb <- anosim(dist.comb, meta.comb$mfd_hab1, permutations = 9999)

anosim.comb

### PERMANOVA natural vs agriculture
## Test dispersion
betadisp.comb.treatment <- betadisper(as.dist(dist.comb), meta.comb$treatment)
anova(betadisp.comb.treatment)

plot(betadisp.comb.treatment)

betadisp.comb.habitat <- betadisper(as.dist(dist.comb), meta.comb$mfd_hab1)
anova(betadisp.comb.habitat)

plot(betadisp.comb.habitat)

## Perform PERMANOVA
set.seed(123)
permanova.comb <- vegan::adonis2(dist.comb ~ treatment + mfd_hab1, by = "terms",
                                 data = meta.comb, permutations = 9999)

df.permanova.comb <- permanova.comb %>%
  data.frame() %>%
  rownames_to_column("Variable") %>%
  mutate(across(Variable, ~str_replace(., "mfd_hab1", "mfdo1"))) %>%
  mutate(across(R2, ~ as.character(round(. * 100, digits = 3)))) %>%
  mutate(pval = ifelse(Pr..F. < 0.001, "<0.001", round(Pr..F., 3))) %>%
  select(Variable, R2, pval)

## Extract the scores for plotting with ggplot
PCOA.comb.scores <- PCOA.comb$vectors %>%
  as.data.frame() %>%
  select(1:6) %>%
  rename_with(., ~str_replace(., "Axis.", "PCO")) %>%
  cbind(fieldsample_barcode = colnames(dist.comb)) %>%
  relocate(fieldsample_barcode, .before = "PCO1") %>%
  left_join(meta.comb) %>%
  group_by(mfd_hab1) %>%
  mutate(PCO1.ctr = mean(PCO1), PCO2.ctr = mean(PCO2)) %>%
  ungroup() %>%
  mutate(label = str_c("__ANOSIM:__ ", "*R* = ", round(anosim.comb$statistic, 3), ", *p* = ", anosim.comb$signif,
                       "<br>", "__PERMANOVA:__ ", "*R^2^* = ", df.permanova.comb$R2[1], ", *p* = ", df.permanova.comb$pval[1]))


## Create object for labelling
label.comb <- PCOA.comb.scores %>%
  select(PCO1.ctr, PCO2.ctr, mfd_hab1) %>%
  distinct()

### Visualise results
## Plot of all samples in the subset - colored by MFDO1
p.pcoa.1v2.comb <- ggplot(data = PCOA.comb.scores, aes(x = PCO1, y = PCO2)) +
  geom_segment(aes(x = PCO1.ctr, y = PCO2.ctr, xend = PCO1, yend = PCO2, color = mfd_hab1), 
               lwd = 1, alpha = 0.5) +
  geom_point(aes(shape = treatment, fill = mfd_hab1), size = 4, alpha = 1) +
  geom_point(aes(x = PCO1.ctr, y = PCO2.ctr), size = 1, alpha = 1, color = "black") +
  scale_shape_manual(values = c(24, 21)) +
  # ggrepel::geom_label_repel(data = label.comb, aes(x = PCO1.ctr, y = PCO2.ctr, label = mfd_hab1),
  #                           max.overlaps = Inf) +
  # stat_ellipse(aes(group = mfd_hab1, linetype = mfd_hab1)) +
  theme_minimal(base_size = 12) +
  # scale_fill_manual(values = mfdo1.palette) +
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 19),
        plot.margin = margin(0,0,0,0),
        strip.text = ggtext::element_markdown(hjust = 0)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 22, size = 5, alpha = 1)),
         color = "none") +
  labs(fill = "MFDO1",
       shape = "Treatment",
       linetype = "Treatment") +
  xlab("PCO1 - 49.6%") +
  ylab("PCO2 - 16.6%") +
  facet_wrap(vars(label), ncol = 1)

p.pcoa.1v2.comb


### Combine plots
ggarranged <- map + p.pcoa.1v2.agri +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggarranged

### Save plots
png(file = 'output/MFD_agri_test.png',
    width = 1900,
    height = 1000) 
ggarranged
dev.off()

pdf(file = 'output/MFD_agri_test.pdf',
    width = 19,
    height = 12)
ggarranged
dev.off()

tiff(file = 'output/MFD_agri_test.tiff',
     width = 1900,
     height = 1200)
ggarranged
dev.off()

jpeg(file = 'output/MFD_agri_test.jpeg',
     width = 1900,
     height = 1200)
ggarranged
dev.off()

ggsave("output/MFD_agri_test.svg", plot = ggarranged, width = 19, height = 12, units = "in", dpi = "retina")


save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_agri_test.RData"))

