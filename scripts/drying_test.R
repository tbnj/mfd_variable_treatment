library(readxl)
library(tidyverse)
library(ampvis2)
library(rstatix)
library(patchwork)

setwd("/mfd_variable_treatment")

### Load data
metadata <- read_excel("data/2024-11-04_drying_test_metadata.xlsx", sheet = 1) %>%
  # filter(!Treatment == "Fresh") %>%
  mutate(fieldsample_barcode = case_when(SampleID == "333" ~ "MFD00578",
                                         SampleID == "353" ~ "MFD02866"),
         across(Treatment, ~factor(., levels = c("Fresh", "Frozen", "25°C", "40°C", "60°C", "80°C"))),
         across(SeqID, ~factor(.)),
         across(SampleID, ~case_when(SampleID == "333" ~ "Inorganic",
                                     SampleID == "353" ~ "Organic")),
         across(SampleID, ~factor(.)))


## Filter out controls
metadata.samples <- metadata %>%
  filter(!is.na(SampleID))

## Load OTU table
otu <- data.table::fread("data/2025-03-06_MFD_V4_ASV_drying.csv")

### Create ampvis object and remove chloroplast sequences
ampvis <- amp_load(otutable = otu,
                   metadata =  metadata.samples) %>%
  amp_subset_taxa(tax_vector = c("Eukaryota", "Chloroplast", "Mitochondria"), remove = TRUE)

## Evaluate number of reads per sample
reads <- ampvis$abund %>%
  reframe(reads = colSums(.)) %>%
  mutate(SeqID = colnames(ampvis$abund), .after = reads) %>%
  left_join(metadata %>% select(SeqID, SampleName)) %>%
  arrange(reads)

## Rarefaction
ampvis.sub <- ampvis %>%
  amp_subset_samples(!is.na(SampleID)) %>%
  amp_subset_samples(minreads = 10000) 

ampvis.sub %>% 
  amp_rarecurve(facet_by = "Treatment", color = "SampleID") + 
  labs(y = "Number of observed ASVs")

## Rarefaction level
rarefy <- ampvis.sub$abund %>%
  reframe(reads = colSums(.)) %>%
  pull(reads) %>%
  min()

## Random subsampling without replacement
ampvis.sub.ra <- ampvis.sub %>%
  amp_rarefy(rarefy = rarefy)


### Community differences
p.heatmap <- ampvis.sub.ra %>% 
  amp_heatmap(group_by = "Treatment",
              facet_by = "SampleID",
              tax_aggregate = "Phylum",
              tax_show = 15,
              measure = mean,
              plot_colorscale = "sqrt",
              plot_legendbreaks = c(0.5,5,15,30,50),
              plot_values = TRUE) +
  theme_bw(base_size = 16) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p.heatmap


### Exploratory ordinations
## PCA
p.pca <- ampvis.sub.ra %>%
  amp_ordinate(filter_species = 0.01,
               distmeasure = "bray",
               type = "PCA",
               constrain = "Treatment",
               transform = "hellinger",
               sample_color_by = "Treatment",
               sample_shape_by = "SampleID",
               sample_trajectory_group = "SampleID",
               sample_trajectory = "SampleID",
               sample_colorframe = F) +
  scale_color_brewer(palette = "Dark2") + 
  theme_bw(base_size = 16) +
  guides(fill = guide_legend(nrow = 1, position = "bottom"),
         shape = guide_legend(title = "Sample"))

p.pca

## RDA
p.rda <- ampvis.sub.ra %>%
  amp_ordinate(filter_species = 0.01,
               distmeasure = "bray",
               type = "RDA",
               constrain = "Treatment",
               transform = "hellinger",
               sample_color_by = "Treatment",
               sample_shape_by = "SampleID",
               sample_trajectory_group = "SampleID",
               sample_trajectory = "SampleID",
               sample_colorframe = F) +
  scale_color_brewer(palette = "Dark2") + 
  theme_bw(base_size = 16) +
  guides(fill = guide_legend(nrow = 1, position = "bottom"),
         shape = guide_legend(title = "Sample"))

p.rda


### Evaluate if dryig has an effect on the amount of extracted DNA
ext.summary <- metadata.samples %>%
  group_by(Treatment) %>%
  summarise(median = round(median(ExtractionConc), 2),
            mean = round(mean(ExtractionConc), 2),
            sd = round(sd(ExtractionConc), 2))

## Test for normality
metadata.samples %>%
  shapiro_test(ExtractionConc)

metadata.samples %>%
  group_by(SampleID) %>%
  shapiro_test(ExtractionConc)

metadata.samples %>%
  group_by(Treatment) %>%
  shapiro_test(ExtractionConc)

## Test for homoscedasticity
metadata.samples %>%
  levene_test(ExtractionConc ~ Treatment)

metadata.samples %>%
  group_by(SampleID) %>%
  levene_test(ExtractionConc ~ Treatment)

## Non-parametric test for difference
metadata.samples %>%
  kruskal_test(ExtractionConc ~ Treatment)

metadata.samples %>%
  group_by(SampleID) %>%
  kruskal_test(ExtractionConc ~ Treatment)

## Non-parametric multiple comparisons
metadata.samples %>%
  wilcox_test(ExtractionConc ~ Treatment, alternative = "two.sided", 
              p.adjust.method = "bonferroni", detailed = T, ref.group = "Frozen") %>%
  mutate(comp = str_c(group1, "-", group2), .keep = "unused", .before = "p.adj") %>%
  select(comp, everything())



### Evaluate if dryig has an effect on alpha diversity
alpha.div <- ampvis.sub.ra %>%
  amp_alphadiv(rarefy = NULL,
               measure = c("uniqueotus","shannon","simpson"),
               richness = T) %>%
  rename(uniqueASVs = uniqueOTUs) %>%
  select(SeqID, SampleID, Treatment, uniqueASVs, Shannon, Simpson, Chao1, ACE) %>%
  pivot_longer(!c(SeqID, SampleID, Treatment), names_to = "index", values_to = "Diversity")

## Make plot of observed richness
p.alpha.div <- alpha.div %>%
  filter(index %in% c("uniqueASVs")) %>%
  ggplot(aes(x = Treatment, 
             y = Diversity, 
             fill = Treatment,
             shape = SampleID)) + 
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(21, 22), guide = "none") +
  # geom_boxplot(outlier.size = 0) +
  geom_point(position = position_jitterdodge(.75), size = 3) +
  scale_y_continuous(limits = c(0,10000)) +
  facet_grid(cols = vars(Treatment), scales = "free_x") +
  labs(x = "", 
       y = "Observed ASVs") +
  theme_bw(base_size = 16) +
  guides(fill = guide_legend(ncol = 1, position = "right", override.aes = list(shape = 21)))


## Alpha Diversity: Observed ASVs
alpha.div.test <- alpha.div %>%
  filter(index %in% c("uniqueASVs"))

## Summarise observed ASV richness
alpha.summary <- alpha.div.test %>%
  group_by(Treatment) %>%
  summarise(median = round(median(Diversity), 2),
            mean = round(mean(Diversity), 2),
            sd = round(sd(Diversity), 2))

## Test for normality
alpha.div.test %>%
  shapiro_test(Diversity)

alpha.div.test %>%
  group_by(SampleID) %>%
  shapiro_test(Diversity)

alpha.div.test %>%
  group_by(Treatment) %>%
  shapiro_test(Diversity)

## Test for homoscedasticity
alpha.div.test %>%
  levene_test(Diversity ~ SampleID)

alpha.div.test %>%
  levene_test(Diversity ~ Treatment)

alpha.div.test %>%
  group_by(SampleID) %>%
  levene_test(Diversity ~ Treatment)

## Non-parametric test for difference
alpha.div.test %>%
  kruskal_test(Diversity ~ Treatment)

alpha.div.test %>%
  group_by(SampleID) %>%
  kruskal_test(Diversity ~ Treatment)

## Non-parametric multiple comparisons
alpha.div.test %>%
  wilcox_test(Diversity ~ Treatment, alternative = "two.sided", 
              p.adjust.method = "bonferroni", detailed = T, ref.group = "Frozen") %>%
  mutate(comp = str_c(group1, "-", group2), .keep = "unused", .before = "p.adj") %>%
  select(comp, everything())



### Evaluate if dryig has an effect on beta diversity
## Hellinger-transformed BC dissimilarity
bc.dis <- ampvis.sub.ra$abund %>%
  t() %>%
  vegan::decostand(method = "hellinger") %>%
  vegan::vegdist(method = "bray")

## Transform to similarity and make long format
bc.sim.long <- bc.dis %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(across(where(is.numeric), ~1-.)) %>%
  rownames_to_column(var = "SeqID") %>%
  pivot_longer(!SeqID, names_to = "comp", values_to = "Diversity") %>%
  filter(!SeqID == comp) %>%
  left_join(metadata.samples %>% 
              select(SeqID, SampleID, Treatment)) %>%
  rename(sample_id = SampleID,
         sample_Treatment = Treatment) %>%
  left_join(metadata.samples %>% select(SeqID, SampleID, Treatment), by = c("comp" = "SeqID")) %>%
  rename(comp_id = SampleID,
         comp_Treatment = Treatment) %>%
  mutate(sample = case_when(sample_id == comp_id ~ "Within",
                            TRUE ~ "Between"),
         method = case_when(sample_Treatment == comp_Treatment ~ "Within",
                            TRUE ~ "Between"))

## Subset to within-sample comparisons
beta.div <- bc.sim.long %>%
  filter(sample == "Within",
         method == "Within") %>%
  mutate(facet = case_when(sample_id == "Inorganic" ~ "Inorganic",
                           TRUE ~ "Organic")) %>%
  rbind(., bc.sim.long %>%
          filter(sample == "Between",
                 method == "Within") %>%
          mutate(facet = "Between samples")) %>%
  mutate(across(facet, ~factor(., levels = c("Inorganic", "Organic", "Between samples")))) %>%
  mutate(tmp1 = str_remove(SeqID, "MQ210406-"),
         tmp2 = str_remove(comp, "MQ210406-"),
         across(tmp1:tmp2, ~as.numeric(.)),
         tmp3 = case_when(tmp1 < tmp2 ~ tmp1,
                          tmp1 > tmp2 ~ tmp2),
         tmp4 = case_when(tmp1 > tmp2 ~ tmp1,
                          tmp1 < tmp2 ~ tmp2),
         across(SeqID, ~str_c("MQ210406-", tmp3)),
         across(comp, ~str_c("MQ210406-", tmp4))) %>%
  select(-c(tmp1:tmp4)) %>%
  distinct()

## Make plot of BC similarity
p.beta.div <- beta.div %>%
  ggplot(aes(x = sample_Treatment, 
             y = Diversity, 
             fill = comp_Treatment, 
             shape = facet)) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_brewer(palette = "Dark2") +
  geom_point(position = position_jitterdodge(0.75), size = 3) +
  theme_bw(base_size = 16) +
  scale_y_continuous(limits = c(0,1)) +
  facet_grid(cols = vars(sample_Treatment), scales = "free_x") +
  labs(x = "",
       y = "BC similarity") +
  guides(fill = guide_legend(ncol = 1, position = "right", title = "Treatment",
                             override.aes = list(shape = 21)),
         shape = guide_legend(ncol = 1, position = "right", title = "Sample"))

## Summarise BC similarity
bc.sim.long.sum <- bc.sim.long %>%
  group_by(sample_Treatment, comp_Treatment) %>%
  summarise(median = round(median(Diversity), 2),
            mean = round(mean(Diversity), 2),
            sd = round(sd(Diversity), 2))

## Make comparison heatmap
p.beta.comp <- bc.sim.long.sum %>%
  ggplot() +
  geom_tile(aes(x = sample_Treatment, y = comp_Treatment, fill = mean)) +
  geom_text(aes(x = sample_Treatment, y = comp_Treatment, label = mean)) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(3, "PuOr")), trans = "sqrt") +
  theme_bw(base_size = 16) +
  guides(fill = guide_colourbar(position = "right", title = "BC similarity")) +
  labs(x = "",
       y = "")

## Summarise observed ASV richness
beta.div.test <- beta.div %>%
  filter(sample == "Within",
         method == "Within")

beta.summary <- beta.div.test %>%
  group_by(sample_Treatment) %>%
  summarise(median = round(median(Diversity), 3),
            mean = round(mean(Diversity), 3),
            sd = round(sd(Diversity), 3))

## Test for normality
beta.div.test %>%
  shapiro_test(Diversity)

beta.div.test %>%
  group_by(sample_id) %>%
  shapiro_test(Diversity)

beta.div.test %>%
  group_by(sample_Treatment) %>%
  shapiro_test(Diversity)

## Test for homoscedasticity
beta.div.test %>%
  levene_test(Diversity ~ sample_id)

beta.div.test %>%
  levene_test(Diversity ~ sample_Treatment)

beta.div.test %>%
  group_by(sample_id) %>%
  levene_test(Diversity ~ sample_Treatment)

## Non-parametric test for difference
beta.div.test %>%
  kruskal_test(Diversity ~ sample_Treatment)

beta.div.test %>%
  group_by(sample_id) %>%
  kruskal_test(Diversity ~ sample_Treatment)

## Non-parametric multiple comparisons
beta.div.test %>%
  wilcox_test(Diversity ~ sample_Treatment, alternative = "two.sided", 
              p.adjust.method = "bonferroni", detailed = T, ref.group = "Frozen") %>%
  mutate(comp = str_c(group1, "-", group2), .keep = "unused", .before = "p.adj") %>%
  select(comp, everything())


### Pairwise contrasts analysis 
meta.sub <- ampvis.sub$metadata %>%
  select(SampleID, Treatment) %>%
  mutate(across(Treatment, ~factor(., levels = rev(c("Frozen", "Fresh", "25°C", "40°C", "60°C", "80°C")))))

### Perform overall PERMANOVA
set.seed(123)
permanova.all <- vegan::adonis2(bc.dis ~ SampleID + Treatment, by = "terms", 
                                data = meta.sub, permutations = 9999)

permanova.all %>%
  data.frame() %>%
  rownames_to_column("Variable") %>%
  mutate(across(R2, ~ as.character(round(. * 100, digits = 3)))) %>%
  mutate(pval = ifelse(Pr..F. < 0.001, "<0.001", round(Pr..F., 3))) %>%
  select(Variable, R2, pval)

## Create contrast objects to partiotion of variance
C(meta.sub$Treatment, contr.sum)

treatment.mm <- meta.sub %>%
  select(-Treatment) %>%
  cbind(as.data.frame(model.matrix( ~ C(meta.sub$Treatment, contr.sum)))) %>%
  select(-"(Intercept)") %>%
  rename("Frozen-80°C" = 2,
         "Frozen-60°C" = 3,
         "Frozen-40°C" = 4,
         "Frozen-25°C" = 5,
         "Frozen-Fresh" = 6) %>%
  select(1,6,5,4,3,2)

### Perform contrasts PERMANOVA
set.seed(123)
permanova.contrast <- vegan::adonis2(bc.dis~., data = treatment.mm, by = "terms", permutations = 9999)

permanova.contrast %>%
  data.frame() %>%
  rownames_to_column("Variable") %>%
  mutate(across(R2, ~ as.character(round(. * 100, digits = 3)))) %>%
  mutate(pval = ifelse(Pr..F. < 0.001, "<0.001", round(Pr..F., 3))) %>%
  select(Variable, R2, pval)


### Arrange plots
p1 <- p.heatmap + p.beta.comp + 
  plot_layout(guides = "keep", ncol = 1) &
  theme(legend.position = "right")

p2 <- p.alpha.div + p.beta.div + 
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "right",
        legend.box = "vertical")

ggarranged <- p1 | p2

### Save plots
png(file = 'output/MFD_drying_test.png',
    width = 1900,
    height = 1000) 
ggarranged
dev.off()

pdf(file = 'output/MFD_drying_test.pdf',
    width = 19,
    height = 12)
ggarranged
dev.off()

tiff(file = 'output/MFD_drying_test.tiff',
     width = 1900,
     height = 1200)
ggarranged
dev.off()

jpeg(file = 'output/MFD_drying_test.jpeg',
     width = 1900,
     height = 1200)
ggarranged
dev.off()

ggsave("output/MFD_drying_test.svg", plot = ggarranged, width = 19, height = 12, units = "in", dpi = "retina")


save.image(paste0('output/', format(Sys.time(), "%Y-%m-%d"), "_drying_test.RData"))

