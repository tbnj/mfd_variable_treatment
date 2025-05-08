library(tidyverse)

setwd("/mfd_variable_treatment")

### Import metadata
## Lab metadata
lab.meta <- read_csv("data/2025-04-24_corrected_combined_metadata.csv") %>%
  select(fieldsample_barcode, BioSample, extraction_method) %>%
  distinct()
  
## Seq metadata
seq.meta <- read_csv("data/2023-10-11_samples_minimal_metadata_collapsed.csv")

## MFD sample metadata
mfd.db <- readxl::read_excel("data/2025-04-14_mfd_db.xlsx") %>%
  select(fieldsample_barcode, BioSample, project_id, sampling_date, longitude, latitude, starts_with("mfd_"))

## Combine metadata
comb.meta <- seq.meta %>%
  left_join(lab.meta) %>%
  left_join(mfd.db) %>%
  distinct()

## Extract sample IDs and sequencing names
samples <- comb.meta %>%
  filter(!is.na(project_id)) %>%
  select(fieldsample_barcode) %>%
  distinct()

names <- comb.meta %>%
  filter(!is.na(project_id)) %>%
  select(flat_name) %>%
  distinct()

## Write list of linked samples
data.table::fwrite(samples, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_linked_samples.csv"))

## Filter combined metadata to linked samples
comb.dedup <- comb.meta %>%
  filter(fieldsample_barcode %in% samples$fieldsample_barcode)


### Summarise data based on project, ontology and extraction method
comb.sum <- comb.meta %>%
  select(flat_name, project_id, mfd_sampletype, extraction_method) %>%
  filter(!is.na(mfd_sampletype)) %>%
  distinct() %>%
  group_by(project_id, mfd_sampletype, extraction_method) %>%
  reframe(n_sum = n()) %>%
  mutate(sampling_method = 
           case_when(project_id %in% c("P06_1") ~ "AU-BIOMAP",
                     project_id %in% c("P01_1", "P01_2", "P02_1", "P02_2") ~ "AU-BIOWIDE",
                     project_id %in% c("P12_1", "P12_2") ~ "AU-SUB-SEDIMENT",
                     project_id %in% c("P03_1") ~ "KU-FORENSIC",
                     project_id %in% c("P09_2") ~ "AU-PONDERFULL",
                     project_id %in% c("P11_2") ~ "AU-SEDIMENT",
                     project_id %in% c("P12_5") ~ "GOV-SEDIMENT",
                     project_id %in% c("P06_2") ~ "MFD-DEEP",
                     project_id %in% c("P04_3", "P04_5") ~ "MFD-DRY",
                     project_id %in% c("P18_1", "P18_2") ~ "MFD-HARBOUR",
                     project_id %in% c("P05_1", "P05_2", "P09_1",
                                       "P10_1", "P10_2", "P10_3", "P11_1") ~ "MFD-SEDIMENT",
                     project_id %in% c("P11_3") & mfd_sampletype == "Water" ~ "MFD-WATER",
                     project_id %in% c("P11_3") & mfd_sampletype == "Sediment" ~ "MFD-SEDIEMENT",
                     project_id %in% c("P04_2", "P04_4", "P04_6", "P04_7", "P08_1",
                                       "P08_2", "P08_3", "P08_6",
                                       "P08_7", "P08_8", "P17_1") ~ "MFD-SOIL",
                     project_id %in% c("P08_5") & mfd_sampletype == "Soil" ~ "MFD-SOIL",
                     project_id %in% c("P08_5") & mfd_sampletype == "Sediment" ~ "MFD-SEDIMENT",
                     project_id %in% c("P04_8") ~ "MFD-SOIL*",
                     project_id %in% c("P09_3", "P09_4") ~ "MFD-SEDIMENT",
                     project_id %in% c("P13_3") ~ "AAU-AD",
                     project_id %in% c("P13_2") ~ "AAU-AS",
                     project_id %in% c("P16_1") ~ "AAU-DRINKING",
                     project_id %in% c("P16_2") ~ "AAU-SANDFILTER",
                     project_id %in% c("P16_3") & extraction_method == "PowerSoil-Pro-HT" ~ "DTU-SANDFILTER",
                     project_id %in% c("P16_3") & extraction_method == "NucliSens-miniMAG" ~ "DTU-DRINKINGWATER",
                     project_id %in% c("P16_3") & extraction_method == "PowerWater" ~ "DTU-DRINKINGWATER",
                     project_id %in% c("P13_1") ~ "AAU-INFLUENT",
                     project_id %in% c("P16_4") ~ "MFD-DEEP",
                     project_id %in% c("P06_3") ~ "MFD-DEEP",
                     project_id %in% c("P25_1") ~ "AAU-SCRAPEOFF",
                     project_id %in% c("P21_1") ~ "AU-SEEP",
                     project_id %in% c("P19_1") ~ "AAU-ARCHAELOGICAL",
                     project_id %in% c("P20_1") & mfd_sampletype == "Soil" ~ "MFD-SOIL",
                     project_id %in% c("P20_1") & mfd_sampletype == "Other" ~ "AAU-LANDFILL",
                     project_id %in% c("P20_1") & mfd_sampletype == "Water" ~ "AAU-LANDFILL",
                     project_id %in% c("P20_1") & mfd_sampletype == "Sediment" ~ "MFD-ENRICHMENT"), 
         .before = "extraction_method") %>%
  arrange(desc(n_sum)) 

comb.sum$n_sum %>% sum()

# Write to output
data.table::fwrite(comb.sum, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_summary_table.csv"))

# Import version of table with sampling method assigned 
# methods <- readxl::read_excel("output/2024-11-20-summary-table-method.xlsx") %>%
#   select(project_id, sampling_method, mfd_sampletype, mfd_areatype) %>%
#   filter(project_id %in% comb.sum$project_id)

# Subset based on sample type
sum.sediment <- comb.sum %>%
  filter(mfd_sampletype == "Sediment") %>%
  arrange(desc(n_sum))

sum.water <- comb.sum %>%
  filter(mfd_sampletype == "Water") %>%
  arrange(desc(n_sum))

sum.soil <- comb.sum %>%
  filter(mfd_sampletype == "Soil") %>%
  arrange(desc(n_sum))

sum.other <- comb.sum %>%
  filter(mfd_sampletype == "Other") %>%
  arrange(desc(n_sum))

### Filter based on methods used
filter <- c("P06_1", "P06_2", "P06_3", "P12_1", "P12_2", 
            "P16_4", "P19_1", "P20_1", "P25_1",
            "P11_1", "P16_1", "P16_2", "P16_3")

# Filter summary table
res.sum <- comb.sum %>%
  filter(!project_id %in% filter)

res.filt <- comb.sum %>%
  filter(project_id %in% filter)

projects <- res.sum %>%
  select(project_id) %>%
  distinct()

# Write to output
data.table::fwrite(projects, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_filtered_projects.csv"))

# Subset based on sample type
res.sediment <- res.sum  %>%
  filter(mfd_sampletype == "Sediment") %>%
  arrange(desc(n_sum))

res.water <- res.sum  %>%
  filter(mfd_sampletype == "Water") %>%
  arrange(desc(n_sum))

res.soil <- res.sum  %>%
  filter(mfd_sampletype == "Soil") %>%
  arrange(desc(n_sum))

res.other <- res.sum  %>%
  filter(mfd_sampletype == "Other") %>%
  arrange(desc(n_sum))

# Result of filtering
sum(sum.sediment$n_sum)
sum(res.sediment$n_sum)

sum(sum.water$n_sum)
sum(res.water$n_sum)

sum(sum.soil$n_sum)
sum(res.soil$n_sum)

sum(sum.other$n_sum)
sum(res.other$n_sum)

sum(sum.sediment$n_sum, sum.soil$n_sum, sum.water$n_sum, sum.other$n_sum)
sum(res.sediment$n_sum, res.soil$n_sum, res.water$n_sum, res.other$n_sum)

setdiff(sum.soil$project_id, res.soil$project_id)
setdiff(sum.sediment$project_id, res.sediment$project_id)
setdiff(sum.water$project_id, res.water$project_id)
setdiff(sum.other$project_id, res.other$project_id)

## Write to output
data.table::fwrite(sum.soil, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_summary_soil.csv"))
data.table::fwrite(sum.sediment, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_summary_sediment.csv"))
data.table::fwrite(sum.water, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_summary_water.csv"))
data.table::fwrite(sum.other, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_summary_other.csv"))

data.table::fwrite(res.soil, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_summary_soil_filtered.csv"))
data.table::fwrite(res.sediment, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_summary_sediment_filtered.csv"))
data.table::fwrite(res.water, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_summary_water_filtered.csv"))
data.table::fwrite(res.other, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("output/", format(Sys.time(), "%Y-%m-%d"), "_summary_other_filtered.csv"))



