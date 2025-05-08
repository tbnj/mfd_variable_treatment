library(tidyverse)

## Load data
lab.meta <- read_csv("data/2023-09-22_corrected_combined_metadata.csv")

mfd.db <- readxl::read_excel("data/2025-04-14_mfd_db.xlsx")

seq.meta <- read_csv("data/2023-10-11_samples_minimal_metadata_collapsed.csv") %>%
  select(fieldsample_barcode, flat_name, before_total_reads, before_total_bases, duplication_rate)

samples <- mfd.db %>%
  filter(!is.na(BioSample)) %>%
  pull(fieldsample_barcode)

comb.meta <- seq.meta %>%
  left_join(lab.meta) %>%
  left_join(mfd.db) %>%
  distinct()

colnames <- c(colnames(lab.meta), "BioSample", "project_id")

## Clean
data <- comb.meta %>%
  # select(fieldsample_barcode, flat_name, project_id, mfd_sampletype, extraction_method, library_method) %>%
  filter(!is.na(mfd_sampletype)) %>%
  distinct() %>%
  mutate(across(extraction_method, ~case_when(project_id %in% c("P20_1") & mfd_sampletype == "Water" ~ "FastDNA-spin-soil",
                                              project_id %in% c("P20_1") & mfd_sampletype == "Sediment" ~ "PowerSoil-Pro-HT",
                                              project_id %in% c("P20_1") & mfd_sampletype == "Other" ~ "PowerSoil-Pro-HT",
                                              project_id == "P13_3" ~ "PowerSoil-Pro",
                                              TRUE~extraction_method))) %>%
  relocate(BioSample, .after = "fieldsample_barcode") %>%
  left_join(seq.meta) %>%
  arrange(fieldsample_barcode) %>%
  mutate(across(extraction_method, ~str_replace(., "FastDNA-spin-soil", "FastDNA-SPIN-Soil"))) %>%
  mutate(across(extraction_method, ~str_replace(., "Phenol-Chloroform-Isoamylacohol", "Phenol-Chloroform-Isoamylalcohol"))) %>%
  mutate(across(extraction_method, ~case_when(project_id %in% c("P11_1", "P13_1", "P13_2") ~ "FastDNA-SPIN-Soil",
                                              project_id %in% c("P06_3", "P16_4") ~ "PowerLyzer-PowerSoil",
                                              project_id == "P16_3" & extraction_method == "PowerWater" ~ "NucliSens-miniMAG",
                                              TRUE~extraction_method))) %>%
  select(all_of(colnames)) %>%
  relocate(c(BioSample, project_id), .after = "fieldsample_barcode")

data.table::fwrite(data, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE,
                   paste0("data/", format(Sys.time(), "%Y-%m-%d"), "_corrected_combined_metadata.csv"))



