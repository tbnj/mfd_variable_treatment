# mfd_variable_treatment
The scripts in this repository are part of the [Microflora Danica project](https://github.com/cmc-aau/mfd_wiki/wiki). 
The scripts are used to clean metadata, make summaries of sampling and extraction methodology of the individual sampling projects and to evaluate the effect of the variable treatment using both V4 amplicons and the 16S fragments extracted from the metagenomes. 

Be advised, that for the metagenomic-derived data, the term "OTU" is only used due to format requirements by ampvis2, and they do not represent classical OTUs. 
The generated profiles can be thought of as taxonomic bins. 

## Scripts
### Metadata
`scripts/meta_clean.R` cleans the metadata for typos and inconsistent use of terms. 


`scripts/variable_sampling.R` performs summaries of methodology across sampling projects. 

### V4 amplicon 16S data 
`scripts/drying_test.R` recreates the analysis on the effect of drying of the soils with V4 amplicons sequenced in the Illumina Miseq platform. The script also renders the figures used in the manuscript. 

### Metagenomic 16S data 
`scripts/agri_test.R` recreates the analysis on the effect of drying of the soils with the metagenomci-derived 16S gene fragments sequenced in the Illumina Novaseq6000. The script also renders the figures used in the manuscript. 

## Data
The scripts rely on data files available from the MFD [github](https://github.com/cmc-aau/mfd_metadata) and the MFD Zenodo [repo](https://zenodo.org/records/12605769). 
