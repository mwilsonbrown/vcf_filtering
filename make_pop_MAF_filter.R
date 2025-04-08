# Make population mapping filter based on individuals in final VCF
# mwilsonbrown
# 2025-04-03
#
# usage
# Rscript assign_species.R prefix
# args[1] is bcftools sample list
# args[2] is full path to metadata
# args[3] is output prefix
args = commandArgs(trailingOnly=TRUE)

######## VARIABLES--------
prefix <- args[3]
filedir <- paste0("/mnt/scratch/wils1582/", prefix, "_filtering/")
######## ASSIGN SPECIES
library(dplyr)
#read in data
inds <- read.delim(args[1], header = F)
wgs <- read.delim(args[2], header = T)

#join data
all <- left_join(inds, wgs, join_by("V1" == "vcf_sample_name"))

# needs to be here for BCFtools +split
all$col2 <- "-"

# make column codes
all <- all %>% mutate(sp_code = case_when(species %in% c("Capsella bursa-pastoris") ~ "CBP",
                                   species %in% c("Capsella rubella") ~ "CR",
                                   species %in% c("Capsella grandiflora") ~ "CG"))
# subset data
cutted <- all[,c("V1", "col2", "sp_code")]
#file output name
out <- paste0(prefix, "_species.txt")

# write table  
write.table(cutted, file = out,
            sep = "\t", row.names = F, quote = F, col.names = F)
