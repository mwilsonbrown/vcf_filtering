
library(dplyr)

inds <- read.delim("~/Desktop/final_CBP_CRCGCONP_filtered_inds.txt", header = F)
wgs <- read.delim("~/Documents/PhD/Research/capsella_sample_info/generated_mkwb/Capsella_vcf_metadata.txt")

all <- left_join(inds, wgs, join_by("V1" == "vcf_sample_name"))
all$col2 <- "-"

# make column codes
all <- all %>% mutate(sp_code = case_when(species %in% c("Capsella bursa-pastoris") ~ "CBP",
                                   species %in% c("Capsella rubella") ~ "CR",
                                   species %in% c("Capsella orientalis") ~ "CO",
                                   species %in% c("Capsella grandiflora") ~ "CG",
                                   species %in% c("Neslia paniculata") ~ "NP",))

cutted <- all[,c("V1", "col2", "sp_code")]

write.table(cutted, file = "maf_Capsella_sppop_split.txt",
            sep = "\t",
            ``)
