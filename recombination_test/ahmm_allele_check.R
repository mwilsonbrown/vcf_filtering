# Checking out allele counts in ancestry hmm file

setwd("~/Documents/PhD/Research/capsella_introgression/")

# read in sample map
ahmm_sample_map <- read.delim("hmm_sample_mapping.txt", sep = " ",
                              header = F)

# joining with human readable information
# load capsella whole genome sequencing information
wgs <- read.csv("~/Documents/PhD/Research/capsella_sample_info/generated_mkwb/Capsella_vcf_metadata.txt", sep = "\t", header = T)
#load C. bursa-pastoris population information
k3pops <- read.csv("~/Documents/PhD/Research/capsella_population_structure/cbp_pop_str.txt", header = T, sep = "\t")

# mutate another column that separates NY and NJ from the rest of the populations
ny_names <- wgs[which(wgs$citation == "R.Panko"),"vcf_sample_name"]

k3pops <- k3pops %>% mutate(k3pop_sm = case_when(vcf_sample_name %in% ny_names ~ "NYC",
                                                 .default = k3population))


sample_info <- left_join(ahmm_sample_map, k3pops, join_by("V1" == "vcf_sample_name"))
input <- read.delim("~/Documents/PhD/Research/exploratory_cbp_analyses/local_ancestry_inference/posterior_data/ahmm_pruned_py3.input",
                    sep = "\t",
                    header = F)
# select tested columns
queried <- input[,8:ncol(input)]
queried %>% rename_with(~ paste0(.x, c("A", "B"), recycle0 = TRUE)
)

 c("A", "B")

# make longer
counts <- pivot_longer()

#allele A
counts_A <- queried[ , c(TRUE,FALSE) ] #select odd columns

counts_B <- queried[ , c(FALSE,TRUE) ] #select odd columns

counts_A <- counts_A %>% rename_with(~ paste0(.x, "_A", recycle0 = TRUE)
)

counts_B <- counts_B %>% rename_with(~ paste0(.x, "_B", recycle0 = TRUE)
)

tog <- cbind(counts_A, counts_B)

t <- tog %>% pivot_longer(cols = everything(),
                     names_to = c(".value", "num"),

