# VCF Filtering plots
# mwilsonbrown
# 2024-06-20
# 
# Altered from:
# cjfiscus; 2020-09-12
#
#
######## VARIABLES--------
prefix <- "CBP2_CR_CG_on_CBP-CG"
filedir <- paste0("/mnt/scratch/wils1582/", prefix, "_filtering/")
#### LIBBABIES----
library(dplyr)
library(ggplot2)
############ SELFING SPECIES VARIANT HETEROZYGOSITY --------
# Result from 002a_het_sites.R
# <PREFIX>_frqx.txt is before heterozygosity filtering, but formatted the way I want
# <PREFIX>_hetmin.txt is sites that get removed by the filter

# load data -- over 15 million rows so it takes a minute
# this is a distribution before filtering to 5% heterozygosity
het_frq <- read.csv(paste0(filedir, prefix, "_frqx.txt"),
                   sep = "\t", header = T)
# prepare for plotting
long<- tidyr::pivot_longer(het_frq, cols = c("FRAC_HET", "FRAC_HOMREF", "FRAC_HOMALT"), names_to = "measure")

# plot heterozygosity before filtering
p1<-ggplot(long, aes(x=value, color=measure), alpha=0.5) +
  geom_density() +
  theme_classic() +
  theme(legend.position="top")

out<-paste0(filedir,prefix, "_frqx.jpeg")
ggsave(out, p1, height=3, width=5)

# make dataset after filtering
after_filt <- het_frq[het_frq$FRAC_HET < 0.05,]

long2<- tidyr::pivot_longer(after_filt, cols = c("FRAC_HET", "FRAC_HOMREF", "FRAC_HOMALT"), names_to = "measure")

# plotdis
p2<-ggplot(long2, aes(x=value, color=measure), alpha=0.5) +
  geom_density() +
  theme_classic() +
  theme(legend.position="top")
#out<-paste0(args[3], "_frqx.jpeg")
out<-paste0(filedir,prefix, "_frqx2.jpeg")
ggsave(out, p2, height=3, width=5)

############# VARIANT DEPTH --------
df <- read.delim(paste0(filedir, prefix, "_depth.txt"),
                      header = F, sep = " ")

df$V3 <- as.numeric(df$V3)
thresh <- c(quantile(df$V3, na.rm = TRUE)[4] + 1.5*IQR(df$V3, na.rm = TRUE), quantile(df$V3, na.rm = TRUE)[2] - 1.5*IQR(df$V3, na.rm = TRUE))

dp <- ggplot() + geom_density(data = df, aes(x=V3)) +
  xlim(0,(thresh[1]+1000)) + geom_vline(aes(xintercept = thresh[1]), color="red") + 
  geom_vline(aes(xintercept = thresh[2]), color = "blue") + 
  theme_classic() + xlab("depth") + ggtitle("Capsella bursa-pastoris depth filter - x-axis truncated")
out<-paste0(filedir,prefix, "_depth.jpeg")
ggsave(out, dp, height=3, width=5)

############## VARIANT MISSINGNESS --------
# load data for site missingness
vmiss <- read.csv(paste0(filedir, prefix, ".vmiss"), 
                   sep = "\t", header = T)
# filter of 5% missing
vmiss_filt <- vmiss[which(vmiss$F_MISS < 0.05),]

#plot missing-ness by scaffold
vmiss_p1 <- ggplot() + geom_boxplot(data = vmiss, aes(x=X.CHROM, y = F_MISS))
ggsave(paste0(filedir, prefix, "missing_variants.jpeg"), vmiss_p1, height=4, width=10)

# same plot but after filtering
vmiss_p2 <- ggplot() + geom_boxplot(data = vmiss_filt, aes(x=X.CHROM, y = F_MISS))
ggsave(paste0(filedir, prefix, "missing_variants_filt.jpeg"), vmiss_p2, height=4, width=10)

######################## INDIVIDUAL HETEROZYGOSITY AND MISSINGNESS --------
scount <- read.delim(paste0(filedir, prefix, ".scount"))
vcf_dat <- read.delim("~/capsella_sample_info/generated_mkwb/Capsella_vcf_metadata.txt")

# calulate heterozygosity
scount$prop_het <- scount$HET_CT/(scount$HOM_REF_CT + scount$HOM_ALT_CT + scount$HET_CT)

# join with species data
scount <- left_join(scount, vcf_dat, join_by("X.IID" == "vcf_sample_name"))

# sort by het
scount<-scount[order(scount$prop_het, decreasing=F),]
scount$X.IID<-factor(scount$X.IID, levels=scount$X.IID)

# plot
p1<-ggplot(scount, aes(x=X.IID, y=prop_het, fill = species)) + geom_bar(stat="identity") +
  theme_bw() + coord_flip() +
  theme(axis.text.x= element_text(size = 7))

ggsave(paste0(filedir, prefix, "_het.jpeg"), p1, height=22, width=8)

# Based on the figure, I would remove the highly heterozygous C. rubella sample,
# and the fairly heterozygous C. bursa-pastoris sample that is in between C. grandiflora samples

