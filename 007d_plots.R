# VCF Filtering plots
# mwilsonbrown
# 2024-06-20
# 
# Altered from:
# cjfiscus; 2020-09-12
#
#

############ Genotype filtering of highly heterozygous sites-------
# Result from 007a_het_sites.R
# load genotype count data
# usage
# Rscript plot_frqx.R plink.gcount sites.txt 0.1 prefix
# args[1] is plink2 gcount file with cols=chrom,pos,ref,alt,homref,refalt,homalt1
# args[2] is min het threshold for filtered list
# args[3] is output prefix

# load data -- over 15 million rows so it takes a minute
# this is a distribution before filtering to 5% heterozygosity
gcount <- read.csv("~/Documents/PhD/Research/capsella_population_structure/vcf_filtering/NPCRCG_CBP/NPCRCG_CBP2.gcount",
                   sep = "\t", header = T)

# caluclate fractions
gcount$TOTAL<-df$HOM_REF_CT+df$HET_REF_ALT_CTS+df$HOM_ALT1_CT #sum allele freq
gcount$FRAC_HET<-df$HET_REF_ALT_CTS/df$TOTAL # calculate fractions
gcount$FRAC_HOMREF<-df$HOM_REF_CT/df$TOTAL 
gcount$FRAC_HOMALT<-df$HOM_ALT1_CT/df$TOTAL

# prepare for plotting
long<- tidyr::pivot_longer(sub, cols = c("FRAC_HET", "FRAC_HOMREF", "FRAC_HOMALT"), names_to = "measure")

# plot heterozygosity before filtering
p1<-ggplot(long, aes(x=value, color=measure), alpha=0.5) +
  geom_density() +
  theme_classic() +
  theme(legend.position="top")
#out<-paste0(args[3], "_frqx.jpeg")
out<-paste0("~/Documents/PhD/Research/vcf_filtering/plots/",SET, "_frqx.jpeg")
ggsave(out, p1, height=3, width=5)

# make dataset after filtering
gcount2 <- gcount[gcount$FRAC_HET < 0.2,]

long2<- tidyr::pivot_longer(sub1, cols = c("FRAC_HET", "FRAC_HOMREF", "FRAC_HOMALT"), names_to = "measure")

# plotdis
p2<-ggplot(long2, aes(x=value, color=measure), alpha=0.5) +
  geom_density() +
  theme_classic() +
  theme(legend.position="top")
#out<-paste0(args[3], "_frqx.jpeg")
out<-paste0("~/Documents/PhD/Research/capsella_population_structure/vcf_filtering/CBPCRCG_removed", "_frqx.jpeg")
ggsave(out, p2, height=3, width=5)

############# plot depth----------------

p2 <- ggplot() + geom_density(data = cbp_msu, aes(x=V3)) +
  xlim(0,7000) + geom_vline(aes(xintercept = thres_cbp), color="red") + 
  geom_vline(aes(xintercept = l_thresh), color = "blue") + 
  theme_classic() + ggtitle("Capsella bursa-pastoris depth filter - x-axis truncated")


############## Distribution of Missing Genotypes-----------
# load data for site missingness
vmiss <- read.csv("~/Documents/PhD/Research/vcf_filtering/large_files/NPCRCG_CBP2.vmiss", 
                   sep = "\t", header = T)

# load data for site missingness after removing individuals that failed
vmiss_prop <- read.csv("~/Documents/PhD/Research/vcf_filtering/large_files/rmInd_NPCRCGCBP2.vmiss", 
                       sep = "\t", header = T)

# filter of 20% missing
vmiss_filt <- vmiss[which(vmiss$F_MISS < 0.2),]

vmiss_pfilt <- vmiss_prop[which(vmiss_prop$F_MISS < 0.2),]


#plot missing-ness by scaffold
vmiss_p1 <- ggplot() + geom_boxplot(data = vmiss, aes(x=X.CHROM, y = F_MISS))
ggsave("~/Documents/PhD/Research/vcf_filtering/large_files/missing_variants.jpeg", vmiss_p1, height=4, width=6)

# same plot but after filtering
vmiss_p2 <- ggplot() + geom_boxplot(data = vmiss_filt, aes(x=X.CHROM, y = F_MISS))
ggsave("~/Documents/PhD/Research/vcf_filtering/large_files/missing_variants_filt.jpeg", vmiss_p2, height=4, width=6)

# plot with inds removed
vmiss_p3 <- ggplot() + geom_boxplot(data = vmiss_pfilt, aes(x=X.CHROM, y = F_MISS))
ggsave("~/Documents/PhD/Research/vcf_filtering/large_files/missing_variants_rmInd.jpeg", vmiss_p3, height=4, width=6)


# plot missing-ness by species
# join with data
# smiss2 <- left_join(smiss, vcf_dat[,c("vcf_sample_name","sample_name","species","citation")],
#                     join_by("IID" == "vcf_sample_name"))
# 
# smiss2[which(is.na(smiss2$species)),"species"] <- "Capsella bursa-pastoris"
# 
# 
# plot
# smiss_p1 <- ggplot() + geom_violin(data = smiss2, aes(x=species, y = F_MISS,color=species),
#                                    drop = F)
# ggsave("./missing_sample_species.jpeg", smiss_p1, height=4, width=6)

######## Individuals by heterozygosity--------
scount <- read.delim("~/Documents/PhD/Research/capsella_population_structure/vcf_filtering/NPCRCG_CBP/NPCRCG_CBP2.scount")

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
  theme(axis.text.x= element_text(size = 0.5))

ggsave("~/Documents/PhD/Research/capsella_population_structure/vcf_filtering/NPCRCG_CBP/NPCRCG_CBP_het.jpeg", p1, height=20, width=8)
