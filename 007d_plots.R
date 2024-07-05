# VCF Filtering plots
# mwilsonbrown
# 2024-06-20
# 
# Altered from:
# cjfiscus; 2020-09-12
#
#
############ Genotype filtering of highly heterozygous sites-------

# load genotype count data
# usage
# Rscript plot_frqx.R plink.gcount sites.txt 0.1 prefix
# args[1] is plink2 gcount file with cols=chrom,pos,ref,alt,homref,refalt,homalt1
# args[2] is min het threshold for filtered list
# args[3] is output prefix
args = commandArgs(trailingOnly=TRUE)

# # load data -- over 15 million rows so it takes a minute
# # this is a distribution before filtering to 5% heterozygosity
# gcount <- read.csv("~/Documents/PhD/Research/capsella_population_structure/vcf_filtering/NPCRCG_CBP/NPCRCG_CBP2.gcount", 
#                    sep = "\t", header = T)
# 
# # caluclate fractions
# gcount <- 
# # prepare for plotting
# long<- tidyr::pivot_longer(sub, cols = c("FRAC_HET", "FRAC_HOMREF", "FRAC_HOMALT"), names_to = "measure")
# 
# # plotdis
# p1<-ggplot(long, aes(x=value, color=measure), alpha=0.5) +
#   geom_density() +
#   theme_classic() +
#   theme(legend.position="top")
# #out<-paste0(args[3], "_frqx.jpeg")
# out<-paste0("~/Documents/PhD/Research/capsella_population_structure/vcf_filtering/CBPCRCG", "_frqx.jpeg")
# ggsave(out, p1, height=3, width=5)
# 
# long2<- tidyr::pivot_longer(sub1, cols = c("FRAC_HET", "FRAC_HOMREF", "FRAC_HOMALT"), names_to = "measure")
# 
# # plotdis
# p2<-ggplot(long2, aes(x=value, color=measure), alpha=0.5) +
#   geom_density() +
#   theme_classic() +
#   theme(legend.position="top")
# #out<-paste0(args[3], "_frqx.jpeg")
# out<-paste0("~/Documents/PhD/Research/capsella_population_structure/vcf_filtering/CBPCRCG_removed", "_frqx.jpeg")
# ggsave(out, p2, height=3, width=5)

############# plot depth----------------

# # plot with threshold 
# p1<-ggplot(df, aes(x=V3)) + geom_density() +
#   geom_vline(aes(xintercept=thres), color="red") +
#   theme_classic() +
#   xlab("INFO/DP")
# #out<-paste0(args[2], "_depth.jpeg")
# out<-paste0("CBPCRCG", "_depth.jpeg")
# ggsave(out, p1, height=5, width=5)


############## Distribution of Missing Genotypes-----------
# load data for missingness
#smiss <- read.csv("~/Documents/PhD/Research/capsella_population_structure/vcf_filtering/NPCRCG_CBP/NPCRCG.smiss", 
#                  sep = "\t", header = T)
vmiss <- read.csv("~/Documents/PhD/Research/capsella_population_structure/vcf_filtering/NPCRCG_CBP/NPCRCG.vmiss", 
                   sep = "\t", header = T)

# load vcf sample data
vcf_dat <- read.csv("~/Documents/PhD/Research/capsella_sample_info/generated_mkwb/Capsella_vcf_metadata.txt", 
                    sep = "\t", header = T)

# current filter
vmiss_old_filt <- vmiss[which(vmiss$F_MISS < 0.05),]

# new proposed filter of 20% missing
vmiss_new_filt <- vmiss[which(vmiss$F_MISS < 0.2),]


#plot missing-ness by scaffold
vmiss_p1 <- ggplot() + geom_boxplot(data = vmiss, aes(x=X.CHROM, y = F_MISS))
ggsave("./missing_variants2.jpeg", vmiss_p1, height=4, width=6)

# same plot but after filtering
vmiss_p2 <- ggplot() + geom_boxplot(data = vmiss_old_filt, aes(x=X.CHROM, y = F_MISS))
ggsave("./missing_variants_old_filt2.jpeg", vmiss_p2, height=4, width=6)

vmiss_p3 <- ggplot() + geom_boxplot(data = vmiss_new_filt, aes(x=X.CHROM, y = F_MISS))
ggsave("./missing_variants_new_filt2.jpeg", vmiss_p3, height=4, width=6)

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
