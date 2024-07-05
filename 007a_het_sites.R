# Site Allele Frequencies 
# mwilsonbrown
# 2024-06-20
# 
# Altered from:
# plot plink.frqx; cjfiscus; 2020-09-12
#
# usage
# Rscript plot_frqx.R plink.gcount sites.txt 0.1 prefix
# args[1] is plink2 gcount file with cols=chrom,pos,ref,alt,homref,refalt,homalt1
# args[2] is min het threshold for filtered list
# args[3] is output prefix
args = commandArgs(trailingOnly=TRUE)

# load libs
library(data.table)
#library(ggplot2)
#library(tidyr)
#library(reshape2)
#p_load(data.table,ggplot2, reshape2)

# read in .acount from plink2
df<-fread(args[1])
#df <- fread("~/Downloads/CBPCRCG.frqx")
df$TOTAL<-df$HOM_REF_CT+df$HET_REF_ALT_CTS+df$HOM_ALT1_CT #sum allele freq
df$FRAC_HET<-df$HET_REF_ALT_CTS/df$TOTAL # calculate fractions
df$FRAC_HOMREF<-df$HOM_REF_CT/df$TOTAL 
df$FRAC_HOMALT<-df$HOM_ALT1_CT/df$TOTAL

# read in sites from bcftools query
#sites<-fread(args[2])
#df1 <- fread("~/Downloads/CBPCRCG_sites.txt")
#names(sites)<-c("CHR", "POS")

# get allele frequency fractions
#sub<-df[,c("FRAC_HET", "FRAC_HOMREF", "FRAC_HOMALT")]
sub<-df[,c("#CHROM","POS","FRAC_HET", "FRAC_HOMREF", "FRAC_HOMALT")]

# combine data and write out
#sub<-as.data.frame(cbind(df1, sub))
out<-paste0(args[3], "_frqx.txt")
fwrite(sub, out, quote=F, sep="\t")
#write.table(sub, file = paste0(args[3], "_frqx.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# get rows with heterozygosity exceeding threshold set in arguments,
# these individuals will get removed
sub1<-sub[sub$FRAC_HET > as.numeric(args[2]),]
out<-paste0(args[3], "_hetmin",".txt")
fwrite(sub1, out, quote=F, sep="\t")
#write.table(sub1, file = paste0(args[3], "_hetmin.txt"), quote = F, sep = "\t", row.names = F, col.names = T)


# prepare for plotting
#long<- tidyr::pivot_longer(sub, cols = c("FRAC_HET", "FRAC_A1", "FRAC_A2"), names_to = "measure")

# plotdis
# p1<-ggplot(long, aes(x=value, color=measure), alpha=0.5) + 
#   geom_density() +
#   theme_classic() +
#   theme(legend.position="top")
# out<-paste0(args[4], "_frqx.jpeg")
# ggsave(out, p1, height=3, width=5)
