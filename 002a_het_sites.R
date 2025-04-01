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
#library(data.table)

# read in .acount from plink2
#df<-fread(args[1])
df <- read.delim(args[1])

df$TOTAL<-df$HOM_REF_CT+df$HET_REF_ALT_CTS+df$HOM_ALT1_CT #sum allele freq
df$FRAC_HET<-df$HET_REF_ALT_CTS/df$TOTAL # calculate fractions
df$FRAC_HOMREF<-df$HOM_REF_CT/df$TOTAL 
df$FRAC_HOMALT<-df$HOM_ALT1_CT/df$TOTAL


# get allele frequency fractions
sub<-df[,c("X.CHROM","POS","FRAC_HET", "FRAC_HOMREF", "FRAC_HOMALT")]

#write out
out<-paste0(args[3], "_frqx.txt")
#fwrite(sub, out, quote=F, sep="\t")
write.table(sub, file = out, quote = F, sep = "\t", row.names = F, col.names = T)

# get rows with heterozygosity exceeding threshold set in arguments,
# these individuals will get removed
sub1<-sub[which(sub$FRAC_HET > as.numeric(args[2])),]
out1<-paste0(args[3], "_hetmin",".txt")
#fwrite(sub1, out, quote=F, sep="\t")
write.table(sub1, file = out1, quote = F, sep = "\t", row.names = F, col.names = T)
