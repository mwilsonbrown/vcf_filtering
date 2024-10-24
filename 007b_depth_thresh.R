#!/usr/bin/env Rscript
# plot depth from bcftools
# cjfiscus
# 2020-09-16
# usage
# Rscript plot_dp.R depth.txt out_prefix
# args[1] is file produced by bcftools query -f '%CHROM %POS %DP\n' test.vcf
# args[2] is output prefix
args = commandArgs(trailingOnly=TRUE)

# read in data
df<-read.table(args[1])

# hard cut threshold
thresh <- c(quantile(df$V3)[4] + 1.5*IQR(df$V3), quantile(df$V3)[2] - 1.5*IQR(df$V3))
#u_thresh <-quantile(df$V3)[4] + 1.5*IQR(df$V3)
#l_thresh <- quantile(df$V3)[2] - 1.5*IQR(df$V3)

# write out depth value
out<-paste0(args[2], "_depth_cutval.txt")
write.table(as.numeric(thresh), out, quote=F, row.names=F, col.names=F, sep = "\n")