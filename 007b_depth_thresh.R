#!/usr/bin/env Rscript
# plot depth from bcftools
# M. Wilson Brown
# November 13, 2024
# Adapted from cjfiscus
# usage
# Rscript plot_dp.R depth.txt out_prefix
# args[1] is file produced by bcftools query -f '%CHROM %POS %DP\n' test.vcf
# args[2] is output prefix
args = commandArgs(trailingOnly=TRUE)

# read in data
df<-read.table(args[1])

# correct class of column
df$V3 <- as.numeric(df$V3)

# hard cut threshold
thresh <- c(quantile(df$V3, na.rm = TRUE)[4] + 1.5*IQR(df$V3, na.rm = TRUE), quantile(df$V3, na.rm = TRUE)[2] - 1.5*IQR(df$V3, na.rm = TRUE))
#u_thresh <-quantile(df$V3)[4] + 1.5*IQR(df$V3)
#l_thresh <- quantile(df$V3)[2] - 1.5*IQR(df$V3)

# write out depth value
out<-paste0(args[2], "_depth_cutval.txt")
write.table(as.numeric(thresh), out, quote=F, row.names=F, col.names=F, sep = "\n")
