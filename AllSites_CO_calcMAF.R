# Calculating multispecies minor allele frequency

# Because I filtered to biallelic sites, there should only be 2 alleles--REF and ALT
# Therefore, the minor allele is just the less common one

# args[1] is working directory
args = commandArgs(trailingOnly=TRUE)

workdir <- as.character(args[1])
prefix <- as.character(args[2])
# will work on non hard-coding the input filenames later

# Read in allele frquency file caluculated with PLINK
frq <- read.delim(paste0(workdir,"/", prefix, ".afreq"))

# Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
var <- frq[which(frq$REF_FREQ != 1 & frq$ALT_FREQ !=1),]

# get sites where the reference allele is the minor allele and below the threshold (they sum to 1, so if it is below the threshold, it is also by definition, the minor allele)
below_minor <- var[which(var$REF_FREQ < 0.025 | var$ALT_FREQ < 0.025),]

# Write to file
write.table(below_minor[1:2], file = paste0(workdir,"/remove_MAF.txt"), sep = "\t", quote = F, row.names = F, col.names = F)

# select invariant sites for easy VCF splitting
invar <- frq[which(frq$REF_FREQ == 1 | frq$ALT_FREQ ==1),]

# write table
write.table(invar[1:2], file = paste0(workdir,"/invariant_sites.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
