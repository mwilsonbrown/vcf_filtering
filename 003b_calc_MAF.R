# Calculating multispecies minor allele frequency


# Because I filtered to biallelic sites, there should only be 2 alleles--REF and ALT
# Therefore, the minor allele is just the less common one
# I also decided that I should only remove minor alleles that are below the minor allele threshold in all three species
# because those are the one that are actually rare and may represent errors

# args[1] is working directory
args = commandArgs(trailingOnly=TRUE)

workdir <- as.character(args[1])
# will work on non hard-coding the input filenames later

# Read in allele frquency file caluculated with PLINK
cbp_frq <- read.delim(paste0(workdir,"/CBP_split.afreq"))
# do the same for Capsella rubella
cr_frq <- read.delim(paste0(workdir,"/CR_split.afreq"))
# and the same for Capsella grandiflora but the MAF cut off is more stringent (0.05)
cg_frq <- read.delim(paste0(workdir,"/CG_split.afreq"))


# Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
cbp_var <- cbp_frq[which(cbp_frq$REF_FREQ != 1 & cbp_frq$ALT_FREQ !=1),]
# Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
cr_var <- cr_frq[which(cr_frq$REF_FREQ != 1 & cr_frq$ALT_FREQ !=1),]
# Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
cg_var <- cg_frq[which(cg_frq$REF_FREQ != 1 & cg_frq$ALT_FREQ !=1),]


# actually, also for this to work, the SAME ALLELE needs to be below the MAF threshold in all three species
# For example: Chr9:139995, Reference allele A is below the MAF threshold in C. grandiflora 
# but alternate allele T is below the MAF threshold in C. bursa-pastoris and C. rubella
# What I want is a list where the reference allele is below the MAF threshold for all pops and the ALT allele is below MAF threshold for all pops

# get sites where the reference allele is the minor allele and below the threshold (they sum to 1, so if it is below the threshold, it is also by definition, the minor allele)
cbp_refm <- cbp_var[which(cbp_var$REF_FREQ < 0.025),]
cr_refm <- cr_var[which(cr_var$REF_FREQ < 0.025),]
cg_refm <- cg_var[which(cg_var$REF_FREQ < 0.05),]

# Then find the shared sites among all three species where the reference allele does not meet the minor allele freq threshold
refm_shared <- reduce(list(cbp_refm[1:2], cr_refm[1:2], cg_refm[1:2]),
                    inner_join)
# Okay, now we can do the same for the sites where the ALT allele does not meet the MAF threshold
cbp_altm <- cbp_var[which(cbp_var$ALT_FREQ < 0.025),]
cr_altm <- cr_var[which(cr_var$ALT_FREQ < 0.025),]
cg_altm <- cg_var[which(cg_var$ALT_FREQ < 0.05),]

#shared 
altm_shared <- reduce(list(cbp_altm[1:2], cr_altm[1:2], cg_altm[1:2]),
                      inner_join)

# concatenate reference and alt removals together
rm_sites <- rbind(altm_shared, refm_shared)

# should be mathematically impossible for there to be duplicates here so no need to worry about that.
# Write to file
write.table(rm_sites, file = paste0(workdir,"/remove_MAF.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
