# Calculating multispecies minor allele frequency
# Read in allele frquency file caluculated with PLINK
cbp_frq <- read.delim("./CBP_split.afreq")

# Because I filtered to biallelic sites, there should only be 2 alleles--REF and ALT
# Therefore, the minor allele is just the less common one

# also, I dont really care about all the other minor alelles, just the ones that are very very rare

cbp_rare <- cbp_frq[which]


# Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
cbp_var <- cbp_frq[which(cbp_frq$REF_FREQ != 1 & cbp_frq$ALT_FREQ !=1),]
# of the variant sites, make a list of those where either the reference or alternate allele frequency is below the threshold
cbp_rm <- cbp_var[which(cbp_var$REF_FREQ < 0.025 | cbp_var$ALT_FREQ < 0.025),]

# do the same for Capsella rubella
cr_frq <- read.delim("./CR_split.afreq")
# Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
cr_var <- cr_frq[which(cr_frq$REF_FREQ != 1 & cr_frq$ALT_FREQ !=1),]
# of the variant sites, make a list of those where either the reference or alternate allele frequency is below the threshold
cr_rm <- cr_var[which(cr_var$REF_FREQ < 0.025 | cr_var$ALT_FREQ < 0.025),]

# and the same for Capsella grandiflora but the MAF cut off is more stringent (0.05)
cg_frq <- read.delim("./CG_split.afreq")
# Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
cg_var <- cg_frq[which(cg_frq$REF_FREQ != 1 & cg_frq$ALT_FREQ !=1),]
# of the variant sites, make a list of those where either the reference or alternate allele frequency is below the threshold
cg_rm <- cg_var[which(cg_var$REF_FREQ < 0.05 | cg_var$ALT_FREQ < 0.05),]

# concatenate all removals together
rm_sites <- rbind(cbp_rm[1:3], cr_rm[1:3], cg_rm[1:3])

# remove the duplicates
rm_sites_nodups <- rm_sites[!duplicated(rm_sites),]

# Write to file
write.table(rm_sites_nodups[1:2], file = "remove_MAF.txt", sep = "\t", header = F, quote = F, row.names = F, col.names = F)
