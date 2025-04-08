#!/bin/bash --login
#
#SBATCH --job-name=IndHetMAF
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-5:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A.out
# VCF filtering
# Maya Wilson Brown
#
# April 3, 2025

# Variables
PREFIX=CBP2_CR_CG_on_CBP-CG
WORKDIR=/mnt/scratch/wils1582/"$PREFIX"_filtering

# change directory
cd $WORKDIR

# purge modules
module purge
# load modules
module load BCFtools/1.18-GCC-12.3.0
module load PLINK/2.00a3.7-gfbf-2023a
module load R/4.3.2-gfbf-2023a

################## INDIVIDUAL + VARIANT TYPE FILTER
# remove individuals with excess heterozygoisty; individuals are identified visually; keep only SNPs 
bcftools view -s ^SRR5803007.sam,ERR636135.sam -v snps "$WORKDIR"/"$PREFIX"_filter2.vcf.gz \
	-Oz -o "$WORKDIR"/"$PREFIX"_qualfiltered.vcf.gz
tabix "$WORKDIR"/"$PREFIX"_qualfiltered.vcf.gz

echo "Quality filters complete"

# Log variant and sample info
echo "High quality SNPs" >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_qualfiltered.vcf.gz | wc -l) \
	>> "$PREFIX"_log.txt
echo "Sample count" >> "$PREFIX"_log.txt
echo $(bcftools query -l "$PREFIX"_qualfiltered.vcf.gz | wc -l) \
	>> "$PREFIX"_log.txt

################## SPLIT SPECIES TO ID MINOR ALLELE FRQ THRESHOLDS
## Make species group file
bcftools query -l "$PREFIX"_qualfiltered.vcf.gz > "$PREFIX"_samples.txt

Rscript make_pop_MAF_filter.R "$PREFIX"_samples.txt ~/capsella_sample_info/generated_mkwb/Capsella_vcf_metadata.txt "$PREFIX"

## Filter VCF on species minor allele freq
# first, split the VCF into species groups
bcftools +split --groups-file "$WORKDIR"/"$PREFIX"_species.txt "$PREFIX"_qualfiltered.vcf.gz

# Splitting the VCF by species introduces invariant sites
# Calculate allele frequencies with PLINK2 for each species VCF
plink2 --vcf CBP.vcf \
  --freq cols=-maybeprovref \
  --allow-extra-chr \
  --double-id \
  --out CBP_split
  
plink2 --vcf CR.vcf \
  --freq cols=-maybeprovref \
  --allow-extra-chr \
  --double-id \
  --out CR_split
  
plink2 --vcf CG.vcf \
  --freq cols=-maybeprovref \
  --allow-extra-chr \
  --double-id \
  --out CG_split
## In R
## Read in allele frquency file caluculated with PLINK
#cbp_frq <- read.delim("./CBP_split.afreq")
## Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
#cbp_var <- cbp_frq[which(cbp_frq$REF_FREQ != 1 & cbp_frq$ALT_FREQ !=1),]
## of the variant sites, make a list of those where either the reference or alternate allele frequecy is below the threshold
#cbp_rm <- cbp_var[which(cbp_var$REF_FREQ < 0.025 | cbp_var$ALT_FREQ < 0.025),]
#
## do the same for Capsella rubella
#cr_frq <- read.delim("./CR_split.afreq")
## Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
#cr_var <- cr_frq[which(cr_frq$REF_FREQ != 1 & cr_frq$ALT_FREQ !=1),]
## of the variant sites, make a list of those where either the reference or alternate allele frequecy is below the threshold
#cr_rm <- cr_var[which(cr_var$REF_FREQ < 0.025 | cr_var$ALT_FREQ < 0.025),]
#
## and the same for Capsella grandiflora but the MAF cut off is more stringent (0.05)
#cg_frq <- read.delim("./CG_split.afreq")
## Only select variant sites (neither the reference nor alternate allele had a frequency of 1)
#cg_var <- cg_frq[which(cg_frq$REF_FREQ != 1 & cg_frq$ALT_FREQ !=1),]
## of the variant sites, make a list of those where either the reference or alternate allele frequecy is below the threshold
#cg_rm <- cg_var[which(cg_var$REF_FREQ < 0.05 | cg_var$ALT_FREQ < 0.05),]
#
## concatenate all removals together
#rm_sites <- rbind(cbp_rm[1:3], cr_rm[1:3], cg_rm[1:3])
#
## remove the duplicates
#rm_sites_nodups <- rm_sites[!duplicated(rm_sites),]
#
## Write to file
#write.table(rm_sites_nodups[1:2], file = "remove_MAF.txt", sep = "\t", header = F, quote = F, row.names = F, col.names = F)
#
#  
## Then, have to recalculate the Allele count, Allele Frequency, and Allele Number on each separately
## before filtering on MAF
## Not filtering NP or C. orientalis because there are not that many samples
## I could probably pipe the fill into the view command but it's fine
#bcftools +fill-tags CBP.vcf -Oz -o CBP_recalc.vcf.gz -- -t AC,AF,AN
#bcftools +fill-tags CR.vcf -Oz -o CR_recalc.vcf.gz -- -t AC,AF,AN
#bcftools +fill-tags CG.vcf -Oz -o CG_recalc.vcf.gz -- -t AC,AF,AN
#
## Filter on MAF for each
## The definition of AC and AF in BCFTools is for the alternate allele only
## I only want to filter on Minor Allele Frequecy on sites where there is a minor allele at all
#
## I cannot figure out how to fix this right now so I think I am just going to use PLINK2
#plink2 --vcf CBP.vcf \
#--freq cols=chrom,pos,id,ref,alt,reffreq,altfreq,nobs \
#--set-all-var-ids @:# \
#--allow-extra-chr \
#--out CBP_split
## C. grandiflora threshold higher bc outcrosser and eQTLs from EJ paper
## Not filtering NP or C. orientalis because there are not that many samples
#bcftools view --min-af 0.025:minor CBP_recalc.vcf.gz -o CBP_maf.vcf.gz -Oz
#bcftools view --min-af 0.025:minor CR_recalc.vcf.gz -o CR_maf.vcf.gz -Oz
#bcftools view --min-af 0.05:minor CG_recalc.vcf.gz -o CG_maf.vcf.gz -Oz
#
## bgzip the other two vcfs
#bgzip NP.vcf
#bgzip CO.vcf
#
## Index all of them
#tabix CBP_maf.vcf.gz 
#tabix CR_maf.vcf.gz 
#tabix CG_maf.vcf.gz 
#tabix CO.vcf.gz 
#tabix NP.vcf.gz
## merge all the VCFs back together; one final time, we recalculate the tags
#bcftools merge CBP_maf.vcf.gz CR_maf.vcf.gz CG_maf.vcf.gz CO.vcf.gz NP.vcf.gz -Ou | bcftools +fill-tags -Oz -o CBP_CRCGCONP_maf_final_filtered.vcf.gz -- -t AC,AF,AN
#
