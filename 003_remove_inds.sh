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
FINALDIR=/mnt/research/josephslab/Maya/capsella/vcf/filtered

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

Rscript ~/vcf_filtering/make_pop_MAF_filter.R "$PREFIX"_samples.txt ~/capsella_sample_info/generated_mkwb/Capsella_vcf_metadata.txt "$PREFIX"

## Filter VCF on species minor allele freq
# first, split the VCF into species groups
bcftools +split --groups-file "$WORKDIR"/"$PREFIX"_species.txt "$PREFIX"_qualfiltered.vcf.gz -o $WORKDIR

# Splitting the VCF by species introduces invariant sites
# Calculate allele frequencies with PLINK2 for each species VCF
plink2 --vcf CBP.vcf \
  --freq cols=chrom,pos,ref,alt,reffreq,altfreq,nobs \
  --allow-extra-chr \
  --double-id \
  --out CBP_split
  
plink2 --vcf CR.vcf \
  --freq cols=chrom,pos,ref,alt,reffreq,altfreq,nobs \
  --allow-extra-chr \
  --double-id \
  --out CR_split
  
plink2 --vcf CG.vcf \
  --freq cols=chrom,pos,ref,alt,reffreq,altfreq,nobs \
  --allow-extra-chr \
  --double-id \
  --out CG_split

## C. grandiflora MAF threshold higher bc outcrosser and eQTLs from EJ paper
Rscript ~/vcf_filtering/003b_multipop_calc_MAF.R $WORKDIR

# Use bcftools to filter on MAF
bcftools view --targets-file ^remove_MAF.txt "$WORKDIR"/"$PREFIX"_qualfiltered.vcf.gz \
	-Oz -o "$FINALDIR"/"$PREFIX"_qualMAF_snps_only.vcf.gz
tabix "$FINALDIR"/"$PREFIX"_qualMAF_snps_only.vcf.gz