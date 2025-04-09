#!/bin/bash --login
#
#SBATCH --job-name=IndHetSNPs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-3:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A.out
# VCF filtering
# Maya Wilson Brown
#
# April 9, 2025

# Variables
PREFIX=CO_on_CBP-CG
WORKDIR=/mnt/scratch/wils1582/"$PREFIX"_filtering

# change directory
cd $WORKDIR

# purge modules
module purge
# load modules
module load BCFtools/1.18-GCC-12.3.0

################## INDIVIDUAL + VARIANT TYPE FILTER
# keep only SNPs 
bcftools view -v snps "$WORKDIR"/"$PREFIX"_filter2.vcf.gz \
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
