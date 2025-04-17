#!/bin/bash --login
#
#SBATCH --job-name=COonCBP_mainFilters
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-30:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A.out
# I want to benchmark how many snps remain after joining the C. orientalis VCF with the CG/CR/CBP VCF
# just to know how much information loss I expect.
#
# Variables
WORKDIR=
PREFIX1=CO_on_CBP-CG
PREFIX2=CBP2_CR_CG_on_CBP-CG
VCF1=/mnt/scratch/wils1582/"$PREFIX1"_filtering/"$PREFIX1"_qualfiltered.vcf.gz
VCF2=/mnt/scratch/wils1582/"$PREFIX2"_filtering/"$PREFIX2"_qualfiltered.vcf.gz

# change directory
cd $WORKDIR

# purge modules
module purge
# load modules
module load BCFtools/1.18-GCC-12.3.0
module load PLINK/2.00a3.7-gfbf-2023a
module load R/4.3.2-gfbf-2023a

# Join the vcfs
bcftools merge --merge snp-ins-del $VCF1 $VCF2 -Oz -o merged_snps_only_CBP_CR_CG_CO_on_CBP-CG.vcf.gz

# use PLINK2 to calculate site missing values
plink2 --vcf merged_snps_only_CBP_CR_CG_CO_on_CBP-CG.vcf.gz \
 --missing variant-only vcols=chrom,pos,nmiss,nobs,fmiss \
 --allow-extra-chr \
 --double-id \
 --out merged_snps_only_CBP_CR_CG_CO_on_CBP-CG

# just documenting if I want to do the same on the VCF with invariant CO sites (I think that I do)
bcftools merge --merge snp-ins-del /mnt/scratch/wils1582/"$PREFIX1"_filtering/"$PREFIX1"_filter2.vcf.gz $VCF2 -Oz -o merged_CBP_CR_CG_CO_on_CBP-CG.vcf.gz

plink2 --vcf merged_CBP_CR_CG_CO_on_CBP-CG.vcf.gz \
 --missing variant-only vcols=chrom,pos,nmiss,nobs,fmiss \
 --allow-extra-chr \
 --double-id \
 --out merged_CBP_CR_CG_CO_on_CBP-CG
