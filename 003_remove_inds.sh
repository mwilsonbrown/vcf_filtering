#!/bin/bash
#
#SBATCH --job-name=RemoveIndiviuals
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-10:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A.out
# VCF filtering
# Maya Wilson Brown
#
# Modified from C. Fiscus Capsella vcf filtering

# Update directories to not run into space issues
WORKDIR=/mnt/scratch/wils1582/CBP_CRCGCONP_filtering


# make sure working directory exists
#mkdir -p "$WORKDIR"

# change directory
cd "$WORKDIR"/

# purge modules
module purge
# load modules
module load BCFtools/1.19-GCC-13.2.0

## Other vars
PREFIX=CBP_CRCGCONP


bcftools view -s ^SRR5803007.sam,ERR636135.sam "$WORKDIR"/"$PREFIX"_filter3.vcf.gz \
	-Oz -o /mnt/research/josephslab/Maya/"$PREFIX"_filtered.vcf.gz

tabix /mnt/research/josephslab/Maya/"$PREFIX"_filtered.vcf.gz

# Make species group file
#cp ~/vcf_filtering/individuals/final_CBP_CRCGCONP_filtered_inds.txt ./

# add Rscript make_pop_MAF_filter.R here

## Filter VCF on species minor allele freq
# first, split the VCF into species groups
bcftools +split --groups-file /mnt/home/wils1582/vcf_filtering/individuals/ \
  /mnt/research/josephslab/Maya/"$PREFIX"_filtered.vcf.gz
  
# Then, have to recalculate the Allele count, Allele Frequency, and Allele Number on each separately
# before filtering on MAF
# Not filtering NP or C. orientalis because there are not that many samples
# I could probably pipe the fill into the view command but it's fine
bcftools +fill-tags CBP.vcf -Oz -o CBP_recalc.vcf.gz -- -t AC,AF,AN
bcftools +fill-tags CR.vcf -Oz -o CR_recalc.vcf.gz -- -t AC,AF,AN
bcftools +fill-tags CG.vcf -Oz -o CG_recalc.vcf.gz -- -t AC,AF,AN

# Filter on MAF for each
# C. grandiflora threshold higher bc outcrosser and eQTLs from EJ paper
# Not filtering NP or C. orientalis because there are not that many samples
bcftools view --min-af 0.025:minor CBP_recalc.vcf.gz -o CBP_maf.vcf.gz -Oz | tabix
bcftools view --min-af 0.025:minor CR_recalc.vcf.gz -o CR_maf.vcf.gz -Oz | tabix
bcftools view --min-af 0.05:minor CG_recalc.vcf.gz -o CG_maf.vcf.gz -Oz | tabix

# bgzip the other two vcfs
bgzip NP.vcf
bgzip CO.vcf

# Index all of them
tabix CBP_maf.vcf.gz 
tabix CR_maf.vcf.gz 
tabix CG_maf.vcf.gz 
tabix CO.vcf.gz 
tabix NP.vcf.gz
# merge all the VCFs back together; one final time, we recalculate the tags
bcftools merge CBP_maf.vcf.gz CR_maf.vcf.gz CG_maf.vcf.gz CO.vcf.gz NP.vcf.gz -Ou | bcftools +fill-tags -Oz -o CBP_CRCGCONP_maf_final_filtered.vcf.gz -- -t AC,AF,AN

