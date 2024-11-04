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
module load BCFtools/1.18-GCC-12.3.0
module load tabixpp/1.1.2-GCC-12.3.0

## Other vars
PREFIX=CBP_CRCGCONP


bcftools view -s ^SRR5803007.sam,ERR636135.sam "$WORKDIR"/"$PREFIX"_filter3.vcf.gz \
	-Oz -o /mnt/research/josephslab/Maya/"$PREFIX"_filtered.vcf.gz

tabix /mnt/research/josephslab/Maya/"$PREFIX"_filtered.vcf.gz

# Make species group file
cp ~/vcf_filtering/individuals/final_CBP_CRCGCONP_filtered_inds.txt ./



# set group file Var ID
## Filter VCF on species minor allele freq
bcftools +split --groups-file $GROUPS
