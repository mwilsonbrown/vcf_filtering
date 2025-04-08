#!/bin/bash --login
#SBATCH --job-name=CO_prefilter
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
#
# Post-alignment; pre-filtering protocol
# M. Wilson Brown
# Goal: to organize and streamline vcf data before processing through filtering, check indidivuals included

# Variables
PREFIX=CO_on_CBP-CG
RAW_VCF=/mnt/scratch/wils1582/BIG_vcfs/Orientalis_on_Grandiflora.all.vcf
WORKDIR=/mnt/scratch/wils1582/"$PREFIX"_filtering
OUTDIR="$WORKDIR"/sample_files
#ALLSITES_VCF=/mnt/research/josephslab/Adrian/CBP_NYC_JLv4/CBP_JLv4_v_CBP.merged.v.all.vcf

# move to working directory
mkdir -p "$OUTDIR"
cd "$WORKDIR"

# Modules
module purge
module load BCFtools/1.19-GCC-13.2.0
module load  R/4.3.3-gfbf-2023b

# get list of samples in each VCF
bcftools query -l $RAW_VCF > "$OUTDIR"/pre_samples.txt

# Run R script to generate files for subsetting vcf
Rscript /mnt/home/wils1582/vcf_filtering/001a_check_raw_vcf.R "$OUTDIR"/pre_samples.txt "$PREFIX" "$OUTDIR"

# reneame samples
bcftools reheader -s "$OUTDIR"/"$PREFIX"_new_names.txt $RAW_VCF -o "$PREFIX".vcf
