# Post-alignment; pre-filtering protocol
# M. Wilson Brown
# Goal: to organize and streamline vcf data before processing through filtering, check indidivuals included

# Variables
PREFIX=CBP2_CR_CG_on_CBP-CG
VAR_VCF=/mnt/scratch/wils1582/BIG_vcfs/merged_final_CBP2_CR_CG_on_CBP-CG.v.vcf
OUTDIR=/mnt/home/wils1582/vcf_filtering/individuals
WORKDIR=/mnt/scratch/wils1582
#ALLSITES_VCF=/mnt/research/josephslab/Adrian/CBP_NYC_JLv4/CBP_JLv4_v_CBP.merged.v.all.vcf

# move to working directory
mkdir "$WORKDIR"/"$PREFIX"_filtering
cd "$WORKDIR"/"$PREFIX"_filtering

# Modules
module purge
module load BCFtools/1.19-GCC-13.2.0
module load  R/4.3.3-gfbf-2023b

# get list of samples in each VCF
bcftools query -l $VAR_VCF > pre_samples.txt

# Run R script to generate files for subsetting vcf
Rscript 001a_check_raw_vcf.R pre_samples.txt "$PREFIX"
