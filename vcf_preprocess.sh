# Post-alignment; pre-filtering protocol
# M. Wilson Brown
# Goal: to organize and streamline vcf data before processing through filtering

# Variables
VAR_VCF=/mnt/research/josephslab/Adrian/CBP_NYC_JLv4/CG_subgenome/July14.final_called.v.vcf.gz
ALLSITES_VCF=/mnt/research/josephslab/Adrian/CBP_NYC_JLv4/CBP_JLv4_v_CBP.merged.v.all.vcf
OUTDIR=/mnt/home/wils1582/vcf_filtering/individuals

# move to working directory
cd ~/vcf_filtering

# Modules
module purge
module load BCFtools/1.19-GCC-13.2.0
module load  R/4.3.3-gfbf-2023b

# get list of samples in each VCF
bcftools query -l $VAR_VCF > samples_variant_only.txt
bcftools query -l $ALLSITES_VCF > samples_all_sites.txt

# Run R script to generate files for subsetting vcf
Rscript 001_check_raw_vcf.R samples_variant_only.txt samples_all_sites.txt

# clean up working directry
rm samples_all_sites.txt
rm samples_variant_only.txt

