#!/bin/bash --login
#
#SBATCH --job-name=AllSites_Orientalis
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=18:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A.out
# Filter allsites VCF
# May 2, 2025
# Maya Wilson Brown
#
# Goal: Filter Capsella orientalis on CG subgenome while keeping invariant sites to later join with other Capsella species
#
######################## SETUP
# Variables
PREFIX=CO_AllSites_on_CBP-CG
RAW_VCF=/mnt/scratch/wils1582/BIG_vcfs/Orientalis_on_Grandiflora.all.vcf
WORKDIR=/mnt/scratch/wils1582/"$PREFIX"_filtering
OUTDIR="$WORKDIR"/sample_files
#ALLSITES_VCF=/mnt/research/josephslab/Adrian/CBP_NYC_JLv4/CBP_JLv4_v_CBP.merged.v.all.vcf

# move to working directory
mkdir -p "$OUTDIR"
cd "$WORKDIR"

# Modules
# purge modules
module purge
# load modules
module load BCFtools/1.18-GCC-12.3.0
module load PLINK/2.00a3.7-gfbf-2023a
module load R/4.3.2-gfbf-2023a

#### PIPELINE #####
### GATK best practices hard filters
#bcftools filter -e'QD < 2 | FS > 60 | SOR > 3 | MQ < 40 | MQRankSum < -12.5 | ReadPosRankSum < -8.0' \
#-Ou "$RAW_VCF" \
#  | bcftools view -f.,PASS -Oz -o "$PREFIX"_filter1.vcf.gz
#
#echo "Filter1 complete"
#
### start filtering report
#touch "$WORKDIR"/"$PREFIX"_log.txt
#
#echo "$PREFIX" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo "GATK best practices filter" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter1.vcf.gz | wc -l) \
#>> "$WORKDIR"/"$PREFIX"_log.txt
#echo "Sample count" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -l "$PREFIX"_filter1.vcf.gz | wc -l) \
#>> "$WORKDIR"/"$PREFIX"_log.txt

## require 3 reads to call and not missing quality, hard filter; set genotypes that do not pass to missing
## then, dump sites with more than 10% missing calls and Quality under 20; maximum alleles is 2
#bcftools filter -e 'QUAL="." || FMT/DP<3' --set-GTs . "$PREFIX"_filter1.vcf.gz -Ou | \
#  bcftools view -i 'F_MISSING<0.1 || QUAL>20' \
#  -M2 -Oz -o "$PREFIX"_temp2.vcf
#
## log progress
#echo "temp2 complete"
#echo "3 reads called and missing sites" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp2.vcf | wc -l) \
#>> "$WORKDIR"/"$PREFIX"_log.txt
#
########################### filter sites with > 5% het & > 5% missing data 
### calculate proportion het per site with plink
#plink2 --vcf "$PREFIX"_temp2.vcf \
#--geno-counts cols=chrom,pos,ref,alt,homref,refalt,homalt1 \
#--allow-extra-chr \
#--double-id \
#--out "$PREFIX"
#
#echo "plink genotype caluclation complete"
#
### Produce list of sites to filter by heterozygosity
#Rscript ~/vcf_filtering/002a_het_sites.R "$PREFIX".gcount 0.05 "$PREFIX"
#
## highly het sites to  remove file
#cut -f1,2 "$PREFIX"_hetmin.txt | tail -n+2 > remove.txt
#
### remove highly heterozygous sites
#bcftools view --targets-file ^remove.txt "$PREFIX"_temp2.vcf -o "$PREFIX"_temp3.vcf
#
## log
#echo "heterozygisty filter complete"
#
#echo "Heterozygosity site filter" >> "$PREFIX"_log.txt
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp3.vcf | wc -l) \
#>> "$PREFIX"_log.txt
#
## calculate site missingness (with plink2 for plotting)
#plink2 --vcf "$PREFIX"_temp3.vcf --missing variant-only vcols=chrom,pos,nmiss,nobs,fmiss \
# --allow-extra-chr \
# --double-id \
# --out "$PREFIX"
#
## Keep sites where at least 95% is not missing
#bcftools view --include 'F_MISSING<0.05' "$PREFIX"_temp3.vcf -Oz -o "$PREFIX"_temp4.vcf.gz
#
#echo "missing data filter complete"
#
#echo "missing site call filter" >> "$PREFIX"_log.txt
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp4.vcf.gz | wc -l) \
#	>> "$PREFIX"_log.txt
#
###########
## top cut depth at value of > Q3 + 1.5IQR; calculate depth per site and plot
#bcftools query -f '%CHROM %POS %DP\n' "$PREFIX"_temp4.vcf.gz \
#	> "$WORKDIR"/depth.txt

Rscript ~/vcf_filtering/002b_depth_thresh.R "$WORKDIR"/depth.txt "$WORKDIR"/"$PREFIX"

echo "depth calc complete"

# filter on depth
MAXDP=$(head -n1 "$PREFIX"_depth_cutval.txt)
MINDP=$(tail -n1 "$PREFIX"_depth_cutval.txt)

bcftools view -i "INFO/DP < "$MAXDP" & INFO/DP > "$MINDP"" "$PREFIX"_temp4.vcf.gz  -Oz -o "$PREFIX"_filter3.vcf.gz

## Write to filtering report
echo "Depth filter complete"
echo "depth filter " >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter3.vcf.gz | wc -l) \
	>> "$WORKDIR"/"$PREFIX"_log.txt

##########
# calculate stats per ind
plink2 --vcf "$PREFIX"_filter3.vcf.gz \
	--sample-counts cols=homref,het,homalt \
	--allow-extra-chr \
        --double-id \
        --out "$PREFIX"

#################### Filter on Minor Allele Frequency in variant sites
# # change directory
# cd /mnt/scratch/wils1582/allSites_msu_filtering_2/
# 
# ### Varibles
# #VCF="$WORKDIR"/"$PREFIX"_filter4.vcf #previous prefix
# PREFIX=CBP_allSites_msu_maf #new prefix
# 
# # purge modules
# module purge
# # load modules
# module load  BCFtools/1.19-GCC-13.2.0
# 
# # Generate log file
# 
# # merge the invariant and variant sites
# bcftools concat --allow-overlaps cbp_allsites_variant_maf_filt.bcf cbp_allsites_invariant.bcf -Ou -o "$PREFIX".vcf
# # check total number of sites after concatenating
# echo "Final site count" >> "$PREFIX"_log.txt
# bcftools query -f '%CHROM\t%POS\n' "$PREFIX".vcf | wc -l >> "$PREFIX"_log.txt
