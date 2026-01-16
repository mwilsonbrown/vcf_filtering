#!/bin/bash --login
#
#SBATCH --job-name=mafAllSites
#SBATCH --nodes=3
#SBATCH --time=1-18:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A.out
# Filter allsites VCF
# January 14, 2025
# Maya Wilson Brown
#
# Updated Goal: combine the whole VCF filtering procedure into a single script
# Successfully generated filtered VCF but script is untested as single runthrough
# as of Jan 31, 2025
#
######################## SETUP
# purge modules
module purge
# load modules
module load BCFtools/1.18-GCC-12.3.0
module load R/4.3.2-gfbf-2023a

## Other vars
PREFIX=CBP_allSites_msu225
VCF=/mnt/scratch/wils1582/CBP_allSites_msu2_raw.vcf.gz
WORKDIR=/mnt/scratch/wils1582/allSites_msu_filtering
INFILES=/mnt/home/wils1582/vcf_filtering

# In the case that I do not have execute permissions for my own github repo,
# copy those files to the working directory
mkdir -p $WORKDIR

cp "$INFILES"/*.R $WORKDIR
cd $WORKDIR


#### PIPELINE #####
# ### GATK best practices hard filters
# # Mark sites
# bcftools filter -e 'QD < 2 | FS > 60 | SOR > 3 | MQ < 40 | MQRankSum < -12.5 | ReadPosRankSum < -8.0' \
# "$VCF" > "$PREFIX"_temp.vcf
# 
# #filter sites that PASS
# bcftools view -f.,PASS "$PREFIX"_temp.vcf -Oz -o "$PREFIX"_filter1.vcf.gz
# 
# echo "Filter1 complete"
# 
# ## start filtering report
# touch "$WORKDIR"/"$PREFIX"_log.txt
#
#echo "$PREFIX" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo "GATK best practices filter" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter1.vcf.gz | wc -l) \
#	>> "$WORKDIR"/"$PREFIX"_log.txt
#echo "Sample count" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -l "$PREFIX"_filter1.vcf.gz | wc -l) \
#	>> "$WORKDIR"/"$PREFIX"_log.txt
#
## remove no quality sites
#bcftools view -e 'QUAL="."' "$PREFIX"_filter1.vcf.gz -Oz -o "$PREFIX"_temp1.vcf.gz
#
## setting regions in next step may require tabix index of vcf
#tabix "$PREFIX"_temp1.vcf.gz
#
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp1.vcf.gz | wc -l) \
#       >> "$WORKDIR"/"$PREFIX"_log.txt
#
### require 3 reads to call and Quality over 20, hard filter; set genotypes that do not pass to missing;
### then, dump sites with more than 10% missing calls
#bcftools filter -e 'QUAL<20 || FMT/DP<3' --set-GTs . "$PREFIX"_temp1.vcf.gz -Ou | \
#  bcftools view -i 'F_MISSING<0.1' -M2 -Oz -o "$PREFIX"_temp2.vcf
#
## log progress
#echo "temp2 complete"
#echo "3 reads called and missing sites" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp2.vcf | wc -l) \
#>> "$WORKDIR"/"$PREFIX"_log.txt
#
###########
## top cut depth at value of > Q3 + 1.5IQR; calculate depth per site and plot
#bcftools query -f '%CHROM %POS %DP\n' "$PREFIX"_temp2.vcf \
#	> "$WORKDIR"/depth.txt
#
#Rscript 002b_depth_thresh.R "$WORKDIR"/depth.txt "$WORKDIR"/"$PREFIX"
#
#echo "depth calc complete"
#
# filter on depth
MAXDP=$(head -n1 "$PREFIX"_depth_cutval.txt)
MINDP=$(tail -n1 "$PREFIX"_depth_cutval.txt)

bcftools view -i "INFO/DP < "$MAXDP" & INFO/DP > "$MINDP"" "$PREFIX"_temp2.vcf.gz  -Oz -o "$PREFIX"_filter3.vcf.gz

## Write to filtering report
echo "Depth filter complete"
echo "depth filter " >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter3.vcf.gz | wc -l) \
	>> "$WORKDIR"/"$PREFIX"_log.txt

# split invariant and potentially variant sites

# separate out invariant sites with unseen alternative allele
bcftools view -i 'ALT="."' "$PREFIX"_filter3.vcf.gz -Oz -o "$PREFIX"_noALT.vcf.gz

# do the oppposite and keep sites with seen alleles of varying allele frequency
bcftools view -e 'ALT="."' "$PREFIX"_filter3.vcf.gz -Oz -o "$PREFIX"_yesALT.vcf.gz
