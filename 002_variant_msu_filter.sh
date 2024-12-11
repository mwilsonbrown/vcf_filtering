#!/bin/bash --login
#SBATCH --job-name=VarMSUVCF_filter
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-50:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A_%a.out
# VCF filtering of MSU all sites VCF alignment to make a normal
# variant only VCF
# December 11, 2024
# Maya Wilson Brown
#
#
# change directory
cd /mnt/scratch/wils1582/cbp_variant_msu_filtering/

### Varibles
VCF=/mnt/scratch/wils1582/allSites_msu_filtering_2/CBP_AllSites_msu_filter3.vcf.gz
PREFIX=CBP_variant_msu
FILEDIR=/mnt/home/wils1582/vcf_filtering

# purge modules
module purge
# load modules
module load BCFtools/1.18-GCC-12.3.0
module load PLINK/2.00a3.7-gfbf-2023a
module load R/4.3.2-gfbf-2023a

####### SETUP
# I do not have read write permissons for my github repo in home
# so I will copy needed files to the working directory
cp "$FILEDIR"/individuals/* ./
cp "$FILEDIR"/*.R ./

#### PIPELINE #####

# start fitering report
touch "$FILEDIR"/"$PREFIX"_log.txt

echo "$PREFIX" >> "$PREFIX"_log.txt
echo date >> "$PREFIX"_log.txt

## GATK best practices hard filters
bcftools filter -e'QD < 2 | FS > 60 | SOR > 3 | MQ < 40 | MQRankSum < -12.5 | ReadPosRankSum < -8.0' \
"$VCF" \
> "$PREFIX"_temp.vcf
bcftools view -f.,PASS "$PREFIX"_temp.vcf \
-Oz -o "$PREFIX"_filter1.vcf.gz

echo "Filter1 complete"

echo "GATK best practices filter" >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter1.vcf.gz | wc -l) \
>> "$PREFIX"_log.txt
echo "Sample count" >> "$PREFIX"_log.txt
echo $(bcftools query -l "$PREFIX"_filter1.vcf.gz | wc -l) \
>> "$PREFIX"_log.txt

# remove temp
rm "$PREFIX"_temp.vcf

## require 3 reads to call and keep only biallelic sites; dump entirely missing sites
bcftools filter -e'FMT/DP<3' -S . "$PREFIX"_filter1.vcf.gz | bcftools view -i 'F_MISSING<1' -m2 -M2 -Oz -o "$PREFIX"_temp2.vcf

# log progress
echo "temp2 complete"

echo "3 reads called and biallelic sites" >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp2.vcf | wc -l) \
>> "$PREFIX"_log.txt

##########
# filter sites with > 5% het & > 5% missing data
## calculate proportion het per site with plink
plink2 --vcf "$PREFIX"_temp2.vcf \
--geno-counts cols=chrom,pos,ref,alt,homref,refalt,homalt1 \
--allow-extra-chr \
--double-id \
--out "$PREFIX"

echo "plink genotype caluclation complete"

## list of site ids for following Rscript
bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp2.vcf > "$PREFIX"_sites.txt

## Produce list of sites to filter by allele frquency
Rscript 007a_het_sites.R "$PREFIX".gcount 0.05 "$PREFIX"

# highly het sites to  remove file
cut -f1,2 "$PREFIX"_hetmin.txt | tail -n+2 > remove.txt

## remove highly heterozygous sites and sites with 5% missing calls
bcftools view --targets-file ^remove.txt "$PREFIX"_temp2.vcf -o "$PREFIX"_temp3.vcf

# log
echo "heterozygisty filter complete"

echo "Heterozygosity site filter" >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp3.vcf | wc -l) \
>> "$PREFIX"_log.txt

# calculate site missingness
plink2 --vcf "$PREFIX"_temp3.vcf --missing variant-only vcols=chrom,pos,nmiss,nobs,fmiss \
 --allow-extra-chr \
 --double-id \
 --out "$PREFIX"

# Keep sites where at least 95% is not missing
bcftools view --include 'F_MISSING<0.05' "$PREFIX"_temp3.vcf -Oz -o "$PREFIX"_filter2.vcf.gz

echo "missing data filter complete"

echo "missing site call filter" >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter2.vcf.gz | wc -l) \
>> "$PREFIX"_log.txt

##########
# top cut depth at value of > Q3 + 1.5IQR
## calculate depth per site and plot
bcftools query -f '%CHROM %POS %DP\n' "$PREFIX"_filter2.vcf.gz \
> depth.txt
Rscript 007b_depth_thresh.R depth.txt "$PREFIX"

echo "depth calc complete"

## filter on depth
MAXDP=$(head -n1 "$PREFIX"_depth_cutval.txt)
MINDP=$(tail -n1 "$PREFIX"_depth_cutval.txt)

bcftools view -i "INFO/DP < "$MAXDP" & INFO/DP > "$MINDP"" "$PREFIX"_filter2.vcf.gz \
	-Oz -o "$PREFIX"_filter3.vcf.gz

## Write to filtering report
echo "Depth filter complete"
echo "depth filter " >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter3.vcf.gz | wc -l) \
	>> "$PREFIX"_log.txt

##########
# calculate stats per ind
plink2 --vcf "$PREFIX"_filter3.vcf.gz \
	--sample-counts cols=homref,het,homalt \
	--allow-extra-chr \
        --double-id \
        --out "$PREFIX"
        
#### Then move on to individual and population filters
