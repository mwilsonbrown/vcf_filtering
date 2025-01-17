#!/bin/bash --login
#
#SBATCH --job-name=mafAllSites
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-2:00:00
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
#
#
# purge modules
module purge
# load modules
module load BCFtools/1.18-GCC-12.3.0
module load PLINK/2.00a3.7-gfbf-2023a
module load R/4.3.2-gfbf-2023a

## Other vars
PREFIX=allSites_CBP_msu
VCF=/mnt/research/josephslab/Adrian/CBP_NYC_JLv4/CBP_JLv4_v_CBP.msu.merged.v.all.vcf
WORKDIR=/mnt/scratch/wils1582
INFILES=/mnt/home/wils1582/vcf_filtering

# In the case that I do not have execute permissions for my own github repo
# copy those files to the current directory
cp "$INFILES"/*.R "$WORKDIR"
cp "$INFILES"/individuals/*.txt "$WORKDIR"

cd "$WORKDIR"
#### PIPELINE #####
### GATK best practices hard filters
bcftools filter -e'QD < 2 | FS > 60 | SOR > 3 | MQ < 40 | MQRankSum < -12.5 | ReadPosRankSum < -8.0' \
"$VCF" \
> "$PREFIX"_temp.vcf
bcftools view -f.,PASS "$PREFIX"_temp.vcf \
-Oz -o "$PREFIX"_filter1.vcf.gz

echo "Filter1 complete"

## start fitering report
touch "$WORKDIR"/"$PREFIX"_log.txt

echo "$PREFIX" >> "$WORKDIR"/"$PREFIX"_log.txt
echo "GATK best practices filter" >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter1.vcf.gz | wc -l) \
>> "$WORKDIR"/"$PREFIX"_log.txt
echo "Sample count" >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -l "$PREFIX"_filter1.vcf.gz | wc -l) \
>> "$WORKDIR"/"$PREFIX"_log.txt

## require 3 reads to call and keep only biallelic sites; dump entirely missing sites
bcftools filter -e'FMT/DP<3' -S . "$PREFIX"_filter1.vcf.gz | bcftools view -i 'F_MISSING<1' -m2 -M2 -Oz -o "$PREFIX"_temp2.vcf

# log progress
echo "temp2 complete"

echo "3 reads called and biallelic sites" >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp2.vcf | wc -l) \
>> "$WORKDIR"/"$PREFIX"_log.txt

########### filter sites with > 5% het & > 5% missing data
## calculate proportion het per site with plink
plink2 --vcf "$PREFIX"_temp2.vcf \
--geno-counts cols=chrom,pos,ref,alt,homref,refalt,homalt1 \
--allow-extra-chr \
--double-id \
--out "$PREFIX"

echo "plink genotype caluclation complete"
#
# Most invariant sites are heterozygous calls so we do not filter those sites out

# calculate site missingness
plink2 --vcf "$PREFIX"_temp2.vcf --missing variant-only vcols=chrom,pos,nmiss,nobs,fmiss \
  --allow-extra-chr \
  --double-id \
  --out "$PREFIX"

# Keep sites where at least 95% is not missing
bcftools view --include 'F_MISSING<0.05' "$PREFIX"_temp2.vcf -Oz -o "$PREFIX"_filter2.vcf.gz

echo "missing data filter complete"

echo "missing site call filter" >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter2.vcf.gz | wc -l) \
	>> "$WORKDIR"/"$PREFIX"_log.txt

##########
# top cut depth at value of > Q3 + 1.5IQR
## calculate depth per site and plot
bcftools query -f '%CHROM %POS %DP\n' "$PREFIX"_filter2.vcf.gz \
	> "$WORKDIR"/depth.txt
#Rscript "$WORKDIR"/007b_depth_thresh.R "$WORKDIR"/depth.txt "$WORKDIR"/"$PREFIX"
Rscript 007b_depth_thresh.R "$WORKDIR"/depth.txt "$WORKDIR"/"$PREFIX"

echo "depth calc complete"

## filter on depth
MAXDP=$(head -n1 "$WORKDIR"/"$PREFIX"_depth_cutval.txt)

bcftools view -i "INFO/DP < $MAXDP" "$PREFIX"_filter2.vcf.gz \
	-Oz -o "$PREFIX"_filter3.vcf.gz

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
# change directory
cd /mnt/scratch/wils1582/allSites_msu_filtering_2/

### Varibles
VCF=/mnt/scratch/wils1582/allSites_msu_filtering_2/CBP_AllSites_msu_filter3.vcf.gz
PREFIX=CBP_AllSites_msu_maf

# purge modules
module purge
# load modules
module load  BCFtools/1.19-GCC-13.2.0

# Generate log file
touch ./"$PREFIX"_log.txt

# Select invariant sites that are all missing or reference (alternative allele count is zero)
bcftools view -i 'INFO/AC == 0' $VCF -Ou | bcftools sort - -Ou -o cbp_allsites_invariant_ref_miss.bcf

# check the number of variants (I expect ~31k based on troubleshooting)
echo "Invariant sites that are all reference or missing calls" >> "$PREFIX"_log.txt
bcftools query -f '%CHROM\t%POS\n' cbp_allsites_invariant_ref_miss.bcf | wc -l >> "$PREFIX"_log.txt

# select variant sites (i.e. the alternative allele count is greater than 0)
bcftools view -i 'INFO/AC > 0' $VCF -Ou -o cbp_allsites_variant.bcf

# log number of variant sites
echo "Variant sites" >> "$PREFIX"_log.txt
bcftools query -f '%CHROM\t%POS\n' cbp_allsites_variant.bcf | wc -l >> "$PREFIX"_log.txt

# Filter for MAF on variant sites only
bcftools view -e 'MAF < 0.05' cbp_allsites_variant.bcf -Ou | bcftools sort - -Ou -o cbp_allsites_variant_maf_filt.bcf

# log number of variant sites after filtering
echo "Variant sites after MAF filter" >> "$PREFIX"_log.txt
bcftools query -f '%CHROM\t%POS\n' cbp_allsites_variant_maf_filt.bcf | wc -l >> "$PREFIX"_log.txt

# tabix index both intermediate files
tabix cbp_allsites_variant_maf_filt.bcf
tabix cbp_allsites_invariant_ref_miss.bcf

# merge the invariant and variant sites
bcftools concat --allow-overlaps \
  cbp_allsites_variant_maf_filt.bcf cbp_allsites_invariant_ref_miss.bcf \
  -Oz -o "$PREFIX".vcf.gz 

# check total number of sites after concatenating
echo "Final site count"
bcftools query -f '%CHROM\t%POS\n' "$PREFIX".vcf.gz | wc -l >> "$PREFIX"_log.txt


