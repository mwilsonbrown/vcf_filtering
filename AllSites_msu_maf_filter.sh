#!/bin/bash --login
#
#SBATCH --job-name=mafAllSites
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
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
# Successfully generated filtered VCF but script is untested as singal runthrough
# as of Jan 31, 2025
#
######################## SETUP
# purge modules
module purge
# load modules
module load BCFtools/1.18-GCC-12.3.0
module load R/4.3.2-gfbf-2023a

## Other vars
PREFIX=CBP_allSites_msu
VCF=/mnt/research/josephslab/Adrian/CBP_NYC_JLv4/CBP_JLv4_v_CBP.msu.merged.v.all.vcf
WORKDIR=/mnt/scratch/wils1582/allSites_msu_filtering
INFILES=/mnt/home/wils1582/vcf_filtering

# In the case that I do not have execute permissions for my own github repo,
# copy those files to the working directory
cp "$INFILES"/*.R "$WORKDIR"
cd "$WORKDIR"

#### PIPELINE #####
### GATK best practices hard filters
# Mark sites
bcftools filter -e 'QD < 2 | FS > 60 | SOR > 3 | MQ < 40 | MQRankSum < -12.5 | ReadPosRankSum < -8.0' \
"$VCF" > "$PREFIX"_temp.vcf

#filter sites that PASS
bcftools view -f.,PASS "$PREFIX"_temp.vcf -Oz -o "$PREFIX"_filter1.vcf.gz

echo "Filter1 complete"

## start filtering report
touch "$WORKDIR"/"$PREFIX"_log.txt

echo "$PREFIX" >> "$WORKDIR"/"$PREFIX"_log.txt
echo "GATK best practices filter" >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter1.vcf.gz | wc -l) \
>> "$WORKDIR"/"$PREFIX"_log.txt
echo "Sample count" >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -l "$PREFIX"_filter1.vcf.gz | wc -l) \
>> "$WORKDIR"/"$PREFIX"_log.txt

# setting regions may require tabix index of vcf
tabix "$PREFIX"_filter1.vcf.gz

## require 3 reads to call and Quality over 20, hard filter; set genotypes that do not pass to missing;
## then, dump sites with more than 10% missing calls
bcftools filter -e 'QUAL<20 || FMT/DP<3' --set-GTs . "$PREFIX"_filter1.vcf.gz -Ou | \
  bcftools view -i 'F_MISSING<0.1' -e 'QUAL="."' -r jlSCF_1,jlSCF_2,jlSCF_3,jlSCF_4,jlSCF_5,jlSCF_6,jlSCF_7,jlSCF_8,jlSCF_9,jlSCF_10,jlSCF_11,jlSCF_12,jlSCF_13,jlSCF_14,jlSCF_15,jlSCF_16 \
  -M2 -Oz -o "$PREFIX"_temp2.vcf

# log progress
echo "temp2 complete"
echo "3 reads called and missing sites" >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp2.vcf | wc -l) \
>> "$WORKDIR"/"$PREFIX"_log.txt

# Most invariant sites are heterozygous calls so we do not filter those sites out

##########
# top cut depth at value of > Q3 + 1.5IQR; calculate depth per site and plot
bcftools query -f '%CHROM %POS %DP\n' "$PREFIX"_temp2.vcf \
	> "$WORKDIR"/depth.txt

Rscript 007b_depth_thresh.R "$WORKDIR"/depth.txt "$WORKDIR"/"$PREFIX"

echo "depth calc complete"

bgzip "$PREFIX"_temp2.vcf
# filter on depth
MAXDP=$(head -n1 "$PREFIX"_depth_cutval.txt)
MINDP=$(tail -n1 "$PREFIX"_depth_cutval.txt)

bcftools view -i "INFO/DP < "$MAXDP" & INFO/DP > "$MINDP"" -r jlSCF_1,jlSCF_2,jlSCF_3,jlSCF_4,jlSCF_5,jlSCF_6,jlSCF_7,jlSCF_8,jlSCF_9,jlSCF_10,jlSCF_11,jlSCF_12,jlSCF_13,jlSCF_14,jlSCF_15,jlSCF_16 "$PREFIX"_temp2.vcf.gz \
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

bcftools view -e 'QUAL="."' CBP_allSites_msu_filter3.vcf.gz -Ou -o "$PREFIX"_filter4.vcf

echo "N0 Quality info filter"
echo "has quality info filter " >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter4.vcf | wc -l) \
       >> "$WORKDIR"/"$PREFIX"_log.txt

bcftools view -i 'ALT="."' "$PREFIX"_filter4.vcf -v snps -Ou -o "$PREFIX"_filter5.vcf


#################### Filter on Minor Allele Frequency in variant sites
# change directory
cd /mnt/scratch/wils1582/allSites_msu_filtering_2/

### Varibles
#VCF="$WORKDIR"/"$PREFIX"_filter4.vcf #previous prefix
PREFIX=CBP_allSites_msu_maf #new prefix

# purge modules
module purge
# load modules
module load  BCFtools/1.19-GCC-13.2.0

# Generate log file

# merge the invariant and variant sites
bcftools concat --allow-overlaps cbp_allsites_variant_maf_filt.bcf cbp_allsites_invariant.bcf -Ou -o "$PREFIX".vcf
# check total number of sites after concatenating
echo "Final site count" >> "$PREFIX"_log.txt
bcftools query -f '%CHROM\t%POS\n' "$PREFIX".vcf | wc -l >> "$PREFIX"_log.txt