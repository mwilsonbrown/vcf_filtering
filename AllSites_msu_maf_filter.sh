#!/bin/bash --login
#
#SBATCH --job-name=FilterOrientalis
#SBATCH --nodes=3
#SBATCH --time=10:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/%A.out
# Filter allsites VCF
# January 14, 2025
# Maya Wilson Brown
#
# Updated Goal: combine the whole VCF filtering procedure into a single script
#
######################## SETUP
# purge modules
module purge
# load modules
module load BCFtools/1.18-GCC-12.3.0
module load R/4.3.2-gfbf-2023a

## Other vars
PREFIX=COonCBP-CR
VCF=/mnt/scratch/wils1582/Orientalis_on_Grandiflora.all.vcf
WORKDIR=/mnt/scratch/wils1582/CO_filtering
INFILES=/mnt/home/wils1582/vcf_filtering

# In the case that I do not have execute permissions for my own github repo,
# copy those files to the working directory
#mkdir -p $WORKDIR

#cp "$INFILES"/*.R $WORKDIR
cd $WORKDIR


###### PIPELINE #####
##### GATK best practices hard filters
### Mark sites
#bcftools filter -e 'QD < 2 | FS > 60 | SOR > 3 | MQ < 40 | MQRankSum < -12.5 | ReadPosRankSum < -8.0' \
#"$VCF" > "$PREFIX"_temp.vcf
#
##filter sites that PASS
#bcftools view -f.,PASS "$PREFIX"_temp.vcf -Oz -o "$PREFIX"_filter1.vcf.gz
#
#echo "Filter1 complete"
#
## start filtering report
#touch "$WORKDIR"/"$PREFIX"_log.txt
#
#echo "$PREFIX" >> "$WORKDIR"/"$PREFIX"_log.txt
#
#echo "GATK best practices filter" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter1.vcf.gz | wc -l) \
#	>> "$WORKDIR"/"$PREFIX"_log.txt
#
#echo "Sample count" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -l "$PREFIX"_filter1.vcf.gz | wc -l) \
#	>> "$WORKDIR"/"$PREFIX"_log.txt
#
## remove 'no quality' sites
#bcftools view -e 'QUAL="."' "$PREFIX"_filter1.vcf.gz -Oz -o "$PREFIX"_temp1.vcf.gz
#
## record sites after filtering for no quality
#echo "Sites without Quality scores removed" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp1.vcf.gz | wc -l) \
#      >> "$WORKDIR"/"$PREFIX"_log.txt
#
### require 3 reads to call and Quality over 20, hard filter; set genotypes that do not pass to missing;
### then, dump sites with more than 10% missing calls;
### keep sites with max 2 alleles
#bcftools filter -e 'QUAL<20 || FMT/DP<3' --set-GTs . "$PREFIX"_temp1.vcf.gz -Ou | \
# bcftools view -i 'F_MISSING<0.1' -M2 -Oz -o "$PREFIX"_temp2.vcf
#
## log progress
#echo "temp2 complete"
#echo "3 reads called, site quality, and max 2 alleles" >> "$WORKDIR"/"$PREFIX"_log.txt
#echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp2.vcf | wc -l) \
#>> "$WORKDIR"/"$PREFIX"_log.txt
#
###########
## top cut depth at value of > Q3 + 1.5IQR; calculate depth per site and plot
#bcftools query -f '%CHROM %POS %DP\n' "$PREFIX"_temp2.vcf \
#> "$WORKDIR"/depth.txt
#
#Rscript 002b_depth_thresh.R "$WORKDIR"/depth.txt "$WORKDIR"/"$PREFIX"
#
#echo "depth calc complete"

# filter on depth
MAXDP=$(head -n1 "$PREFIX"_depth_cutval.txt)
MINDP=$(tail -n1 "$PREFIX"_depth_cutval.txt)

bcftools view -i "INFO/DP < "$MAXDP" & INFO/DP > "$MINDP"" "$PREFIX"_temp2.vcf  -Oz -o "$PREFIX"_filter3.vcf.gz

# Write to filtering report
echo "Depth filter complete"
echo "depth filter " >> "$WORKDIR"/"$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_filter3.vcf.gz | wc -l) \
	>> "$WORKDIR"/"$PREFIX"_log.txt



#Part II:
# Goals: To go from single species AllSites base filters (ends with filter3)
# thorugh more site level filters and filtering for minor allele frquency
# The resulting the input and output VCFs from this script are both ALLSITES VCFs
#
# Started: December 1, 2025
# Updated: Mar 4, 2026
######################## SETUP
# purge modules
module purge
# load modules
module load R/4.3.2-gfbf-2023a
module load PLINK/2.00a3.7-gfbf-2023a
module load BCFtools/1.18-GCC-12.3.0 #this version of bcftools works with the other modules
## Other vars
#PREFIX=CBP_allSites_msu225
#WORKDIR=/mnt/scratch/wils1582/allSites_msu_filtering

# In the case that I do not have execute permissions for my own github repo,
# copy those files to the working directory
#cd "$WORKDIR"

###### split invariant and potentially variant sites
# separate out invariant sites with unseen alternative allele
bcftools view -i 'ALT="."' "$PREFIX"_filter3.vcf.gz -Oz -o "$PREFIX"_noALT.vcf.gz
# index for concatentation
bcftools index "$PREFIX"_noALT.vcf.gz
# do the oppposite and keep sites with seen alleles of varying allele frequency
bcftools view -e 'ALT="."' -v snps "$PREFIX"_filter3.vcf.gz -Oz -o "$PREFIX"_yesALT_snps.vcf.gz

# Use plink to calculate genotype counts for site heterozygosity filter
plink2 --vcf "$PREFIX"_yesALT_snps.vcf.gz \
  --geno-counts cols=chrom,pos,ref,alt,homref,refalt,homalt1 \
  --allow-extra-chr \
  --double-id \
  --out "$PREFIX"

echo "plink genotype caluclation complete"

## Produce list of sites to filter by heterozygosity
Rscript 002a_het_sites.R "$PREFIX".gcount 0.05 "$PREFIX"

# highly het sites to remove file
cut -f1,2 "$PREFIX"_hetmin.txt | tail -n+2 > remove_hetsites.txt

## remove highly heterozygous sites
bcftools view --targets-file ^remove_hetsites.txt "$PREFIX"_yesALT_snps.vcf.gz -o "$PREFIX"_temp3.vcf

# log
echo "site heterozygisty filter complete"

echo "Heterozygosity site filter" >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_temp3.vcf | wc -l) \
>> "$PREFIX"_log.txt

#### record allele frequncies
plink2 --vcf "$PREFIX"_temp3.vcf \
  --freq cols=chrom,pos,ref,alt,reffreq,altfreq,nobs \
  --allow-extra-chr \
  --double-id \
  --out "$PREFIX"

## Calculate minor allele frequency
# note that script (of course) only considers variant sites for MAF caluclation
# generates list of sites to remove because they do not reach MAF threshold (remove_MAF.txt) and 
# generates list of sites with known ALT allele but are invariant (invariant_sites.txt)
Rscript AllSites_CO_calcMAF.R $WORKDIR $PREFIX

## append invariant sites to low MAF sites file
cat invariant_sites.txt >> remove_MAF.txt

# Use bcftools to filter on MAF
bcftools view --targets-file ^remove_MAF.txt "$WORKDIR"/"$PREFIX"_temp3.vcf \
	-Oz -o "$WORKDIR"/"$PREFIX"_maf_snps.vcf.gz
# index for concaetnation later
bcftools index "$WORKDIR"/"$PREFIX"_maf_snps.vcf.gz

# log
echo "minor allele frquency filter complete"
echo "Variant site MAF filter" >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_maf_snps.vcf.gz | wc -l) \
>> "$PREFIX"_log.txt

####################### Generate invariant sites VCF
# select sites that have a known ALT allele but allele freq=1
bcftools view --targets-file invariant_sites.txt "$WORKDIR"/"$PREFIX"_temp3.vcf \
        -Oz -o "$WORKDIR"/"$PREFIX"_temp_invariant.vcf.gz
# index for concatentation
bcftools index "$WORKDIR"/"$PREFIX"_temp_invariant.vcf.gz

# Join invariant sites (with and without known alt allele) with variant sites
bcftools concat --allow-overlaps --rm-dups exact \
	"$WORKDIR"/"$PREFIX"_temp_invariant.vcf.gz \
	"$WORKDIR"/"$PREFIX"_noALT.vcf.gz \
	"$WORKDIR"/"$PREFIX"_maf_snps.vcf.gz \
	-Ou | bcftools sort -Ov -o "$WORKDIR"/"$PREFIX"_AllSites_snps.vcf

# log
echo "Invariant and Variant site combination complete"
echo "Invariant and filtered variant sites combined" >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_AllSites_snps.vcf | wc -l) \
>> "$PREFIX"_log.txt


# rename samples
bcftools query -l "$PREFIX"_AllSites_snps.vcf > vcf_samples.txt
sed 's/_R1_001.fastq.gz.trimmed.fastq.sam//g' vcf_samples.txt > temp_names.txt
sed -i 's/_R1.fastq.gz.sam//g' temp_names.txt
paste vcf_samples.txt temp_names.txt > vcf_new_names.txt

bcftools reheader --samples vcf_new_names.txt "$PREFIX"_AllSites_snps.vcf | bcftools view --threads 30 -Oz -o "$PREFIX"_invariantWithfiltered_snps_rehead.vcf.gz
bcftools reheader --samples vcf_new_names.txt "$PREFIX"_maf_snps.vcf.gz | bcftools view --threads 30 -Oz -o "$PREFIX"_filtered_snps_rehead.vcf.gz

### Individual stats
plink2 --vcf "$PREFIX"_filtered_snps_rehead.vcf.gz \
	--sample-counts cols=homref,homalt,het,hapref,hapalt \
	--missing \
	--allow-extra-chr \
	--double-id \
	--out "$PREFIX"
