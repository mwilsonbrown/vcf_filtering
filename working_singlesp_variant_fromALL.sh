#!/bin/bash --login
#
#SBATCH --job-name=mafCBPvariant
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/%A.out
# Working script 
# To go from singple species AllSites base filters (ends with filter3)
# thorugh more site level filters and filtering for minor allele frquency
# The resulting the input and output VCFs from this script are both ALLSITES VCFs
#
# December 1, 2025
# Maya Wilson Brown
#
# Updated: January 15, 2026
######################## SETUP
# purge modules
module purge
# load modules
module load R/4.3.2-gfbf-2023a
module load PLINK/2.00a3.7-gfbf-2023a
module load BCFtools/1.18-GCC-12.3.0 #this version of bcftools works with the other modules
## Other vars
PREFIX=CBP_allSites_msu225
WORKDIR=/mnt/scratch/wils1582/allSites_msu_filtering

# In the case that I do not have execute permissions for my own github repo,
# copy those files to the working directory
cd "$WORKDIR"

####### split invariant and potentially variant sites
# separate out invariant sites with unseen alternative allele
bcftools view -i 'ALT="."' "$PREFIX"_filter3.vcf.gz -Oz -o "$PREFIX"_noALT.vcf.gz
# index for concatentation
bcftools index "$PREFIX"_noALT.vcf.gz
# do the oppposite and keep sites with seen alleles of varying allele frequency
bcftools view -e 'ALT="."' -v snps "$PREFIX"_filter3.vcf.gz -Oz -o "$PREFIX"_yesALT_snps.vcf.gz

#### record allele frequncies
plink2 --vcf "$PREFIX"_yesALT_snps.vcf.gz \
  --freq cols=chrom,pos,ref,alt,reffreq,altfreq,nobs \
  --allow-extra-chr \
  --double-id \
  --out "$PREFIX"

## Calculate minor allele frequency
# note that script (of course) only considers variant sites for MAF caluclation
# generates list of sites to remove because they do not reach MAF threshold and 
# generates list of sites with known ALT allele but are invariant
Rscript AllSites_CO_calcMAF.R $WORKDIR $PREFIX

# Use bcftools to filter on MAF
bcftools view --targets-file ^remove_MAF.txt "$WORKDIR"/"$PREFIX"_yesALT_snps.vcf.gz \
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
bcftools view --targets-file invariant_sites.txt "$WORKDIR"/"$PREFIX"_yesALT_snps.vcf.gz \
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

echo "Invariant sites combined" >> "$PREFIX"_log.txt
echo $(bcftools query -f'%CHROM %POS\n' "$PREFIX"_AllSites_snps.vcf | wc -l) \
>> "$PREFIX"_log.txt


# rename samples
bcftools query -l "$PREFIX"_AllSites_snps.vcf > vcf_samples.txt
sed 's/_R1_001.fastq.gz.trimmed.fastq.sam//g' vcf_samples.txt > temp_names.txt
sed -i 's/_R1.fastq.gz.sam//g' temp_names.txt
paste vcf_samples.txt temp_names.txt > vcf_new_names.txt
bcftools reheader --samples vcf_new_names.txt "$PREFIX"_AllSites_snps.vcf | bcftools view --threads 30 -Oz -o "$PREFIX"_AllSites_snps_rehead.vcf.gz
