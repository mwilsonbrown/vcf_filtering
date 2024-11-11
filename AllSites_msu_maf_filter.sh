#!/bin/bash
#
#SBATCH --job-name=mafAllSites
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-8:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/slurm-%A.out
# Filter allSites VCF for minor allele frequency in the variant sites only (MAF)
# November 11, 2024
# Maya Wilson Brown
#
#
#
# change directory
cd /mnt/scratch/wils1582/allSites_msu_filtering/

### Varibles
VCF=/mnt/research/josephslab/Maya/CBP_AllSites_msu_qual_filtered.vcf.gz
PREFIX=CBP_AllSites_msu_maf


# purge modules
module purge
# load modules
module load  BCFtools/1.19-GCC-13.2.0

# Generate log file
touch ./"$PREFIX"_log.txt

# Select invariant sites that are all missing or reference
bcftools view -i 'INFO/AC == 0' $VCF -Ou | bcftools sort - -Ou -o cbp_allsites_invariant_ref_miss.bcf

# check the number of variants (I expect ~31k based on troubleshooting)
echo "Invariant sites that are all reference or missing calls" >> "$PREFIX"_log.txt
bcftools query -f '%CHROM\t%POS\n' cbp_allsites_invariant_ref_miss.bcf | wc -l >> "$PREFIX"_log.txt

# select variant sites
bcftools view -i 'INFO/AC > 0' $VCF -Ou -o cbp_allsites_variant.bcf

# log number of variant sites
echo "Variant sites" >> "$PREFIX"_log.txt
bcftools query -f '%CHROM\t%POS\n' cbp_allsites_variant.bcf | wc -l >> "$PREFIX"_log.txt

# Filter for MAF on vairant sites only
bcftools view -e 'MAF < 0.01' cbp_allsites_variant.bcf -Ou | bcftools sort - -Ou -o cbp_allsites_variant_maf_filt.bcf

# log number of variant sites after filtering
echo "Variant sites after MAF filter" >> "$PREFIX"_log.txt
bcftools query -f '%CHROM\t%POS\n' cbp_allsites_variant_maf_filt.bcf | wc -l >> "$PREFIX"_log.txt

# merge the invariant and variant sites
bcftools concat --allow-overlaps --write-index \
  cbp_allsites_variant_maf_filt.bcf cbp_allsites_invariant_ref_miss.bcf \
  -Oz -o "$PREFIX".vcf.gz 

# check total number of sites after concatenating
bcftools query -f '%CHROM\t%POS\n' "$PREFIX".vcf.gz | wc -l >> "$PREFIX"_log.txt


