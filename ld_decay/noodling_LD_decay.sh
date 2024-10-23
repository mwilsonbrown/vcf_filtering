# Making LD decay plots so I can choose biologically informed window sizes when doing windowed  genomic analyses
# M Wilson Brown
# October 16, 2024

# move to scratch
cd /mnt/scratch/wils1582/

# initialize conda environment for python2
source /mnt/home/wils1582/miniconda3/bin/activate py2

module purge
module load PLINK/1.9b_6.21-x86_64 


# Originally, I was doing this analysis on NYC and East Asian Cbp separately, based on previously known populations.
# Here, I do not have previously known populations so I will use all the CBP samples together.
# I will compare the all sites VCF and variant-only vcf but use all the cbp populations to estimate
# as LD is expected to be stronger in local populations (https://www.nature.com/articles/ng813z)
#
# The vcf is already filtered in some ways, but it will not hurt to add additional filters to the PLINK call
# Here, I thin the data to a random subset of 10% of SNPs, calculate the r^2 correlation coefficient,
# gzip the output
# I set the lower window boundary to 100 sites (ie ignore anythig closer than 100 sites), 
# and the upper window bounrary to 1000kb (1 Mb)
#
##### Global sample of cbp from variant only VCF
VCF=/mnt/research/josephslab/Maya/CBP_CRCGCONP_filtered.vcf.gz
# take the first field of the pop file (cbp names), and duplicate it via awk for PLINK
POPFILE=$(cut -f1 /mnt/home/wils1582/capsella_population_structure/cbp_pop_str.txt | | awk -F '\t' '{print $1, $NF}')
$PREFIX=var_cbp_JLv4_LD

# Use PLINK to calculate r^2 between SNPs
plink --vcf $VCF --double-id \
        --keep $POPFILE \
        --allow-extra-chr \
        --set-missing-var-ids @:# \
        --maf 0.025 \
        --geno 0.05 \
        --mind 0.5 \
        --thin 0.1 -r2 gz \
        --ld-window 100 \
        --ld-window-kb 1000 \
        --ld-window-r2 0 \
        --make-bed \
        --out $PREFIX
        
# Then, calculate average LD in genomic bins; written in python 2 so need to activate that conda environment
python ld_decay_cal.py -i "$PREFIX".ld.gz -o "$PREFIX"_chrAll

# 
# 
# module purge
# module load PLINK/1.9b_6.21-x86_64 
# 
# plink --vcf $VCF --double-id \
#         --keep /mnt/home/wils1582/vcf_filtering/ld_decay/easia_vcf_sample2.txt \
#         --allow-extra-chr \
#         --set-missing-var-ids @:# \
#         --maf 0.025 \
#         --geno 0.05 \
#         --mind 0.5 \
#         --thin 0.1 -r2 gz \
#         --ld-window 100 \
#         --ld-window-kb 1000 \
#         --ld-window-r2 0 \
#         --make-bed \
#         --out cbp_LD_JLv4_easia
        
# Then, calculate average LD in genomic bins; written in python 2 so need to activate that conda environment
#python ld_decay_cal.py -i cbp_LD_JLv4_easia.ld.gz -o easia_chrAll
