# Making LD decay plots so I can choose biologically informed window sizes when doing windowed  genomic analyses
# M Wilson Brown
# October 16, 2024

##modules
#module purge
#module load PLINK/2.00a3.7-gfbf-2023a
#
## setting up VCF variables 
#VCF=/mnt/home/wils1582/allSites_CBP_final.filtered.vcf.gz
VCF=/mnt/scratch/wils1582/july15/CBP_CRCG_final_filtered.vcf.gz

# I plan on calculating LD on the different subpopulations separately, so I need the population file for that
# Load population file
#POPFILE=/mnt/home/wils1582/capsella_population_structure/cbp_pop_str.txt

# use awk to split them by the third column
#awk '$3 == N_Europe' $POPFILE > n_europe

# The vcf is already filtered in some ways, but it will not hurt to add additional filters to the PLINK call
# Here, I thin the data to a random subset of 10% of SNPs, calculate the r^2 correlation coefficient,
# gzip the output
# I set the lower window boundary to 100 sites (ie ignore anythig closer than 100 sites), 
# and the upper window bounrary to 1000kb (1 Mb)


module purge
module load PLINK/1.9b_6.21-x86_64 

plink --vcf $VCF --double-id \
        --keep /mnt/home/wils1582/vcf_filtering/ld_decay/easia_vcf_sample2.txt \
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
        --out cbp_LD_JLv4_easia
        
# Then, calculate average LD in genomic bins; written in python 2 so need to activate that conda environment
python ld_decay_cal.py -i cbp_LD_JLv4_easia.ld.gz -o easia_chrAll
