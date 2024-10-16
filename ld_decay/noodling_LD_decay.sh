# Making LD decay plots so I can choose biologically informed window sizes when doing windowed  genomic analyses
# M Wilson Brown
# October 16, 2024

##modules
#module purge
#module load PLINK/2.00a3.7-gfbf-2023a
#
## setting up VCF variables 
#VCF=/mnt/home/wils1582/allSites_CBP_final.filtered.vcf.gz

# I plan on calculating LD on the different subpopulations separately, so I need the population file for that
# Load population file
POPFILE=/mnt/home/wils1582/capsella_population_structure/cbp_pop_str.txt

# use awk to split them by the third column
awk '$3 == N_Europe' $POPFILE > n_europe


module purge
module load PLINK/1.9b_6.21-x86_64 

plink --vcf $VCF --double-id \
        --allow-extra-chr \
        --set-missing-var-ids @:# \
        --thin 0.1 -r2 gz \
        --ld-window 100 \
        --ld-window-kb 1000 \
        --ld-window-r2 0 \
        --make-bed --out cbp_LD_JLv4_all
