# just  noodling around to mkae LD decay plots so I can choose biologically informed window sizes when doing windowed  genomic analyses
# M Wilson Brown
# September 11, 2024

##modules
#module purge
#module load PLINK/2.00a3.7-gfbf-2023a
#
## setting up VCF variables 
VCF=/mnt/home/wils1582/allSites_CBP_final.filtered.vcf.gz
## run PLINK on Capsella bursa-pastoris vcf only on one chromosome
#plink2 --vcf $VCF --double-id --allow-extra-chr \
#	--set-missing-var-ids @:# \
#	--chr 'jlSCF_10' \
#	--thin 0.1 \
#	--r2-unphased gz \
#	--ld-window 100 \
#	--ld-window-kb 1000 \
#	--ld-window-r2 0 \
#	--make-bed --out cbp_LD_JLv4
## The above does not recognizze its own flags and throws an error
# advice online is to just use PLINK1.9 when this happens so I guess I will do that (silly)
module purge
module load PLINK/1.9b_6.21-x86_64 

plink --vcf $VCF --double-id --allow-extra-chr \
        --set-missing-var-ids @:# \
        --thin 0.1 -r2 gz \
        --ld-window 100 \
        --ld-window-kb 1000 \
        --ld-window-r2 0 \
        --make-bed --out cbp_LD_JLv4_all
