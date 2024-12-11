#!/bin/bash --login
#SBATCH --job-name=ReLERNN
#SBATCH --constraint=amd20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wils1582@msu.edu
#SBATCH --output=/mnt/scratch/wils1582/slurm/%x-%A.SLURMout
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
# ReLERNN on VCF
# DEC 2, 2024
# Maya Wilson Brown
# constrained job request to amd20 nodes because that is the node I built tensorflow on
# output information about how this job is running using bash commands
echo "This job is running on $HOSTNAME on `date`"


######## SETUP
export PATH=/mnt/home/wils1582/miniconda3/bin:$PATH
source /mnt/home/wils1582/miniconda3/bin/activate relernn

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lib/:/lib64/:$CONDA_PREFIX/lib/:$CONDA_PREFIX/lib/python3.9/site-packages/tensorrt
export XLA_FLAGS=--xla_gpu_cuda_data_dir=$CONDA_PREFIX/lib

cd /mnt/home/wils1582/ReLERNN/examples/

######## PIPELINE
SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="42"
MU="1e-8"
URTR="1"
DIR="./example_output/"
VCF="./example.vcf"
GENOME="./genome.bed"
MASK="./accessibility_mask.bed"

# Simulate data
${SIMULATE} \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --mask ${MASK} \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --upperRhoThetaRatio ${URTR} \
    --nTrain 13000 \
    --nVali 2000 \
    --nTest 100 \
    --seed ${SEED}

# Train network
${TRAIN} \
    --projectDir ${DIR} \
    --seed ${SEED}

# Predict
${PREDICT} \
    --vcf ${VCF} \
    --projectDir ${DIR} \
    --seed ${SEED}

# Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${DIR} \
    --nSlice 2 \
    --nReps 2 \
    --seed ${SEED}
