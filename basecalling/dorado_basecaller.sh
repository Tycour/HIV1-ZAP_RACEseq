#!/bin/bash -l

#SBATCH --job-name=TPC240218_RACE2_dor
#SBATCH --output=/scratch/prj/ont_zap_seq/Thomas/%u/%j.out
#SBATCH --error=/users/k2149964/error_logs/TPC240218_RACE2_dor_%j.err
#SBATCH --gres=gpu

# Load modules or set environment variables
module load test_switch_kcl
source test_switch
module load dorado/0.5.1-gcc-10.3.0
module load cuda

nvidia-debugdump -l

srun dorado basecaller \
    /users/k2149964/basecalling/dorado/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
    /scratch/prj/ont_zap_seq/Thomas/20240205_1352_MN41657_FAY64429_5b6f8c7c \
    --kit-name SQK-NBD114-24 \
    --barcode-both-ends \
    --recursive \
    > /scratch/prj/ont_zap_seq/Thomas/dorado_called/240218_RACE2.bam

srun dorado demux \
    --emit-fastq \
    --output-dir /scratch/prj/ont_zap_seq/Thomas/dorado_called \
    --kit-name SQK-NBD114-24 \
    /scratch/prj/ont_zap_seq/Thomas/dorado_called/240218_RACE2.bam
