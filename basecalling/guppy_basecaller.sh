#!/bin/bash -l

#SBATCH --job-name=TPC240218_GupRACE2400bps
#SBATCH --output=/scratch/prj/ont_zap_seq/Thomas/%u/%j.out
#SBATCH --error=/users/k2149964/error_logs/TPC240218_GupRACE2400bps_%j.err
#SBATCH --gres=gpu

# Load modules or set environment variables
module load test_switch_kcl
source test_switch
module load ont-guppy/6.4.6-cuda-gcc-10.3.0
module load cuda

# Execute your application or script
time srun guppy_basecaller -i /scratch/prj/ont_zap_seq/Thomas/20240205_1352_MN41657_FAY64429_5b6f8c7c/ \
    -s /scratch/prj/ont_zap_seq/Thomas/240218_Guppy400bps/ \
    -c dna_r10.4.1_e8.2_400bps_sup.cfg \
    --barcode_kits "SQK-NBD114-24" \
    --require_barcodes_both_ends \
    --detect_mid_strand_barcodes \
    --detect_mid_strand_adapter \
    --compress_fastq \
    --recursive \
    --do_read_splitting \
    --enable_trim_barcodes \
    --device auto
