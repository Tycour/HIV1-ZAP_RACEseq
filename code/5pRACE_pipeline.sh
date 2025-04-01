#!/bin/bash

### Define Parameters and default paths
PREFIX="WC"
QSCORE=20
VIRUS="WT"
THREADS=16
DISTANCE=150
HAMMINGDIST=2

# Create necessary directories
mkdir -p reads nanoplot alignments

# Override default arguments with command-line arguments if provided
while getopts "P:q:v:t:r:h:" opt; do
  case $opt in
    P) PREFIX="$OPTARG" ;;
    q) QSCORE="$OPTARG" ;;
    v) VIRUS="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    r) REF_GENOME="$OPTARG" ;;
    h) HAMMINGDIST="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2
        exit 1
    ;;
  esac
done

GENOME_FASTA="/Users/tycour/PycharmProjects/HIV1/venv/data/references/HIV-${VIRUS}.fasta"
REF_GENOME="/Users/tycour/PycharmProjects/HIV1/venv/data/references/HIV-${VIRUS}.mmi"

function merge_reads() {
    mkdir -p reads
    cat *.fastq.gz > reads/"${PREFIX}".fastq.gz
    rm *.fastq.gz
    gunzip -c reads/"${PREFIX}".fastq.gz > reads/"${PREFIX}".fastq
}

function filter_reads() {
    nanofilt -q "${QSCORE}" reads/"${PREFIX}.fastq" > reads/"${PREFIX}_q${QSCORE}.fastq"
}

function quality_control() {
    nanoplot -t "$THREADS" --fastq reads/"${PREFIX}.fastq" -o nanoplot/unfiltered
    nanoplot -t "$THREADS" --fastq reads/"${PREFIX}_q${QSCORE}.fastq" -o nanoplot/q"${QSCORE}"_filtered
}

function process_cDNA_reads() {
    python /Users/tycour/PycharmProjects/HIV1/venv/data/code/search_near_ends.py \
    reads/"${PREFIX}"_q"${QSCORE}".fastq \
    CTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGTACATGGG \
    reads/"${PREFIX}"_q"${QSCORE}"_5p.fastq \
    ${DISTANCE}

    python /Users/tycour/PycharmProjects/HIV1/venv/data/code/search_near_ends.py \
    reads/"${PREFIX}"_q"${QSCORE}"_5p.fastq \
    CCAGCAGCAGATGGGGTGGGAGCAGTAT \
    reads/"${PREFIX}"_q"${QSCORE}"_5p3p.fastq \
    ${DISTANCE}
}

function align_reads() {
    minimap2 -ax splice -t "$THREADS" -uf "$REF_GENOME" reads/"${PREFIX}"_q"${QSCORE}"_5p3p.fastq > alignments/"${PREFIX}"_q"${QSCORE}"_aln.sam
    check_exit_status "Alignment" $?
}

function convert_sort_index_bam() {
    samtools view -@ "$THREADS" -bS alignments/"${PREFIX}"_q"${QSCORE}"_aln.sam | samtools sort -@ "$THREADS" -o alignments/"${PREFIX}"_q"${QSCORE}"_aln_sorted.bam
    samtools index -@ "$THREADS" alignments/"${PREFIX}"_q"${QSCORE}"_aln_sorted.bam
    check_exit_status "Conversion to BAM and sorting" $?
}

function generate_alignment_stats() {
    samtools flagstat alignments/"${PREFIX}"_q"${QSCORE}"_aln_sorted.bam > alignments/"${PREFIX}"_q"${QSCORE}"_alignment_stats.txt
    check_exit_status "Generation of alignment statistics" $?
}

function check_exit_status() {
    local command_name="$1"
    local status="$2"
    if [ "$status" -ne 0 ]; then
        echo "$command_name failed."
        exit 1
    fi
    echo "$command_name completed successfully."
}

function orientReads() {
    python /Users/tycour/PycharmProjects/HIV1/venv/data/code/orientRead.py \
    alignments/"${PREFIX}"_q"${QSCORE}"_aln.sam \
    reads/"${PREFIX}"_q"${QSCORE}"_5p3p.fastq \
    reads/"${PREFIX}"_q"${QSCORE}"_rev.fastq
}

function posReads() {
    python /Users/tycour/PycharmProjects/HIV1/venv/data/code/posRead.py \
    reads/"${PREFIX}"_q"${QSCORE}"_rev.fastq \
    reads/"${PREFIX}"_q"${QSCORE}"_revFiltered.fastq \
    ATACTGCTCCCACCCCATCTGCTGCTGG \
    NNNWNNNWNNNWNNN
}

function extractUMIs() {
    umi_tools extract \
    --stdin=reads/"${PREFIX}"_q"${QSCORE}"_revFiltered.fastq \
    --bc-pattern=NNNNNNNNNNNNNNN \
    --log=processed.log \
    --stdout reads/"${PREFIX}"_q"${QSCORE}"_revUMI.fastq
}

function revFastq() {
    python /Users/tycour/PycharmProjects/HIV1/venv/data/code/reverseFastq.py \
    reads/"${PREFIX}"_q"${QSCORE}"_revUMI.fastq \
    reads/"${PREFIX}"_q"${QSCORE}"_forUMI.fastq
}

function createBAM() {
    minimap2 -ax splice -uf -t "$THREADS" "${REF_GENOME}" reads/"${PREFIX}"_q"${QSCORE}"_forUMI.fastq > alignments/"${PREFIX}"_processed_aln.sam
    samtools view -bS -F 4 alignments/"${PREFIX}"_processed_aln.sam | samtools sort -o alignments/"${PREFIX}"_processed_aln_sorted.bam # Added -F 4 to remove unmapped reads
}

#function UMIcollapse() {
#    umicollapse bam -i alignments/"${PREFIX}"_processed_aln_sorted.bam \
#    -k "$HAMMINGDIST" \
#    -o alignments/"$PREFIX"_dedup.bam
#}

# Change _dedup.bam to _processed_aln_sorted.bam
function findStarts() {
    bedtools bamtobed -i alignments/"$PREFIX"_processed_aln_sorted.bam > alignments/"$PREFIX"_processed_aln_sorted.bed
    awk 'BEGIN {OFS="\t"} {if ($6 == "+") print $1, $2, $2+1; else print $1, $3-1, $3}' alignments/"$PREFIX"_processed_aln_sorted.bed > alignments/"$PREFIX"_5start.bed
    python /Users/tycour/PycharmProjects/HIV1/venv/data/code/countPlotStarts.py \
    alignments/"$PREFIX"_5start.bed \
    "${GENOME_FASTA}" \
    "$PREFIX"_final.csv \
    "$PREFIX"_plot.png \
    --chop_bases 0
}

# Main script execution
merge_reads
filter_reads
#quality_control ## Optional step to QC sequencing results
process_cDNA_reads
align_reads
convert_sort_index_bam
generate_alignment_stats
orientReads
posReads
extractUMIs
revFastq
createBAM
#UMIcollapse ## Optional
findStarts

echo "Script completed. Check output files for results."
