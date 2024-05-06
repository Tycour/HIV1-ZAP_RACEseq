# HIV1-ZAP_RACEseq
Bioinformatics pipeline for transcript cleavage analysis based on 5'RACE NGS/Nanopore data

## Installation
Create a conda environment and required packages

```bash
$ conda env create -f 5pRACE_Analysis.yml
$ conda activate 5pRACE_Analysis
```

## Usage Example

The default structure of the [5pRACE_pipeline.sh](code%2F5pRACE_pipeline.sh) pipeline assumes that basecalled Nanopore reads have been demultiplexed with barcodes trimmed, and generated .fastq.gz files are separated into different folders.

>**The many paths included in the scripts are personal and must be amended by the user.**

Change your working directory to a folder containing the *.fastq.gz* files for a given barcode. Here is an example usage of the pipeline for that same demultiplexed barcode/sample:

```bash
$ ./5pRACE_pipeline.sh -P "your_chosen_prefix" -q 20 -v WT -t 4
```
