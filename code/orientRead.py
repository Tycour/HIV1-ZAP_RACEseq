import argparse
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def is_forward_strand(flag):
    """ Check if the read is aligned to the forward strand """
    return (flag & 16) == 0


def reverse_complement_read(seq_record):
    """ Reverse complement the sequence and reverse the quality scores of a SeqRecord """
    seq_record.seq = seq_record.seq.reverse_complement()  # Reverse complement the sequence
    seq_record.letter_annotations["phred_quality"] = seq_record.letter_annotations["phred_quality"][::-1]  # Reverse the quality scores
    return seq_record


def process_files(samfile_path, fastq_path, output_fastq_path):
    # Read the alignment file and store the orientation of each read
    orientation_dict = {}
    with pysam.AlignmentFile(samfile_path, "r") as samfile:
        for read in samfile:
            read_name = read.query_name
            orientation_dict[read_name] = is_forward_strand(read.flag)

    # Process the FASTQ file
    with open(fastq_path, "r") as fastq_file, open(output_fastq_path, "w") as output_file:
        for record in SeqIO.parse(fastq_file, "fastq"):
            if orientation_dict.get(record.id, False):
                # Reverse complement the read if it's aligned to the forward strand
                record = reverse_complement_read(record)
            SeqIO.write(record, output_file, "fastq")


def main():
    parser = argparse.ArgumentParser(
        description='Reverse complement sequences in a FASTQ file based on SAM/BAM alignment.')
    parser.add_argument('samfile', help='Path to the SAM/BAM alignment file.')
    parser.add_argument('fastq', help='Path to the original FASTQ file.')
    parser.add_argument('output', help='Path to the output FASTQ file.')

    args = parser.parse_args()

    process_files(args.samfile, args.fastq, args.output)


if __name__ == "__main__":
    main()
