import argparse
import re
from Bio import SeqIO

def is_valid_upstream_sequence(sequence, pattern):
    """Check if the upstream sequence matches the degenerate pattern"""
    return re.match(pattern, sequence) is not None

def trim_fastq_sequences(input_fastq, output_fastq, sequence_to_find, degenerate_pattern):
    pattern_regex = re.compile(degenerate_pattern)

    with open(input_fastq, "r") as input_file, open(output_fastq, "w") as output_file:
        for record in SeqIO.parse(input_file, "fastq"):
            sequence = str(record.seq)
            sequence_pos = sequence.find(sequence_to_find)

            if sequence_pos >= 15:  # Check if there are at least 15 bases upstream
                upstream_sequence = sequence[sequence_pos-15:sequence_pos]

                # Check if the upstream sequence matches the degenerate pattern
                if is_valid_upstream_sequence(upstream_sequence, pattern_regex):
                    # Trim the sequence and quality scores
                    record = record[sequence_pos-15:]
                    SeqIO.write(record, output_file, "fastq")

def main():
    parser = argparse.ArgumentParser(description='Process a FASTQ file by trimming reads based on a specified sequence and checking for a degenerate sequence pattern upstream.')
    parser.add_argument('input_fastq', help='Path to the input FASTQ file.')
    parser.add_argument('output_fastq', help='Path to the output FASTQ file.')
    parser.add_argument('sequence', help='Sequence to find in the reads.')
    parser.add_argument('degenerate_pattern', help='Degenerate sequence pattern to check upstream of the found sequence.')

    args = parser.parse_args()

    # Convert the degenerate sequence to a regex pattern
    degenerate_regex_pattern = args.degenerate_pattern.replace('N', '.').replace('W', '[AT]')

    trim_fastq_sequences(args.input_fastq, args.output_fastq, args.sequence, degenerate_regex_pattern)

if __name__ == "__main__":
    main()
