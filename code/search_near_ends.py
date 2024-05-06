from Bio.Seq import Seq
import sys


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def is_sequence_near_ends(seq, target_sequence, max_distance):
    return (seq.startswith(target_sequence) or
            seq.endswith(target_sequence) or
            seq[:max_distance].find(target_sequence) != -1 or
            seq[-max_distance:].find(target_sequence) != -1)


def process_fastq(file_path, sequence, output_path, max_distance):
    reads_found = 0
    total_reads = 0

    with open(file_path, 'r') as fastq, open(output_path, 'w') as output:
        while True:
            header = fastq.readline().strip()
            if not header:
                break
            seq = fastq.readline().strip()
            plus = fastq.readline().strip()
            quality = fastq.readline().strip()

            total_reads += 1  # Increment total_reads here for each complete read

            if sequence in seq or sequence in reverse_complement(seq):
                reads_found += 1
                output.write(f'{header}\n{seq}\n{plus}\n{quality}\n')

    print(f"Total Reads: {total_reads}, Reads Found: {reads_found}")


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python script.py <FASTQ file> <DNA sequence> <Output file> <Max distance>")
        sys.exit(1)

    fastq_file = sys.argv[1]
    dna_sequence = sys.argv[2].upper()
    output_file = sys.argv[3]
    max_distance = int(sys.argv[4])

    process_fastq(fastq_file, dna_sequence, output_file, max_distance)
