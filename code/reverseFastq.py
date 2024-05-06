import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to reverse complement a sequence and invert quality scores
def reverse_complement_and_invert_quality(record):
    seq = record.seq.reverse_complement()
    qual = record.letter_annotations["phred_quality"][::-1]
    inverted_record = SeqRecord(seq, id=record.id, description=record.description)
    inverted_record.letter_annotations["phred_quality"] = qual
    return inverted_record

# Function to remove a specified sequence from the read and quality scores
def remove_sequence(record, sequence_to_remove):
    seq = record.seq
    qual = record.letter_annotations["phred_quality"]
    start_index = seq.find(sequence_to_remove)
    if start_index != -1:
        end_index = start_index + len(sequence_to_remove)
        seq = seq[end_index:]
        qual = qual[end_index:]
    return SeqRecord(seq, id=record.id, description=record.description, letter_annotations={"phred_quality": qual})

def main():
    parser = argparse.ArgumentParser(description="Reverse complement sequences, invert quality scores, and remove specified sequences in a FASTQ file.")
    parser.add_argument("input_file", help="Input FASTQ file")
    parser.add_argument("output_file", help="Output FASTQ file")
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    sequence_to_remove = "CTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGTACATGGG"
    final_sequence_to_remove = "CTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGTAC"

    # Process the input FASTQ file and create the output FASTQ file
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fastq"):
            inverted_record = reverse_complement_and_invert_quality(record)
            final_record = remove_sequence(inverted_record, sequence_to_remove)
            if final_sequence_to_remove not in final_record.seq:
                SeqIO.write(final_record, output_handle, "fastq")

    print("Processing completed. Inverted sequences with specified removal and final check written to", output_file)

if __name__ == "__main__":
    main()
