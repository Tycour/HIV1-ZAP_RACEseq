import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from Bio import SeqIO
import numpy as np

def read_first_reference_genome(fasta_file):
    first_record = next(SeqIO.parse(fasta_file, "fasta"), None)
    if first_record:
        return first_record.id, str(first_record.seq)
    raise ValueError("No sequences found in the FASTA file")

def main(args):
    # Read the first reference genome sequence
    reference_name, reference_sequence = read_first_reference_genome(args.fasta_file)

    # Load the BED data
    bed_data = pd.read_csv(args.bed_file, sep='\t', header=None, names=['chrom', 'start', 'end'])

    # Filter for the specific reference and adjust for 1-based index
    filtered_data = bed_data[bed_data['chrom'] == reference_name]
    filtered_data['start'] += 1

    # Count the number of reads starting at each position
    position_counts = filtered_data['start'].value_counts().sort_index()

    # Add 0 values for missing positions
    max_position = len(reference_sequence)
    all_positions = pd.Series(index=np.arange(1, max_position + 1), data=0).rename('count')
    position_counts = position_counts.add(all_positions, fill_value=0).sort_index()

    # Get the base associated with each position
    position_counts = position_counts.reset_index()
    position_counts.columns = ['position', 'count']
    position_counts['base'] = position_counts['position'].apply(lambda x: reference_sequence[x - 1] if x <= len(reference_sequence) else 'N')

    # Filter out the first X bases as specified in args.chop_bases
    position_counts_filtered = position_counts[position_counts['position'] > args.chop_bases]

    # Normalise the counts to percentages of the total reads
    total_reads = position_counts_filtered['count'].sum()
    position_counts_filtered['percentage'] = (position_counts_filtered['count'] / total_reads) * 100

    # Save the counts, bases, and normalized percentages to a file
    position_counts_filtered.to_csv(args.output_csv, index=False)

    # Plotting the data
    plt.figure(figsize=(12, 6))
    plt.plot(position_counts_filtered['position'], position_counts_filtered['percentage'])

    # Setting a logarithmic scale for y-axis with raw number labels
    plt.yscale('log')
    ax = plt.gca()  # Get current axis
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))

    # Enhancing the plot aesthetics
    plt.title(
        f'Read Start Positions as Percentage Along {reference_name}',
        fontsize=14, fontweight='bold')
    plt.xlabel('Position', fontsize=12)
    plt.ylabel('% Total Reads', fontsize=12)
    plt.grid(True, which="both", linestyle='--', linewidth=0.5)
    plt.tight_layout()

    # Save and show the plot
    plt.savefig(args.output_plot)
    plt.show()

    print(f"Counts of reads starting at each position with associated bases and normalized percentages have been saved to {args.output_csv}")
    print(f"Plot has been saved as '{args.output_plot}'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count and plot read start positions with associated bases from the first sequence in a reference genome, normalized by total read count")
    parser.add_argument("bed_file", help="Path to the BED file")
    parser.add_argument("fasta_file", help="Path to the reference genome FASTA file")
    parser.add_argument("output_csv", help="Output CSV file name")
    parser.add_argument("output_plot", help="Output plot file name")
    parser.add_argument("--chop_bases", type=int, default=0, help="Number of bases to exclude from the beginning of the plot (default 0)")

    args = parser.parse_args()
    main(args)
