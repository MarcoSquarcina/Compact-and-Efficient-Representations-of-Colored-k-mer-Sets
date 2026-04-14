#!/usr/bin/env python3

import sys

if __name__ == '__main__':
    fasta_filename = "test-merged.fa"
    sequence_filename = "test.fa"
    colors_filename = "colors.rle.colors"
    k = 3

    if len(sys.argv) == 5:
        sequence_filename = sys.argv[1]
        colors_filename = sys.argv[2]
        fasta_filename = sys.argv[3]
        k = int(sys.argv[4])
    else:
        print("Usage: python fasta-merger <sequences.fa> <colors.rle> <output.fa> <k>")
        print("Going with defaults...")

    # collecting colors vectors
    print(f"Reading color file \"{colors_filename}\"...")
    colors = []
    runs = []
    try:
        with open(colors_filename, 'r') as colors_rle:
            for line in colors_rle:
                rle = line.split(sep=":")
                colors.append(int(rle[0]))
                if len(rle) == 2:
                    runs.append(int(rle[1]))
                else:
                    runs.append(1)

        print("Merging sequences and colors...")
        partial_sum_kmer = 0
        partial_sum_color = 0
        run_id = 0
        with open(fasta_filename, 'w') as fasta:
            with open(sequence_filename, 'r') as sequences:
                for line in sequences:
                    if line[0] == '>':  # found a header
                        continue  # need to get sequence information
                    else:  # found the sequence
                        sequence = line.strip()
                        sequence_length = len(sequence)
                        kmer_num = sequence_length - k + 1
                        partial_sum_kmer += kmer_num
                        header = ">"
                        while partial_sum_kmer > partial_sum_color:
                            if len(header) == 1:
                                header += f"{colors[run_id]}:{runs[run_id]}"
                            else:
                                header += f":{colors[run_id]}:{runs[run_id]}"
                            partial_sum_color += runs[run_id]
                            run_id += 1
                        fasta.write(f"{header}\n{sequence}\n")
    except FileNotFoundError:
        print("Can't open file!")
        exit(1)

    assert partial_sum_color == partial_sum_kmer
    print(f"Total number of runs: {run_id}")
    print(f"Total number of kmers: {partial_sum_kmer}")
    print(f"Merged file saved in {fasta_filename}.")
