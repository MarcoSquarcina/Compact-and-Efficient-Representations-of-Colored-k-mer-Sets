#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt


def save_to_csv(data, filename):
    with open(filename, 'w') as csv_file:
        csv_file.write("color,frequency")
        for item in data.items():
            csv_file.write(f"{item[0]},{item[1]}\n")


def plot_distribution(data_dict, title="Colors Distribution", xlabel="colors", ylabel="frequency", logscaled=True, filename="plot.pdf"):
    # Extracting keys and values
    x = list(data_dict.keys())
    y = list(data_dict.values())

    # Create the plot
    plt.figure(figsize=(10, 6))  # Set the size of the figure
    # plt.bar(x, y, color='skyblue', edgecolor='black')
    # plt.plot(x, y, marker='o', linestyle='-', color='skyblue', markerfacecolor='black', markeredgewidth=1)
    # plt.plot(x, y, marker=None, linestyle='-', color='skyblue', markerfacecolor='black', markeredgewidth=1)
    plt.scatter(x, y, s=1, color='blue', alpha=0.5)  # Use scatter for better performance, s=1 sets point size

    # Set x and y axes to logarithmic scale
    if logscaled:
        # plt.xscale('log')
        plt.yscale('log')

    # Add titles and labels
    plt.title(title, fontsize=14)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)

    # Rotate x-axis labels if necessary
    plt.xticks(rotation=45, ha="right", fontsize=10)

    # Add gridlines for better readability
    plt.grid(True, axis='y', linestyle='--', alpha=0.6)

    # Display the plot
    plt.tight_layout()
    plt.savefig(f"{filename}.pdf", format="pdf")


if __name__ == '__main__':
    sequence_filename = "test.fa"
    colors_filename = "colors.rle.colors"
    # colors_filename = "colored-graph-41.fasta.ustar.rle.colors"
    k = 3

    if len(sys.argv) == 4:
        sequence_filename = sys.argv[1]
        colors_filename = sys.argv[2]
        k = int(sys.argv[3])
    else:
        print("Usage: python colors-distributions.py <sequences.fa> <colors.rle> <k>")
        print("Going with defaults...")

    # collecting colors vectors
    colors = []
    sequence_lengths = []
    num_runs = 0
    try:
        print(f"Reading color file \"{colors_filename}\"...")
        with open(colors_filename, 'r') as colors_rle:
            for line in colors_rle:
                num_runs += 1
                rle = line.split(sep=":")
                try:
                    value = int(rle[0])
                    if len(rle) == 2:
                        count = int(rle[1])
                    else:
                        count = 1
                    for _ in range(count):
                        colors.append(value)
                except ValueError:
                    print(f"Can't decode: {line}")
                    exit(1)

        print(f"Reading sequence file \"{sequence_filename}\"...", flush=True)
        with open(sequence_filename, 'r') as sequences:
            for line in sequences:
                if line[0] == '>':  # found a header
                    continue  # need to get sequence information
                else:  # found the sequence
                    sequence = line.strip()
                    sequence_length = len(sequence)
                    kmer_num = sequence_length - k + 1
                    sequence_lengths.append(kmer_num)

    except FileNotFoundError:
        print("Can't open this file!")
        exit(1)

    print("Counting colors...", flush=True)
    colors_dist = {}
    for color in colors:
        colors_dist[color] = colors_dist.get(color, 0) + 1

    plot_distribution(colors_dist, title=f"Colors distribution, k={k}", filename=f"colors-distribution-k{k}")
    save_to_csv(colors_dist, f"colors-distribution-k{k}.csv")

    head_tail_dist = {}
    curr_kmer = 0
    singletons = 0
    for length in sequence_lengths:
        if length == 1:
            singletons += 1
            color = colors[curr_kmer]
            head_tail_dist[color] = head_tail_dist.get(color, 0) + 1
        else:
            head_color = colors[curr_kmer]
            head_tail_dist[head_color] = head_tail_dist.get(head_color, 0) + 1
            tail_color = colors[curr_kmer + length - 1]
            head_tail_dist[tail_color] = head_tail_dist.get(tail_color, 0) + 1
        curr_kmer = curr_kmer + length

    plot_distribution(head_tail_dist, title=f"Heads tails distribution, k={k}", filename=f"heads-tails-distribution-k{k}")
    save_to_csv(head_tail_dist, f"head-tails-distribution-k{k}.csv")

    print(f"Singletons: {singletons}")
    print(f"Total number of runs: {num_runs}")
    print(f"Total number of kmers: {len(colors)}", flush=True)
