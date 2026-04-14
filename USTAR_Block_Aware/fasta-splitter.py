import sys

fasta_filename = sys.argv[1]
sequence_filename = "sequences.fa"
headers_filename = "headers.txt"

with open(fasta_filename, 'r') as fasta:
    with open(sequence_filename, 'w') as sequence:
        with open(headers_filename, 'w') as headers:
            for line in fasta:
                if line[0] == '>':
                    headers.write(line)
                else:
                    sequence.write(f">\n{line}")
