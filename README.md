# Compact-and-Efficient-Representations-of-Colored-k-mer-Sets
Color Table Compression is a method inspired by md-Fulgor that compresses the GGCAT color table using similarity-based reordering. It computes MinHash signatures and applies K-Means clustering to group similar colors. Data is encoded with 64-bit bitmasks in a CSR matrix, reducing redundancy and enabling fast, lossless access.
