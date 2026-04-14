# Optimizer & USTAR Block-Aware

This repository contains the implementation of the tools described in the paper "Compact and Efficient Representations of Colored k-mer Sets" by Marco Squarcina. The project was developed at the Department of Information Engineering of the University of Padua.

## Project Description

The analysis of large-scale pangenomic and metagenomic datasets requires compact data structures, such as Colored de Bruijn Graphs (c-dBG). State-of-the-art tools like GGCAT are extremely fast at building these graphs, but they suffer from the "color table explosion" problem, producing extremely large color output files that exceed the size of the genomic files themselves.

This project introduces two main components to solve this bottleneck:
* **Optimizer**: A method inspired by the md-Fulgor tool that drastically compresses the color table produced by GGCAT.
* **USTAR Block-Aware**: A modified version of the USTAR tool that accepts the optimized files to assemble simplitigs, reducing storage space and maintaining color coherence.

## Main Features

### 1. Optimizer
The optimization process eliminates data redundancy through several phases:
* **Sketching**: Computes compact MinHash signatures to capture similarities between sets containing similar colors.
* **Clustering-based Reordering**: Applies the K-Means algorithm to group and reorder data, maximizing data locality.
* **Bitmask Compression**: Replaces explicit representations with 64-bit bitmasks inside a Compressed Sparse Row (CSR) matrix.

### 2. Generated Outputs
The Optimizer produces three optimized binary files designed for fast random access:
* **`.dat` (Sparse Matrix)**: Maps Global IDs (GIDs) to Pattern IDs in CSR format.
* **`.dict.bin` (BitMask Dictionary)**: Contains the bitmask dictionaries divided by blocks. Each 64-bit integer represents the colors present in the set.
* **`.map.bin` (Permutation Map)**: Saves the mapping required to reconstruct the original order of the color sets after the clustering phase.

### 3. USTAR Block-Aware
The modified version of USTAR accesses the optimized `.dat` file and applies a three-level decision process during path extension for simplitig creation:
* **Priority 1 (Exact Match)**: Favors connections between nodes with the exact same color set ID (GID), fully preserving coherence.
* **Priority 2 (Block Match)**: Checks if the IDs belong to the same color block, improving subsequent compression.
* **Priority 3 (Less Connected)**: Selects the least connected neighbor among the available candidates to minimize fragmentation.

## Performance and Results

Tests conducted on large-scale real datasets demonstrate the effectiveness of the pipeline:
* On the *Salmonella Enterica* pangenomic dataset (up to 150,000 genomes), Optimizer achieves a 97% compression rate compared to the standard GGCAT color table.
* On the *Gut Bacteria* metagenomic dataset (up to 30,000 genomes), characterized by high fragmentation, it achieves a 33% compression compared to GGCAT.
* Using the complete pipeline (Optimizer + USTAR Block-Aware) on 150,000 Salmonella genomes allows storing the entire set of information (topology and colors) in just 2,671 MB, achieving a 96% disk space reduction compared to GGCAT and its variants.
* The method outperforms the state-of-the-art md-Fulgor approach for large-scale datasets, offering an additional space reduction of up to 43.9% on the Gut Bacteria dataset.

## Dependencies
* **DKM**: C++ library used to perform K-Means clustering.
* **LZ4**: Used for decompression when reading the files originally written by GGCAT.