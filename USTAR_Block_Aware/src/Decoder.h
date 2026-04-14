//
// Created by enrico on 29/12/22.
// Modified by Marco Squarcina on 14/04/26.
//

#ifndef USTAR_DECODER_H
#define USTAR_DECODER_H

#include <string>
#include <vector>
#include <cstdint>

#include "consts.h"

using namespace std;
class Decoder{
    string fasta_file_name;
    string counts_file_name;
    int kmer_size;
    bool debug;

    vector<uint32_t> counts;

public:
    /**
     * Decode and extract kmers from simplitigs
     * @param fasta_file_name USTAR or UST fasta file
     * @param counts_file_name USTAR or UST colors file
     * @param kmer_size kmer size
     * @param debug debug flag
     */
    Decoder(const string &fasta_file_name, const string &counts_file_name, int kmer_size, bool debug);

    /**
     * Extract kmers and colors and write them to a txt file
     * @param output_file_name where to put kmers and colors
     */
    void extract_kmers_and_counts(const string &output_file_name);

    /**
     * Reads and decode colors
     * @param encoding encoding used by USTAR
     */
    void decode(encoding_t encoding);
};
#endif //USTAR_DECODER_H
