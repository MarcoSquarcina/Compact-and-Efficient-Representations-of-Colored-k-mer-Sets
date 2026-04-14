//
// Created by Marco Squarcina on 14/04/26.
//

#ifndef MDMAP_HPP
#define MDMAP_HPP

#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>

class mdMap {
public:
    // Structure for the .dat file header
    struct Header {
        char magic[8];
        uint32_t num_gids;
        uint32_t num_blocks;
        uint32_t block_size;
        uint64_t num_nnz;
    };

private:
    // CSR (Compressed Sparse Row) data loaded into RAM
    std::vector<uint64_t> offsets;      // Cumulative indices (start/end for each GID)
    std::vector<uint32_t> flat_blocks;  // Flat vector of Block IDs
    std::vector<uint32_t> flat_pids;    // Flat vector of Pattern IDs

    // Dictionary of Bitmasks (BlockID -> [PatternID -> Bitmask])
    std::vector<std::vector<uint64_t>> dictionary;
    
    bool loaded = false;

    // Helper to read VarInt (LEB128) from the stream
    static uint64_t read_varint(std::ifstream& is) {
        uint64_t value = 0;
        uint64_t shift = 0;
        char byte;
        while (is.get(byte)) {
            value |= (uint64_t)(byte & 0x7F) << shift;
            if ((byte & 0x80) == 0) return value;
            shift += 7;
        }
        return 0;
    }

public:
    mdMap() = default;
    
    explicit mdMap(const std::string& base_filename) {
        load(base_filename);
    }

    // Load the .dat and .dict.bin files (new VarInt format)
    bool load(const std::string& base_filename) {
        loaded = false;
        std::string dat_file = base_filename + ".dat";
        std::string dict_file = base_filename + ".dict.bin";

        std::cout << "[mdMap] Loading new VarInt format from: " << base_filename << "..." << std::endl;

        // 1. LOADING .DAT (CSR Structure)
        std::ifstream inDat(dat_file, std::ios::binary);
        if (!inDat) {
            std::cerr << "Error: Cannot open .dat file: " << dat_file << std::endl;
            return false;
        }

        Header h;
        inDat.read((char*)&h, sizeof(Header));

        if (std::strncmp(h.magic, "MDFULVR1", 8) != 0) {
            if (std::strncmp(h.magic, "MDFULCSR", 8) == 0) {
                std::cerr << "Error: Old format detected (MDFULCSR). Please update dataset or use legacy code." << std::endl;
                return false;
            }
            std::cerr << "Error: Invalid Magic Number in .dat file. Expected MDFULVR1." << std::endl;
            return false;
        }

        std::cout << "[mdMap] Header Info: " << h.num_gids << " GIDs, " 
                  << h.num_blocks << " Blocks, " << h.num_nnz << " NNZ." << std::endl;

        // 1.1 Read Offsets (These are fixed at 64 bits, not compressed)
        offsets.resize(h.num_gids + 1);
        inDat.read((char*)offsets.data(), (h.num_gids + 1) * sizeof(uint64_t));

        // 1.2 Read and Decompress Block IDs (VarInt Stream)
        std::cout << "[mdMap] Decompressing Block IDs..." << std::endl;
        flat_blocks.clear();
        flat_blocks.reserve(h.num_nnz);
        for(size_t i = 0; i < h.num_nnz; ++i) {
            flat_blocks.push_back((uint32_t)read_varint(inDat));
        }

        // 1.3 Read and Decompress Pattern IDs (VarInt Stream)
        std::cout << "[mdMap] Decompressing Pattern IDs..." << std::endl;
        flat_pids.clear();
        flat_pids.reserve(h.num_nnz);
        for(size_t i = 0; i < h.num_nnz; ++i) {
            flat_pids.push_back((uint32_t)read_varint(inDat));
        }
        
        inDat.close();

        // 2. LOADING .DICT.BIN (Dictionary of Bitmasks)
        std::ifstream inDict(dict_file, std::ios::binary);
        if (inDict) {
            std::cout << "[mdMap] Loading Bitmask Dictionary..." << std::endl;
            uint32_t num_blocks_dict;
            inDict.read((char*)&num_blocks_dict, sizeof(uint32_t));

            if (num_blocks_dict != h.num_blocks) {
                 std::cerr << "[mdMap] Warning: Mismatch in block count between .dat and .dict" << std::endl;
            }

            if (dictionary.size() < num_blocks_dict) {
                dictionary.resize(num_blocks_dict);
            }

            for(uint32_t b = 0; b < num_blocks_dict; ++b) {
                if (inDict.peek() == EOF) break;

                uint32_t block_id;
                uint32_t num_patterns;
                inDict.read((char*)&block_id, sizeof(uint32_t));
                inDict.read((char*)&num_patterns, sizeof(uint32_t));

                if (block_id < dictionary.size()) {
                    dictionary[block_id].reserve(num_patterns);
                    for(uint32_t p = 0; p < num_patterns; ++p) {
                        uint64_t mask;
                        inDict.read((char*)&mask, sizeof(uint64_t));
                        dictionary[block_id].push_back(mask);
                    }
                } else {
                    inDict.seekg(num_patterns * 8, std::ios::cur);
                }
            }
            inDict.close();
        } else {
            std::cerr << "[mdMap] Warning: .dict.bin not found. Running in ID-only mode." << std::endl;
        }

        loaded = true;
        std::cout << "[mdMap] Load complete." << std::endl;
        return true;
    }

    uint64_t get_mapped_value(uint32_t gid) const {
        if (gid >= offsets.size() - 1) return 0;

        uint64_t start_idx = offsets[gid];
        uint64_t end_idx = offsets[gid + 1];

        if (start_idx >= end_idx) return 0;
        
        if (start_idx >= flat_blocks.size()) return 0;

        uint32_t block = flat_blocks[start_idx];
        uint32_t pid = flat_pids[start_idx];

        return (static_cast<uint64_t>(block) << 32) | static_cast<uint64_t>(pid);
    }

    bool is_loaded() const { return loaded; }
    
    size_t get_num_gids() const { return offsets.size() > 0 ? offsets.size() - 1 : 0; }
};

#endif // MDMAP_HPP