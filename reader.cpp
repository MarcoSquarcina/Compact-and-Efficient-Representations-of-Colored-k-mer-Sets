#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <thread>
#include <chrono>
#include <sys/stat.h>

// --- CONFIGURATION ---
const int SKETCH_SIZE = 64;
const int ESTIMATED_GENOMES = 150000;
const uint64_t MAX_ALLOWED_GENOME_ID = 20000000;

// Hash function
inline uint32_t hash_unitig(uint32_t key, uint32_t seed) {
    key = key ^ seed;
    key = (key + 0x7ed55d16) + (key << 12);
    key = (key ^ 0xc761c23c) ^ (key >> 19);
    key = (key + 0x165667b1) + (key << 5);
    key = (key + 0xd3a2646c) ^ (key << 9);
    key = (key + 0xfd7046c5) + (key << 3);
    key = (key ^ 0xb55a4f09) ^ (key >> 16);
    return key;
}

uint64_t read_varint(const char*& ptr, const char* end) {
    uint64_t value = 0;
    uint64_t shift = 0;
    while (ptr < end) {
        uint8_t byte = *ptr++;
        value |= (uint64_t)(byte & 0x7F) << shift;
        if ((byte & 0x80) == 0) return value;
        shift += 7;
        if (shift > 64) return UINT64_MAX;
    }
    return 0;
}

void write_varint(std::ofstream& outFile, uint64_t value) {
    while (value >= 0x80) {
        outFile.put((char)((value & 0x7F) | 0x80));
        value >>= 7;
    }
    outFile.put((char)value);
}

// Structure for managing physical sorting
struct ChunkInfo {
    uint32_t start_unitig;
    uint64_t offset;
    size_t original_index; // Per debug
};

int main(int argc, char* argv[]) {
    // Bufferless output for real-time debug
    std::setvbuf(stdout, NULL, _IONBF, 0);
    std::setvbuf(stderr, NULL, _IONBF, 0);

    if (argc < 5) {
        std::cerr << "Uso: " << argv[0] << " <input> <output> <tmp_lz4> <tmp_bin>" << std::endl;
        return 1;
    }

    std::string filePath = argv[1];
    std::string outPath = argv[2];
    std::string tmpLz4Path = argv[3];
    std::string tmpBinPath = argv[4];

    std::cout << "--- READER PHYSICAL REORDER MODE ---" << std::endl;
    std::ifstream file(filePath, std::ios::binary);
    if (!file) { std::cerr << "[FATAL] Input non trovato." << std::endl; return 1; }

    std::ofstream outFile(outPath, std::ios::binary);
    std::vector<std::vector<uint32_t>> sketches;
    sketches.reserve(ESTIMATED_GENOMES);
    std::vector<uint32_t> current_unitig_hashes(SKETCH_SIZE);

    // 1. Read Header
    uint64_t index_offset;
    file.seekg(24, std::ios::beg);
    file.read(reinterpret_cast<char*>(&index_offset), 8);

    std::cout << "Index Offset (End of Data): " << index_offset << std::endl;

    // 2. Read Index
    file.seekg(index_offset, std::ios::beg);
    uint64_t entries_count;
    file.read(reinterpret_cast<char*>(&entries_count), 8);

    std::vector<ChunkInfo> chunks;
    chunks.reserve(entries_count);

    for(uint64_t i=0; i<entries_count; ++i) {
        ChunkInfo c;
        file.read(reinterpret_cast<char*>(&c.start_unitig), 4);
        file.read(reinterpret_cast<char*>(&c.offset), 8);
        c.original_index = i;
        chunks.push_back(c);
    }

    std::cout << "Indice letto. Riordinamento chunk per offset fisico..." << std::endl;

    // --- Sort by ascending offset ---
    std::sort(chunks.begin(), chunks.end(), [](const ChunkInfo& a, const ChunkInfo& b) {
        return a.offset < b.offset;
    });
    // ----------------------------------------------

    // 3. Sequential Elaboration
    for (size_t i = 0; i < chunks.size(); ++i) {
        uint64_t currentOffset = chunks[i].offset;

        // Compute the size using the NEXT PHYSICAL chunk (or the end of the data)
        uint64_t nextOffset = (i == chunks.size() - 1) ? index_offset : chunks[i+1].offset;

        if (nextOffset < currentOffset) {
            std::cerr << "[FATAL] Offset non monotoni anche dopo il sort! File corrotto gravemente." << std::endl;
            return 1;
        }

        uint64_t chunkSize = nextOffset - currentOffset;

        if (i < 5) {
             std::cout << "Proc. Chunk #" << chunks[i].original_index
                       << " (ID " << chunks[i].start_unitig << ")"
                       << " Offset: " << currentOffset
                       << " Size: " << chunkSize << std::endl;
        }

        if (chunkSize == 0) continue;

        // --- Reading and Decompression ---
        file.clear();
        file.seekg(currentOffset);
        std::vector<char> compData(chunkSize);
        if (!file.read(compData.data(), chunkSize)) {
             std::cerr << "[FATAL] Errore lettura chunk fisico " << i << std::endl;
             return 1;
        }

        {
            std::ofstream tmp(tmpLz4Path, std::ios::binary);
            tmp.write(compData.data(), chunkSize);
            tmp.flush();
        }

        std::string cmd = "lz4 -d -q -f " + tmpLz4Path + " " + tmpBinPath;
        if (system(cmd.c_str()) != 0) {
            std::cerr << "[FATAL] LZ4 fallito." << std::endl;
            return 1;
        }

        std::ifstream rawFile(tmpBinPath, std::ios::binary | std::ios::ate);
        size_t rawSize = rawFile.tellg();
        rawFile.seekg(0);
        std::vector<char> rawData(rawSize);
        rawFile.read(rawData.data(), rawSize);
        rawFile.close();

        // --- Parsing ---
        const char* ptr = rawData.data();
        const char* end = ptr + rawSize;
        // IMPORTANT: Use the ID stored in the chunk
        uint32_t currentUnitigId = chunks[i].start_unitig;

        while (ptr < end) {
            uint64_t count = read_varint(ptr, end);

            for(int k=0; k < SKETCH_SIZE; ++k) {
                current_unitig_hashes[k] = hash_unitig(currentUnitigId, k * 1000 + 123);
            }
            write_varint(outFile, currentUnitigId);
            write_varint(outFile, count);

            for(uint64_t k=0; k<count; ++k) {
                uint64_t genomeId = read_varint(ptr, end);
                write_varint(outFile, genomeId);

                // Safe Resize
                if (genomeId >= sketches.size()) {
                    if (genomeId > sketches.capacity() + 200000) sketches.reserve(genomeId + 100000);
                    sketches.resize(genomeId + 1, std::vector<uint32_t>(SKETCH_SIZE, UINT32_MAX));
                }
                for(int s=0; s < SKETCH_SIZE; ++s) {
                    if (current_unitig_hashes[s] < sketches[genomeId][s]) {
                        sketches[genomeId][s] = current_unitig_hashes[s];
                    }
                }
            }
            currentUnitigId++;
        }

        if (i % 100 == 0) std::cout << "\rProgress: " << i << "/" << chunks.size() << " chunks physical sorted" << std::flush;
    }

    std::cout << "\nScrittura firme finali..." << std::endl;
    // Remove temp
    std::remove(tmpLz4Path.c_str());
    std::remove(tmpBinPath.c_str());
    outFile.close();

    // Signature writing
    std::ofstream sigFile("signatures.bin", std::ios::binary);
    uint32_t num_genomes = sketches.size();
    uint32_t signature_len = SKETCH_SIZE;
    sigFile.write(reinterpret_cast<char*>(&num_genomes), 4);
    sigFile.write(reinterpret_cast<char*>(&signature_len), 4);
    for (const auto& sig : sketches) {
        sigFile.write(reinterpret_cast<const char*>(sig.data()), signature_len * sizeof(uint32_t));
    }
    sigFile.close();

    std::cout << "SUCCESS." << std::endl;
    return 0;
}
