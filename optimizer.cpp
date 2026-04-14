#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cstring>
#include <cstdint> 
#include <omp.h>
#include <array>
#include <iomanip>
#include <sys/stat.h>

#include "dkm.hpp" 

// --- CONFIGURATION ---
const int SKETCH_SIZE = 64; 
const int BLOCK_WIDTH = 64; 

// CSR Structures
struct BinHeaderCSR {
    char magic[8];        
    uint32_t num_gids;    
    uint32_t num_blocks; 
    uint32_t block_size; 
    uint64_t num_nnz;     
};

// Utils
bool check_file(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

void draw_progress_bar(size_t current, size_t total, const std::string& label) {
    static int last_percent = -1;
    if (total == 0) return;
    int percent = (int)((current * 100.0) / total);
    if (percent != last_percent) {
        last_percent = percent;
        std::cout << "\r" << label << " " << percent << "% " << std::flush;
    }
}

uint64_t read_varint_stream(std::ifstream& is) {
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

void write_varint_stream(std::ofstream& os, uint64_t value) {
    while (value >= 0x80) {
        os.put((char)((value & 0x7F) | 0x80));
        value >>= 7;
    }
    os.put((char)value);
}

// Compute Reordering
std::vector<uint32_t> compute_reordering(const std::string& sigPath, const std::string& mapOutputPath, uint32_t& num_blocks_out) {
    std::cout << "[SCPO] Loading signatures from " << sigPath << "..." << std::endl;
    
    if (!check_file(sigPath)) {
        std::cerr << "[FATAL] Missing signatures file: " << sigPath << std::endl;
        exit(1);
    }

    std::ifstream sigFile(sigPath, std::ios::binary);
    uint32_t num_genomes, sig_len;
    sigFile.read((char*)&num_genomes, 4);
    sigFile.read((char*)&sig_len, 4);
    
    std::cout << "[SCPO] Genomes found: " << num_genomes << std::endl;

    if (num_genomes == 0) {
        std::cerr << "[FATAL] Number of genomes is 0. signatures.bin is empty or corrupted." << std::endl;
        exit(1);
    }

    std::vector<std::array<float, SKETCH_SIZE>> data;
    data.reserve(num_genomes);
    
    for (uint32_t i = 0; i < num_genomes; ++i) {
        std::vector<uint32_t> raw(SKETCH_SIZE);
        sigFile.read((char*)raw.data(), SKETCH_SIZE * 4);
        std::array<float, SKETCH_SIZE> p;
        for(int k=0; k<SKETCH_SIZE; ++k) p[k]=(float)raw[k];
        data.push_back(p);
    }
    sigFile.close();

    uint32_t K = num_genomes / BLOCK_WIDTH;
    if (K < 1) K = 1;
    if (K > 60000) K = 60000;

    std::cout << "[SCPO] K-Means Clustering (K=" << K << ")..." << std::endl;
    dkm::KMeans<float> kmeans(K);
    auto result = kmeans.cluster(data);
    const auto& labels = std::get<1>(result);

    struct GenomeInfo { uint32_t original_id; uint32_t cluster_id; };
    std::vector<GenomeInfo> genomes(num_genomes);
    for(uint32_t i=0; i<num_genomes; ++i) genomes[i] = {i, labels[i]};

    std::sort(genomes.begin(), genomes.end(), [](const GenomeInfo& a, const GenomeInfo& b) {
        if (a.cluster_id != b.cluster_id) return a.cluster_id < b.cluster_id;
        return a.original_id < b.original_id;
    });

    std::vector<uint32_t> permutation(num_genomes);
    std::vector<uint32_t> reverse_map(num_genomes);

    for(uint32_t new_id=0; new_id<num_genomes; ++new_id) {
        permutation[genomes[new_id].original_id] = new_id;
        reverse_map[new_id] = genomes[new_id].original_id;
    }

    std::ofstream mapFile(mapOutputPath, std::ios::binary);
    mapFile.write((char*)&num_genomes, 4);
    mapFile.write((char*)reverse_map.data(), num_genomes * sizeof(uint32_t));
    mapFile.close();
    
    num_blocks_out = (num_genomes + BLOCK_WIDTH - 1) / BLOCK_WIDTH;
    return permutation;
}

struct PartitionBlock {
    std::map<uint64_t, uint32_t> mask_to_pid;
    std::vector<uint64_t> pid_to_mask;
    PartitionBlock() {
        mask_to_pid[0] = 0;
        pid_to_mask.push_back(0);
    }
    uint32_t get_or_create_pid(uint64_t mask) {
        auto it = mask_to_pid.find(mask);
        if (it != mask_to_pid.end()) return it->second;
        uint32_t pid = (uint32_t)pid_to_mask.size();
        mask_to_pid[mask] = pid;
        pid_to_mask.push_back(mask);
        return pid;
    }
};

struct SparseEntry { uint32_t block_idx; uint32_t pid; };

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_compact.bin> <output_base>" << std::endl;
        return 1;
    }
    std::string inputPath = argv[1];
    std::string outputBase = argv[2];

    if (!check_file(inputPath)) {
        std::cerr << "[FATAL] Missing input file: " << inputPath << std::endl;
        return 1;
    }

    // 1. REORDERING
    uint32_t num_blocks = 0;
    std::string mapPath = outputBase + ".map.bin";
    std::vector<uint32_t> permutation = compute_reordering("signatures.bin", mapPath, num_blocks);
    
    std::vector<PartitionBlock> partitions(num_blocks);

    // 2. INPUT READING
    std::ifstream infile(inputPath, std::ios::binary);
    infile.seekg(0, std::ios::end); size_t total_bytes = infile.tellg(); infile.seekg(0, std::ios::beg);

    std::vector<std::vector<SparseEntry>> global_sparse_table;
    uint32_t max_gid = 0;
    int update_ctr = 0;

    std::cout << "[INFO] Bitmask Compression..." << std::endl;

    while (infile.peek() != EOF) {
        if (++update_ctr % 5000 == 0) draw_progress_bar(infile.tellg(), total_bytes, "Processing");

        uint32_t gid = (uint32_t)read_varint_stream(infile);
        uint32_t count = (uint32_t)read_varint_stream(infile);

        std::map<uint32_t, uint64_t> current_row_masks;

        for(uint32_t k=0; k<count; ++k) {
            uint32_t old_id = (uint32_t)read_varint_stream(infile);
            if (old_id < permutation.size()) {
                uint32_t new_id = permutation[old_id]; 
                
                uint32_t block_idx = new_id / BLOCK_WIDTH;
                uint32_t bit_idx   = new_id % BLOCK_WIDTH;

                current_row_masks[block_idx] |= (1ULL << bit_idx);
            }
        }
        
        if (gid > max_gid) max_gid = gid;
        if (global_sparse_table.size() <= gid) global_sparse_table.resize(gid + 1);

        std::vector<SparseEntry> row_entries;
        for (auto const& [b_idx, mask] : current_row_masks) {
            uint32_t pid = partitions[b_idx].get_or_create_pid(mask);
            row_entries.push_back({b_idx, pid});
        }
        global_sparse_table[gid] = row_entries;
    }
    std::cout << std::endl;

    // 3. FILE WRITING
    std::cout << "[INFO] Writing CSR..." << std::endl;
    std::ofstream outBin(outputBase + ".dat", std::ios::binary);
    uint64_t total_nnz = 0;
    for (const auto& r : global_sparse_table) total_nnz += r.size();

    BinHeaderCSR h; 
    memset(&h, 0, sizeof(h)); 
    strncpy(h.magic, "MDFULVR1", 8);
    h.num_gids = global_sparse_table.size(); 
    h.num_blocks = num_blocks; 
    h.block_size = BLOCK_WIDTH;
    h.num_nnz = total_nnz; 
    outBin.write((char*)&h, sizeof(h));

    uint64_t curr_off = 0;
    outBin.write((char*)&curr_off, 8);
    for (const auto& r : global_sparse_table) {
        curr_off += r.size();
        outBin.write((char*)&curr_off, 8);
    }

    for (const auto& r : global_sparse_table) 
        for (const auto& e : r) write_varint_stream(outBin, (uint64_t)e.block_idx); 

    for (const auto& r : global_sparse_table) 
        for (const auto& e : r) write_varint_stream(outBin, (uint64_t)e.pid);
    
    outBin.close();

    // 4. DICTIONARY WRITING
    std::cout << "[INFO] Writing Dictionary..." << std::endl;
    std::ofstream outDict(outputBase + ".dict.bin", std::ios::binary);
    outDict.write((char*)&num_blocks, 4);

    for (uint32_t i=0; i<partitions.size(); ++i) {
        const auto& part = partitions[i];
        uint32_t np = part.pid_to_mask.size();
        outDict.write((char*)&i, 4); 
        outDict.write((char*)&np, 4);
        for (uint64_t mask : part.pid_to_mask) outDict.write((char*)&mask, 8);
    }
    outDict.close();

    std::cout << "\n[COMPLETED] SUCCESS" << std::endl;
    return 0;
}