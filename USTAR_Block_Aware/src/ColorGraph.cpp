//
// Created by enrico on 29/12/22.
// Modified by Marco Squarcina on 14/04/26.
//

#include "ColorGraph.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <unordered_map>
#include <limits>
#include <iomanip>
#include <tuple>
#include <map>

using namespace std;
using namespace std::chrono;

// --- Helper VarInt utilities ---
static inline void write_varint32(ofstream &out, uint32_t value) {
    while (value >= 0x80) {
        out.put(static_cast<char>((value & 0x7F) | 0x80));
        value >>= 7;
    }
    out.put(static_cast<char>(value & 0x7F));
}

static inline void write_varint64(ofstream &out, uint64_t value) {
    while (value >= 0x80) {
        out.put(static_cast<char>((value & 0x7F) | 0x80));
        value >>= 7;
    }
    out.put(static_cast<char>(value & 0x7F));
}

static inline uint64_t zigzag_encode(int64_t value) {
    return (static_cast<uint64_t>(value) << 1) ^ static_cast<uint64_t>(value >> 63);
}

inline char complement(char c){
    switch(c){
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'a': return 't';
        case 't': return 'a';
        case 'c': return 'g';
        case 'g': return 'c';
        default: return 'N';
    }
}

std::string reverse_complement(std::string sequence){
    for(auto &c: sequence) c = complement(c);
    std::reverse(sequence.begin(), sequence.end());
    return sequence;
}

// --- HYBRID RLE FUNCTION ---
// Uses "mapped_colors" to decide whether to extend the run (logical comparison).
// Uses "original_colors" to store the representative value (GID) in the file.
void encode_RLE_hybrid(const std::vector<uint64_t>& mapped_colors, 
                       const std::vector<color_id_t>& original_colors, 
                       std::vector<color_id_t> &out_values, 
                       std::vector<size_t> &out_counts) {
    if (mapped_colors.empty()) return;

    // First element initialization
    uint64_t prev_map = mapped_colors[0];
    color_id_t prev_orig = original_colors[0];
    size_t count = 1;

    // If we are continuing a previous path, try to merge (join)
    // Note: here we check whether the last inserted GID maps to the same block as the new one.
    // For simplicity and safety in this implementation, we always reset on node/path change
    // if the output vector is not empty, but we only check within the current vector.
    

    for(size_t i = 1; i < mapped_colors.size(); i++){
        if(prev_map == mapped_colors[i]) {
            // Same pattern in the DB (Block + PID are equal) , extend the run
            count++;
        } else {
            // Different pattern , save the previous GID and reset
            out_values.push_back(prev_orig);
            out_counts.push_back(count);
            
            prev_map = mapped_colors[i];
            prev_orig = original_colors[i];
            count = 1;
        }
    }
    // Final push for the last run
    out_values.push_back(prev_orig);
    out_counts.push_back(count);
}

// --- Constructors ---
ColorGraph::ColorGraph(std::string sequences_file_name, std::string colors_file_name, int kmer_length, mdMap* map, bool debug_flag) {
    this->sequences_file_name = std::move(sequences_file_name);
    this->colors_file_name = std::move(colors_file_name);
    this->kmer_length = kmer_length;
    this->map_ptr = map;
    this->debug = debug_flag;

    if(kmer_length < 3){
        cerr << "Error ColorGraph(): kmer_length must be at least 3" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "* Building the graph..." << endl;
    build_graph();
    print_stats();
    cout << "* Computing a path cover..." << endl;
    compute_path_cover();
    cout << "* Path cover ready!" << endl;
}

ColorGraph::ColorGraph(const std::vector<std::string>& sequences, const std::vector<std::vector<color_id_t>>& colors,
                       int kmer_length, mdMap* map, bool debug_flag) {
    this->sequences_file_name = "tmp";
    this->colors_file_name = "tmp";
    this->map_ptr = map; 
    this->debug = debug_flag;

    if(kmer_length == 0)
        kmer_length = sequences[0].length() - colors[0].size() + 1;
    this->kmer_length = kmer_length;

    if(kmer_length < 3){
        cerr << "Error ColorGraph(): kmer_length must be at least 3" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "* Building the graph..." << endl;
    build_graph(sequences, colors);
    print_stats();
    cout << "* Computing a path cover..." << endl;
    compute_path_cover();
    cout << "* Path cover ready!" << endl;
}

void ColorGraph::compute_path_cover() {
    auto start = chrono::steady_clock::now();

    // Ordering nodes by degree (base heuristic)
    vector<node_id_t> degrees;
    degrees.reserve(nodes.size());
    for(node_id_t i = 0; i < nodes.size(); i++){
        auto d = nodes_head[nodes[i].colors.front()].size() + nodes_tail[nodes[i].colors.front()].size() +
                 nodes_head[nodes[i].colors.back()].size() + nodes_tail[nodes[i].colors.back()].size();
        degrees.push_back(d);
    }
    auto conn_key = [&degrees](node_id_t a, node_id_t b) {return degrees[a] < degrees[b];};

    vector<node_id_t> order;
    order.reserve(nodes.size());
    for(node_id_t i = 0; i < nodes.size(); i++)
        order.push_back(i);
    sort(order.begin(), order.end(), conn_key);

    auto stop = chrono::steady_clock::now();
    cout << "compute_path_cover() [sort nodes]: " << duration_cast<seconds>(stop - start).count() << " seconds" << endl;

    // Ordering arcs
    for(auto colors_adj: nodes_head){
        auto &neighbours = colors_adj.second;
        neighbours.sort(conn_key);
    }
    for(auto colors_adj: nodes_tail){
        auto &neighbours = colors_adj.second;
        neighbours.sort(conn_key);
    }
    stop = chrono::steady_clock::now();
    cout << "compute_path_cover() [sort arcs]: " << duration_cast<seconds>(stop - start).count() << " seconds" << endl;

    start = chrono::steady_clock::now();
    // Greedy exploration
    for(auto& node_id: order){
        auto &seed = nodes[node_id];
        auto seed_id = node_id;

        if(seed.is_visited()) continue;

        seed.visit();
        Path path(seed_id);

        while(has_next(path)) path.extend(next(path));

        path.reverse();
        while(has_next(path)) path.extend(next(path));

        paths.push_back(path);
    }
    stop = chrono::steady_clock::now();
    cout << "compute_path_cover() [explore]: " << duration_cast<seconds>(stop - start).count() << " seconds" << endl;

    // Added sorting of simplitigs based on colors
    start = chrono::steady_clock::now();
    sort_paths(); 

    start = chrono::steady_clock::now();
    finalize_path_cover();
    stop = chrono::steady_clock::now();
    cout << "compute_path_cover() [finalize]: " << duration_cast<seconds>(stop - start).count() << " seconds" << endl; 
}

void ColorGraph::build_graph() {
    auto start = steady_clock::now();
    vector<color_id_t> colors = decode_RLE_colors();
    ifstream sequences_file(sequences_file_name);
    
    if(!sequences_file.is_open()){
        cerr << "Error build_graph(): cannot open file " << sequences_file_name << endl;
        exit(EXIT_FAILURE);
    }

    cout << "** Reading sequences " << sequences_file_name << endl;

    node_id_t node_id = 0;
    string line;
    while(getline(sequences_file, line)){
        if(line[0] == '>') continue;

        long length = static_cast<long>(line.size());
        if(length < kmer_length){
            cerr << "Error build_graph(): sequence too short!" << endl;
            exit(EXIT_FAILURE);
        }

        long n_kmer = length - kmer_length + 1;
        vector<color_id_t> sequence_colors = vector(colors.begin() + tot_kmers, colors.begin() + tot_kmers + n_kmer);

        nodes[node_id] = Node(line, sequence_colors);
        nodes_head[sequence_colors.front()].push_back(node_id);
        nodes_tail[sequence_colors.back()].push_back(node_id);

        tot_kmers += n_kmer;
        node_id++;
    }
    if(tot_kmers != static_cast<long>(colors.size())){
        cerr << "Error build_graph(): wrong number of colors" << endl;
        exit(EXIT_FAILURE);
    }

    auto stop = chrono::steady_clock::now();
    cout << "build_graph(): " << duration_cast<seconds>(stop - start).count() << " seconds" << endl;
}

void ColorGraph::build_graph(const std::vector<std::string> &sequences, const std::vector<std::vector<color_id_t>> &colors) {
    assert(sequences.size() == colors.size());
    auto start = chrono::steady_clock::now();

    node_id_t node_id = 0;
    for(size_t i = 0; i < sequences.size(); i++){
        const string &sequence = sequences[i];
        const vector<color_id_t> &sequence_colors = colors[i];
        long length = static_cast<long>(sequence.size());

        if(length < kmer_length){
            cerr << "Error build_graph(): sequence too short!" << endl;
            exit(EXIT_FAILURE);
        }

        long n_kmer = length - kmer_length + 1;

        nodes[node_id] = Node(sequence, sequence_colors);
        nodes_head[sequence_colors.front()].push_back(node_id);
        nodes_tail[sequence_colors.back()].push_back(node_id);

        tot_kmers += n_kmer;
        node_id++;
    }
    auto stop = chrono::steady_clock::now();
    cout << "build_graph(): " << duration_cast<seconds>(stop - start).count() << " seconds" << endl;
}

void ColorGraph::print_stats() const{
    cout << "\nGraph stats:\n";
    cout << "   number of kmers:        " << tot_kmers << "\n";
    cout << "   number of sequences:    " << nodes.size() << "\n";
    cout << endl;
}

std::vector<color_id_t> ColorGraph::decode_RLE_colors() {
    ifstream colors_file(colors_file_name);
    if(!colors_file.is_open()){
        cerr << "Error decode_RLE_colors(): cannot open file " << colors_file_name << endl;
        exit(EXIT_FAILURE);
    }

    std::vector<color_id_t> colors;
    string line;
    ulong value = 0;
    ulong count = 0;
    ulong pos = 0;

    cout << "** Reading colors " << colors_file_name << "..." << endl;

    while(colors_file >> line){
        pos = line.find(':');
        if(pos == string::npos) count = 1;
        else count = stoul(line.substr(pos + 1, string::npos));

        value = stoul(line.substr(0, pos));
        for(ulong i = 0; i < count; i++){
            colors.push_back(value);
        }
    }
    return colors;
}

bool ColorGraph::has_next(Path path) {
    color_id_t color;
    auto tail = path.get_tail_node_id();

    if(path.get_tail_orientation() == orientation_t::direct)
        color = nodes[tail].colors.back();
    else
        color = nodes[tail].colors.front();

    // Limit attempts to avoid infinite searches in disconnected components
    size_t max_attempts = nodes.size();
    for(size_t i = 0; i < max_attempts; i++) {
        if (!nodes_head[color].empty()) {
            next_node = nodes_head[color].front();
            if (nodes[next_node].is_visited()) {
                nodes_head[color].pop_front();
                continue;
            }
            next_orientation = orientation_t::direct;
            return true;
        }

        if (!nodes_tail[color].empty()) {
            next_node = nodes_tail[color].front();
            if (nodes[next_node].is_visited()) {
                nodes_tail[color].pop_front();
                continue;
            }
            next_orientation = orientation_t::reverse;
            return true;
        }
        return false;
    }
    return false;
}

oriented_node_t ColorGraph::next(Path path) {
    nodes[next_node].visit();
    return oriented_node_t{next_node, next_orientation};
}

void ColorGraph::finalize_path_cover(){
    // clear values and counts before filling them with the new RLE logic
    values.clear();
    counts.clear();

    for(auto &path: paths){
        for(size_t i = 0; i < path.length(); i++){
            node_id_t node_id = path.get_n_node_id(i);
            Node &node = nodes.at(node_id);
            if(path.get_n_orientation(i) == orientation_t::reverse)
                node.reverse();

            // --- Modified logic ---
            if(map_ptr != nullptr) {
                // 1. Prepare the vector of mapped values (for comparison)
                std::vector<uint64_t> mapped_colors;
                mapped_colors.reserve(node.colors.size());
                
                for(auto gid : node.colors) {
                    mapped_colors.push_back(map_ptr->get_mapped_value(gid));
                }
                
                // 2. Call the hybrid RLE (compare mapped, store original)
                encode_RLE_hybrid(mapped_colors, node.colors, values, counts);
            } 
            else {
                // RLE Standard (legacy, no mapping, just for testing/debugging)
                if(node.colors.empty()) continue;
                
                color_id_t prev = node.colors[0]; 
                size_t cnt = 1;
                
                // Try join with previous path/node end
                if(!counts.empty() && node.colors[0] == values.back()){
                    cnt += counts.back(); 
                    values.pop_back(); 
                    counts.pop_back();
                }
                
                for(size_t k=1; k<node.colors.size(); k++){
                    if(prev == node.colors[k]) cnt++;
                    else { 
                        values.push_back(prev); 
                        counts.push_back(cnt); 
                        prev = node.colors[k]; 
                        cnt=1; 
                    }
                }
                values.push_back(prev); 
                counts.push_back(cnt);
            }
        }
    }
    assert(values.size() == counts.size());
}

void ColorGraph::write_sequences(std::string sequences_filename) {
    cout << "** Writing sequences to " << sequences_filename << endl;
    ofstream sequences_file(sequences_filename);

    for(auto &path: paths){
        for(size_t i = 0; i < path.length(); i++){
            node_id_t node_id = path.get_n_node_id(i);
            Node &node = nodes.at(node_id);
            sequences_file << ">\n" << node.sequence << "\n";
        }
    }
}

void ColorGraph::write_colors(std::string colors_filename){
    cout << "** Writing colors to " << colors_filename << "\n";
    cout << "       number of runs: " << counts.size() << endl;
    
    ofstream colors_file(colors_filename, ios::binary); 
    
    if (map_ptr != nullptr) {
        cout << "[md-Fulgor] Writing COMPRESSED BINARY colors (block-aware PID delta)..." << endl;

        struct BlockStats {
            uint64_t runs = 0;
            uint64_t kmers = 0;
        };

        unordered_map<uint32_t, BlockStats> block_stats;
        if(debug)
            block_stats.reserve(1024);
        streamsize old_precision = 0;
        ios::fmtflags old_flags{};
        if(debug) {
            old_flags = cout.flags();
            old_precision = cout.precision();
        }

        uint64_t total_runs = 0;
        uint64_t total_kmers = 0;
        uint64_t sum_run = 0;
        uint32_t min_run = std::numeric_limits<uint32_t>::max();
        uint32_t max_run = 0;
        uint64_t block_changes = 0;
        uint32_t prev_block = std::numeric_limits<uint32_t>::max();
        uint32_t prev_pid = 0;
        bool have_prev = false;
        uint64_t pid_delta_sum = 0;
        uint64_t pid_delta_count = 0;
        const size_t sample_cap = 8;
        vector<tuple<size_t, uint32_t, uint32_t, uint32_t>> samples;
        if(debug)
            samples.reserve(sample_cap);

        const char magic[8] = {'C','L','R','M','A','P','2','\0'};
        colors_file.write(magic, sizeof(magic));
        write_varint64(colors_file, static_cast<uint64_t>(values.size()));

        for(size_t i = 0; i < values.size(); i++) {
            if(counts[i] > std::numeric_limits<uint32_t>::max()) {
                cerr << "write_colors(): run length exceeds 32-bit limit" << endl;
                exit(EXIT_FAILURE);
            }

            uint64_t composite = map_ptr->get_mapped_value(values[i]);
            uint32_t block = static_cast<uint32_t>(composite >> 32);
            uint32_t pid = static_cast<uint32_t>(composite & 0xFFFFFFFFULL);
            uint32_t run_len = static_cast<uint32_t>(counts[i]);

            bool block_changed = (prev_block == std::numeric_limits<uint32_t>::max()) || (block != prev_block);
            if(block_changed) {
                write_varint32(colors_file, 0);
                write_varint32(colors_file, block);
                write_varint32(colors_file, pid);
            } else {
                int64_t delta = static_cast<int64_t>(pid) - static_cast<int64_t>(prev_pid);
                write_varint64(colors_file, zigzag_encode(delta) + 1);
            }
            write_varint32(colors_file, run_len);

            if(debug) {
                total_runs++;
                total_kmers += run_len;
                sum_run += run_len;
                min_run = std::min(min_run, run_len);
                max_run = std::max(max_run, run_len);

                if(block_changed && have_prev)
                    block_changes++;

                if(have_prev && !block_changed) {
                    uint64_t delta_abs = (pid >= prev_pid) ? (pid - prev_pid) : (prev_pid - pid);
                    pid_delta_sum += delta_abs;
                    pid_delta_count++;
                }

                auto &entry = block_stats[block];
                entry.runs++;
                entry.kmers += run_len;

                if(samples.size() < sample_cap)
                    samples.emplace_back(i, block, pid, run_len);
            }

            prev_block = block;
            prev_pid = pid;
            have_prev = true;
        }

        if(debug) {
            cout << "[md-Fulgor][Debug] Unique blocks: " << block_stats.size() << "\n";
            cout << "[md-Fulgor][Debug] Block transitions: " << block_changes;
            if(total_runs > 1) {
                double pct = (static_cast<double>(block_changes) / static_cast<double>(total_runs - 1)) * 100.0;
                cout << " (~" << fixed << setprecision(2) << pct << "% of consecutive runs)";
            }
            cout << "\n";

            cout << "[md-Fulgor][Debug] Total runs: " << total_runs
                 << ", total kmers spanned: " << total_kmers << "\n";

            if(total_runs > 0) {
                double avg_run = static_cast<double>(sum_run) / static_cast<double>(total_runs);
                cout << "[md-Fulgor][Debug] Run length min/avg/max: " << min_run << " / "
                     << fixed << setprecision(2) << avg_run << " / " << max_run << "\n";
            }

            if(pid_delta_count > 0) {
                double avg_delta = static_cast<double>(pid_delta_sum) / static_cast<double>(pid_delta_count);
                cout << "[md-Fulgor][Debug] Avg consecutive PID delta: " << fixed << setprecision(2) << avg_delta << "\n";
            }

            vector<pair<uint32_t, BlockStats>> top_blocks;
            top_blocks.reserve(block_stats.size());
            for(auto &kv : block_stats)
                top_blocks.emplace_back(kv.first, kv.second);

            sort(top_blocks.begin(), top_blocks.end(), [](const auto &a, const auto &b){
                return a.second.runs > b.second.runs;
            });

            size_t top_limit = min<size_t>(5, top_blocks.size());
            if(top_limit > 0) {
                cout << "[md-Fulgor][Debug] Top blocks by runs:\n";
                for(size_t i = 0; i < top_limit; i++) {
                    cout << "   #" << (i + 1) << " block " << top_blocks[i].first
                         << " -> runs: " << top_blocks[i].second.runs
                         << ", kmers: " << top_blocks[i].second.kmers << "\n";
                }
            }

            if(!samples.empty()) {
                cout << "[md-Fulgor][Debug] Sample runs (idx block pid len):\n";
                for(auto &sample : samples) {
                    cout << "   " << get<0>(sample) << " " << get<1>(sample)
                         << " " << get<2>(sample) << " " << get<3>(sample) << "\n";
                }
            }

            cout.flags(old_flags);
            cout.precision(old_precision);
        }
    } 
    else {
        // Fallback to the classic text-based RLE output (for testing/debugging)
        for(size_t i = 0; i < values.size(); i++) {
            colors_file << values[i];
            if(counts[i] != 1)
                colors_file << ":" << counts[i];
            colors_file << "\n";
        }
    }
    colors_file.close();
}

size_t ColorGraph::get_num_run() {
    return counts.size();
}

double ColorGraph::get_average_run() {
    size_t sum = 0;
    for(auto v: counts)
        sum += v;
    return static_cast<double>(sum) / static_cast<double>(values.size());
}

Node::Node(std::string sequence, std::vector<color_id_t> colors) {
    this->sequence = std::move(sequence);
    this->colors = std::move(colors);
    this->visited = false;
}

bool Node::is_visited() const { return visited; }
void Node::visit() { visited = true; }
void Node::reverse() {
    sequence = reverse_complement(sequence);
    std::reverse(colors.begin(), colors.end());
}
std::string Node::colors_to_string() {
    string colors_str;
    for(color_id_t color: colors) colors_str += to_string(color) + '\n';
    return colors_str;
}
Node::Node() = default;

Path::Path(node_id_t seed) {
    nodes.push_back(seed);
    orientations.push_back(orientation_t::direct);
}
void Path::extend(oriented_node_t node) {
    nodes.push_back(node.node);
    orientations.push_back(node.orientation);
}

// --- Explicit use of orientation_t:: to avoid shadowing ---
void Path::reverse() {
    std::reverse(nodes.begin(), nodes.end());
    std::reverse(orientations.begin(), orientations.end());
    for(auto &o: orientations) {
        o = (o == orientation_t::direct) ? orientation_t::reverse : orientation_t::direct;
    }
}


PathSortKey ColorGraph::get_path_key(Path& path) {
    // If the map is empty, zero
    if (map_ptr == nullptr || path.length() == 0) return {0};

    // Take the first node and its first color as the key for sorting
    node_id_t first_node_id = path.get_n_node_id(0);
    const auto& colors = nodes[first_node_id].colors;

    if (colors.empty()) return {0};

    // Take the first color of the first node
    // Do we need to consider the orientation?
    // In finalize_path_cover, if the orientation is reverse, the colors are inverted.
    // Here we do a quick estimate by taking the "physical" color 0 or last based on the orientation.

    color_id_t gid;
    if (path.get_n_orientation(0) == orientation_t::reverse) {
        gid = colors.back(); // Se è reverse, il "primo" colore è l'ultimo del vettore
    } else {
        gid = colors.front();
    }

    // Directly return the mapped value (Block + PID)
    return { map_ptr->get_mapped_value(gid) };
}

void ColorGraph::sort_paths() {
    cout << "* Reordering paths (Strategy: First Value Alignment)..." << endl;
    auto start = chrono::steady_clock::now();

    // Vector: <Key, OriginalIndex>
    vector<pair<PathSortKey, size_t>> sort_vec;
    sort_vec.reserve(paths.size());

    for (size_t i = 0; i < paths.size(); i++) {
        sort_vec.push_back({ get_path_key(paths[i]), i });
    }

    // Simple Sorting on the Starting Value
    // This groups all paths that START with the same Block and similar PID
    std::sort(sort_vec.begin(), sort_vec.end(), 
        [](const pair<PathSortKey, size_t>& a, const pair<PathSortKey, size_t>& b) {
            return a.first.start_mapped_value < b.first.start_mapped_value;
        }
    );

    // Reconstruction
    vector<Path> new_paths;
    new_paths.reserve(paths.size());
    for (const auto& p : sort_vec) {
        new_paths.push_back(std::move(paths[p.second]));
    }
    paths = std::move(new_paths);

    auto stop = chrono::steady_clock::now();
    cout << "sort_paths(): Reordered " << paths.size() 
         << " paths in " << duration_cast<seconds>(stop - start).count() << " seconds" << endl;
}

orientation_t Path::get_tail_orientation() { return orientations.back(); }
color_id_t Path::get_tail_node_id() { return nodes.back(); }
size_t Path::length() { return nodes.size(); }
orientation_t Path::get_n_orientation(size_t n) { return orientations.at(n); }
color_id_t Path::get_n_node_id(size_t n) { return nodes.at(n); }