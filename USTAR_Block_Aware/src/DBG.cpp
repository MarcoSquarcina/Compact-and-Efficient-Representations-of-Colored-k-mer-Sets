//
// Created by enrico on 20/12/22.
//

#include <fstream>
#include <iostream>
#include <cstring>
#include <algorithm>
#include "DBG.h"
#include "commons.h"

size_t DBG::estimate_n_nodes(){
    // minimum BCALM2 entry
    // >0 LN:i:31 ab:Z:2
    // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    const uintmax_t MINIMUM_ENTRY_SIZE = 18 + kmer_size + 2;
    // auto file_size = std::filesystem::file_size(bcalm_file_name);
    auto file_size = 250000000; // at least 250MB
    return file_size / MINIMUM_ENTRY_SIZE;
}

void DBG::parse_bcalm_file() {
    ifstream bcalm_file;
    bcalm_file.open(bcalm_file_name);

    if(!bcalm_file.good()){
        cerr << "parse_bcalm_file(): Can't access file " << bcalm_file_name << endl;
        exit(EXIT_FAILURE);
    }

    // improve vector push_back() time
    nodes.reserve(estimate_n_nodes());
    size_t nodes_cap = nodes.capacity();
    if(debug)
        cout << "estimated number of unitigs: " << estimate_n_nodes() << endl;

    // start parsing two line at a time
    string line;
    while(getline(bcalm_file, line)){
        // escape comments
        if(line[0] == '#')
            continue;

        size_t serial; // BCALM2 serial
        char dyn_line[MAX_LINE_LEN]; // line after id and length

        // make a new node
        node_t node;

        // check if line fits in dyn_line
        if(line.size() > MAX_LINE_LEN){
            cerr << "parse_bcalm_file(): Lines must be smaller than " << MAX_LINE_LEN << " characters! Found " << line.length() << endl;
            exit(EXIT_FAILURE);
        }

        // ------ parse line ------
        // >25 LN:i:32 ab:Z:14 12   L:-:23:+ L:-:104831:+  L:+:22:-

        // check consistency:
        // must have a def-line
        if(line[0] != '>' || line.find("LN:i:") == string::npos || line.find("ab:Z:") == string::npos){
            cerr << "parse_bcalm_file(): Bad formatted input file: no def-line found!" << endl;
            exit(EXIT_FAILURE);
        }

        // find BCALM2 sequence serial and unitig length
        // format: (ignore 1 char) (read 1 integer) (ignore 5 char) (read 1 integer) (read 1 string)
        sscanf(line.c_str(), "%*c %zd %*5c %d %[^\n]s", &serial, &node.length, dyn_line);

        // check consistency:
        // must have progressive IDs
        if(serial != nodes.size()){
            cerr << "parse_bcalm_file(): Bad formatted input file: lines must have progressive IDs!" << endl;
            exit(EXIT_FAILURE);
        }

        // ------ parse abundances ------
        // dyn_line = "ab:Z:14 12   L:-:23:+ L:-:104831:+  L:+:22:-"
        uint32_t sum_abundance = 0;
        char *token = strtok(dyn_line + 5, " "); // tokenize abundances
        do{
            uint32_t abundance = atoi(token);
            sum_abundance += abundance;
            node.abundances.push_back(abundance);
            token = strtok(nullptr, " "); // next token
        }while(token != nullptr && token[0] != 'L');
        node.average_abundance = sum_abundance / (double) node.abundances.size();
        node.median_abundance = median(node.abundances);

        // ------ parse arcs ------
        // token = "L:-:23:+ L:-:104831:+  L:+:22:-"
        while(token != nullptr){
            arc_t arc{};
            char s1, s2; // left and right signs
            sscanf(token, "%*2c %c %*c %d %*c %c", &s1, &arc.successor, &s2); // L:-:23:+
            arc.forward = (s1 == '+');
            arc.to_forward = (s2 == '+');
            node.arcs.push_back(arc);
            // next arcs
            token = strtok(nullptr, " ");
        }

        // ------ parse sequence line ------
        // TTGAAGGTAACGGATGTTCTAGTTTTTTCTCTTT}
        if(!getline(bcalm_file, line)){
            cerr << "parse_bcalm_file(): expected a sequence here!" << endl;
            exit(EXIT_FAILURE);
        }

        // get the sequence
        node.unitig = line;

        // check consistency:
        // there must be one count for each kmer!
        if((node.unitig.size() - kmer_size + 1) != node.abundances.size()){
            cerr << "parse_bcalm_file(): Bad formatted input file: wrong number of abundances!" << endl;
            cerr << "parse_bcalm_file(): Also make sure that kmer_size=" << kmer_size << endl;
            exit(EXIT_FAILURE);
        }

        // save the node
        nodes.push_back(node);

        if(debug){
            if(nodes_cap != nodes.capacity()){
                cout << "parse_bcalm_file(): nodes capacity changed!\n";
                nodes_cap = nodes.capacity();
            }
        }
    }
    nodes.shrink_to_fit();
    bcalm_file.close();
}

void DBG::parse_ggcat_file() {
    ifstream bcalm_file;
    bcalm_file.open(bcalm_file_name);

    if(!bcalm_file.good()){
        cerr << "parse_ggcat_file(): Can't access file " << bcalm_file_name << endl;
        exit(EXIT_FAILURE);
    }

    // improve vector push_back() time
    nodes.reserve(estimate_n_nodes());
    size_t nodes_cap = nodes.capacity();
    if(debug)
        cout << "estimated number of unitigs: " << estimate_n_nodes() << endl;

    // start parsing two line at a time
    string line;
    while(getline(bcalm_file, line)){
        // escape comments
        if(line[0] == '#')
            continue;

        size_t serial; // BCALM2 serial
        char dyn_line[MAX_LINE_LEN]; // line after id and length

        // make a new node
        node_t node;

        // check if line fits in dyn_line
        if(line.size() > MAX_LINE_LEN){
            cerr << "parse_ggcat_file(): Lines must be smaller than " << MAX_LINE_LEN << " characters! Found " << line.length() << endl;
            exit(EXIT_FAILURE);
        }

        // ------ parse line ------
        // >4 LN:i:21 C:0:1 L:+:0:- L:+:1:+ L:-:2:-

        // check consistency:
        // must have a def-line
        if(line[0] != '>' || line.find("LN:i:") == string::npos || line.find("C:") == string::npos){
            cerr << "parse_ggcat_file(): Bad formatted input file: no def-line found!" << endl;
            exit(EXIT_FAILURE);
        }

        // find BCALM2 sequence serial and unitig length
        // format: (ignore 1 char) (read 1 integer) (ignore 5 char) (read 1 integer) (read 1 string until new line)
        sscanf(line.c_str(), "%*c %zd %*5c %d %[^\n]s", &serial, &node.length, dyn_line);

        // check consistency:
        // must have progressive IDs
        if(serial != nodes.size()){
            cerr << "parse_ggcat_file(): Bad formatted input file: lines must have progressive IDs!" << endl;
            exit(EXIT_FAILURE);
        }

        if(debug)
            cout << "node " << serial << ", size " << node.length << "\n";

        // ------ parse colors ------
        // dyn_line = "C:0:1 L:+:0:- L:+:1:+ L:-:2:-"
        char *token = strtok(dyn_line, " "); // tokenize colors

        uint32_t sum_colors = 0;
        while(token != nullptr && token[0] != 'L'){
            int count, color; // left and right signs
            sscanf(token, "%*2c %x %*c %d", &color, &count); // C:0:1

            if(debug)
                cout << "\tcolor " << color << ", count " << count << "\n";

            // decode RLE
            for(int i = 0; i < count; i++) {
                node.colors.push_back(color);
                sum_colors += sum_colors + color;
            }

            // next color
            token = strtok(nullptr, " ");
        }

        node.average_color = sum_colors / (double) node.colors.size();
        node.median_color = median(node.colors);

        // ------ parse arcs ------
        // token = "L:+:0:- L:+:1:+ L:-:2:-"
        while(token != nullptr){
            arc_t arc{};
            char s1, s2; // left and right signs
            sscanf(token, "%*2c %c %*c %d %*c %c", &s1, &arc.successor, &s2); // L:+:0:-
            arc.forward = (s1 == '+');
            arc.to_forward = (s2 == '+');
            node.arcs.push_back(arc);
            // next arcs
            token = strtok(nullptr, " ");
        }

        // ------ parse sequence line ------
        // TTGAAGGTAACGGATGTTCTAGTTTTTTCTCTTT
        if(!getline(bcalm_file, line)){
            cerr << "parse_ggcat_file(): expected a sequence here!" << endl;
            exit(EXIT_FAILURE);
        }

        // get the sequence
        node.unitig = line;

        // check consistency:
        // there must be one color for each kmer!
        if((node.unitig.size() - kmer_size + 1) != node.colors.size()){
            cerr << "parse_ggcat_file(): Bad formatted input file: wrong number of colors in line " << serial << "!" << endl;
            cerr << "parse_ggcat_file(): Also make sure that kmer_size=" << kmer_size << endl;
            exit(EXIT_FAILURE);
        }

        // save the node
        nodes.push_back(node);

        if(debug){
            if(nodes_cap != nodes.capacity()){
                cout << "parse_ggcat_file(): nodes capacity changed!\n";
                nodes_cap = nodes.capacity();
            }
        }
    }
    nodes.shrink_to_fit();
    bcalm_file.close();
}

DBG::DBG(const string &bcalm_file_name, uint32_t kmer_size, bool debug){
    this->bcalm_file_name = bcalm_file_name;
    this->kmer_size = kmer_size;
    this->debug = debug;

    // build the graph
    parse_ggcat_file();

    // compute graph parameters
    compute_graph_parameters();
}

DBG::~DBG(){
    nodes.clear();
}

DBG::DBG(const vector<string> &simplitigs, const vector<vector<uint32_t>> &simplitigs_colors, bool debug){
    this->bcalm_file_name = "NULL";

    // induced by parameters
    this->kmer_size = simplitigs[0].length() - simplitigs_colors[0].size() + 1;

    this->debug = debug;

    // build the graph
    build_colors_graph(simplitigs, simplitigs_colors);

    // compute graph parameters
    compute_graph_parameters();
}

void DBG::build_colors_graph(const vector<string> &simplitigs, const vector<vector<uint32_t>> &simplitigs_colors){
    // build incident matrix
    map<uint32_t,vector<uint32_t>> color_map_front;
    map<uint32_t,vector<uint32_t>> color_map_back;
    size_t idx = 0;
    cout << "Building edge maps..." << endl;
    for(const auto &node_colors: simplitigs_colors){
        color_map_front[node_colors.front()].push_back(idx);
        color_map_back[node_colors.back()].push_back(idx);
        idx++;
    }

    cout << "Building adjacency list..." << endl;
    nodes.reserve(simplitigs.size());
    for(idx = 0; idx < simplitigs.size(); idx++){
        node_t node;

        node.unitig = simplitigs[idx];
        node.length = simplitigs[idx].length();
        node.colors = simplitigs_colors[idx];

        auto back_color = simplitigs_colors[idx].back();    // -------C
        auto front_color = simplitigs_colors[idx].front();  // C-------

        // arcs +:-
        for(auto succ: color_map_back[back_color]){
            arc_t arc{};
            arc.successor = succ;
            arc.forward = true;
            arc.to_forward = false;
            node.arcs.push_back(arc);
        }
        // arcs +:+
        for(auto succ: color_map_front[back_color]){
            arc_t arc{};
            arc.successor = succ;
            arc.forward = true;
            arc.to_forward = true;
            node.arcs.push_back(arc);
        }
        // arcs -:-
        for(auto succ: color_map_back[front_color]){
            arc_t arc{};
            arc.successor = succ;
            arc.forward = false;
            arc.to_forward = false;
            node.arcs.push_back(arc);
        }
        // arcs -:+
        for(auto succ: color_map_front[front_color]){
            arc_t arc{};
            arc.successor = succ;
            arc.forward = false;
            arc.to_forward = true;
            node.arcs.push_back(arc);
        }
        nodes.push_back(node);
    }
    cout << "Done." << endl;

}

void DBG::compute_graph_parameters(){
    size_t sum_unitig_length = 0;
    double sum_colors = 0;
    for(const auto &node : nodes) {
        n_arcs += node.arcs.size();
        // n_kmers += node.abundances.size();
        n_kmers += node.colors.size();
        sum_unitig_length += node.length;
        sum_colors += node.average_color * (double) node.colors.size();

        if(node.arcs.empty()) n_iso++;

        size_t n_arcs_in = 0;
        size_t n_arcs_out = 0;
        for(const auto &arc : node.arcs)
            if(arc.forward)
                n_arcs_out++;
            else
                n_arcs_in++;

        if(n_arcs_out == 0)
            n_sinks++;
        if(n_arcs_in == 0)
            n_sources++;
    }
    avg_unitig_len = (double) sum_unitig_length / (double) nodes.size();
    avg_colors = sum_colors / (double) n_kmers;
}

void DBG::print_stat() {
    cout << "\n";
    cout << "DBG stats:\n";
    cout << "   number of kmers:            " << n_kmers << "\n";
    cout << "   number of nodes:            " << nodes.size() << "\n";
    cout << "   number of isolated nodes:   " << n_iso << " (" << double (n_iso) / double (nodes.size()) * 100 << "%)\n";
    cout << "   number of sources:          " << n_sources << " (" << double (n_sources) / double (nodes.size()) * 100 << "%)\n";
    cout << "   number of sinks:            " << n_sinks << " (" << double (n_sinks) / double (nodes.size()) * 100 << "%)\n";
    cout << "   number of arcs:             " << n_arcs << "\n";
    cout << "   graph density:              " << double (n_arcs) / double (8 * nodes.size()) * 100 << "%\n";
    cout << "   average unitig length:      " << avg_unitig_len << "\n";
    cout << "   average color:              " << avg_colors << "\n";
    cout << "\n";
}

bool DBG::verify_overlaps() {
    for(const auto &node : nodes){
        for(const auto &arc : node.arcs)
            if(!overlaps(node, arc))
                return false;
    }
    return true;
}

bool DBG::overlaps(const node_t &node, const arc_t &arcs){
    string u1, u2;
    if (arcs.forward) // + --> +/-
        // last kmer_size - 1 characters
        u1 = node.unitig.substr(node.unitig.length() - kmer_size + 1);
    else // - --> +/-
        // first kmer_size - 1 characters reverse-complemented
        u1 = reverse_complement(node.unitig.substr(0, kmer_size - 1));

    if(arcs.to_forward) // +/- --> +
        // first kmer_size - 1 characters
        u2 = nodes[arcs.successor].unitig.substr(0, kmer_size - 1);
    else  // +/- --> -
        // last kmer_size - 1 characters reverse-complemented
        u2 = reverse_complement(nodes[arcs.successor].unitig.substr(nodes[arcs.successor].unitig.length() - kmer_size + 1));

    return u1 == u2;
}

string DBG::reverse_complement(const string &s) {
    string rc(s);

    for(size_t i = 0; i < s.length(); i++) {
        char c;
        switch (s[i]) {
            case 'A':
            case 'a':
                c = 'T';
                break;
            case 'C':
            case 'c':
                c = 'G';
                break;
            case 'T':
            case 't':
                c = 'A';
                break;
            case 'G':
            case 'g':
                c = 'C';
                break;
            default:
                cerr << "complement(): Unknown nucleotide!" << endl;
                exit(EXIT_FAILURE);
        }
        rc[s.length() - 1 - i] = c;
    }
    return rc;
}

void DBG::to_bcalm_file(const string &file_name) {
    ofstream file;
    file.open(file_name);

    int id = 0;
    for(const auto &node : nodes){
        // >3 LN:i:33 ab:Z:2 2 3    L:+:138996:+
        // CAAAACCAGACATAATAAAAATACTAATTAATG
        file << ">" << id++ << " LN:i:" << node.length << " ab:Z:";
        for(auto &ab : node.abundances)
            file << ab << " ";
        for(auto &arcs : node.arcs)
            file << "L:" << (arcs.forward ? "+" : "-") << ":" << arcs.successor << ":" << (arcs.to_forward ? "+" : "-") << " ";
        file << "\n" << node.unitig << "\n";
    }

    file.close();
}

void DBG::to_ggcat_file(const string &file_name) {
    ofstream file;
    file.open(file_name);

    int id = 0;
    for(const auto &node : nodes){
        if(debug)
            printf("node %d\n", id);

        // >0 LN:i:21 C:1:1 L:+:4:-
        // CGTTTTTTTTTTTTTTTTTTT
        file << ">" << id++ << " LN:i:" << node.length << " ";

        int prev_color = -1;
        int count = 0;
        for(auto &color : node.colors) {
            if (debug)
                printf("\tcolor %d\n", color);
            if (prev_color == color)
                count++;
            else {
                if (prev_color != -1)
                    file << "C:" << count << ":" << prev_color << " ";
                count = 1;
                prev_color = color;
            }
        }
        file << "C:" << prev_color << ":" << count << " ";

        for(auto &arcs : node.arcs)
            file << "L:" << (arcs.forward ? "+" : "-") << ":" << arcs.successor << ":" << (arcs.to_forward ? "+" : "-") << " ";
        file << "\n" << node.unitig << "\n";
    }

    file.close();
}


bool DBG::validate(){
    string fasta_dbg = "unitigs.k"+ to_string(kmer_size) +".ustar.fa";
    to_ggcat_file(fasta_dbg);

    ifstream bcalm_dbg, this_dbg;
    bcalm_dbg.open(bcalm_file_name);
    this_dbg.open(fasta_dbg);

    string tok1, tok2;
    while(bcalm_dbg >> tok1 && this_dbg >> tok2)
        if(tok1 != tok2) {
            cerr << "Files differ here: " << tok1 << " != " <<  tok2.length() << endl;
            return false;
        }
    return true;
}

bool DBG::verify_input(){
    bool good = true;
    if (verify_overlaps())
        cout << "YES! DBG is an overlapping graph!\n";
    else {
        cout << "OOPS! DBG is NOT an overlapping graph\n";
        good = false;
    }
    if(validate())
        cout << "YES! DBG is the same as input file!\n";
    else {
        cout << "OOPS! DBG is NOT the same as input file!\n";
        good = false;
    }
    cout << endl;
    return good;
}

void DBG::get_nodes_from(node_idx_t node, vector<bool> &forwards, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards, const vector<bool> &mask) {
    // vectors must be empty
    to_nodes.clear();
    to_forwards.clear();

    for(auto &arc : nodes.at(node).arcs){
        if(mask.at(arc.successor)) continue;
        to_nodes.push_back(arc.successor);
        forwards.push_back(arc.forward);
        to_forwards.push_back(arc.to_forward);
    }
}

void DBG::get_consistent_nodes_from(node_idx_t node, bool forward, vector<node_idx_t> &to_nodes, vector<bool> &to_forwards, const vector<bool> &mask) {
    // vectors must be empty
    to_nodes.clear();
    to_forwards.clear();

    for(auto &arc : nodes.at(node).arcs){
        if(!mask.empty() && mask.at(arc.successor)) continue;
        if(arc.forward == forward) { // consistent nodes only
            to_nodes.push_back(arc.successor);
            to_forwards.push_back(arc.to_forward);
        }
    }
}

string DBG::spell(const vector<node_idx_t> &path_nodes, const vector<bool> &forwards) {
    if(path_nodes.size() != forwards.size()){
        cerr << "spell(): Inconsistent path!" << endl;
        exit(EXIT_FAILURE);
    }
    if(path_nodes.empty()) {
        cerr << "spell(): You're not allowed to spell an empty path!" << endl;
        exit(EXIT_FAILURE);
    }

    string contig;
    // first node as a seed
    if(forwards[0])
        contig = nodes.at(path_nodes[0]).unitig;
    else
        contig = reverse_complement(nodes.at(path_nodes[0]).unitig);

    // extend the seed
    for(size_t i = 1; i < path_nodes.size(); i++){
        if(forwards[i])
            contig += nodes.at(path_nodes[i]).unitig.substr(kmer_size - 1);
        else {
            string unitig = nodes.at(path_nodes[i]).unitig;
            size_t len = unitig.length() - (kmer_size - 1);
            contig += reverse_complement(unitig.substr(0, len));
        }
    }

    return contig;
}

string DBG::spell(const node_idx_t path_nodes, const bool forwards) {
    string contig;
    // first node as a seed
    if(forwards)
        contig = nodes.at(path_nodes).unitig;
    else
        contig = reverse_complement(nodes.at(path_nodes).unitig);
    return contig;
}

void DBG::get_colors(const node_idx_t path_nodes, const bool forwards, vector<uint32_t> &colors) {
    //          3 5
    // forward: A C T T
    //          5 3
    // rev-com: A A G T
    if (forwards) // read forward
        for(uint32_t color : nodes.at(path_nodes).colors)
            colors.push_back(color);
    else // read backward
        for(int k = int (nodes.at(path_nodes).colors.size() - 1); k > -1; k--)
            colors.push_back(nodes.at(path_nodes).colors[k]);
}

void DBG::get_colors(const vector<node_idx_t> &path_nodes, const vector<bool> &forwards, vector<uint32_t> &colors) {
    //          3 5
    // forward: A C T T
    //          5 3
    // rev-com: A A G T
    for (size_t i = 0; i < path_nodes.size(); i++)
        if (forwards[i]) // read forward
            for(uint32_t color : nodes.at(path_nodes[i]).colors)
                colors.push_back(color);
        else // read backward
            for(int k = int (nodes.at(path_nodes[i]).colors.size() - 1); k > -1; k--)
                colors.push_back(nodes.at(path_nodes[i]).colors[k]);
}

bool DBG::check_path_consistency(const vector<node_idx_t> &path_nodes, const vector<bool> &forwards) {
    // one orientation for each node
    if(path_nodes.size() != forwards.size())
        return false;

    for(size_t i = 0; i < path_nodes.size() - 1; i++){
        bool found = false;
        // is there an arc leading to a consistent node?
        for(auto &arc : nodes.at(path_nodes[i]).arcs)
            // same node orientation and successor check
            if(arc.forward == forwards[i] && arc.successor == path_nodes[i + 1])
                found = true;
        if(!found)
            return false;
    }
    return true;
}

uint32_t DBG::get_n_kmers() const {
    return n_kmers;
}

uint32_t DBG::get_n_nodes() const {
    return nodes.size();
}

const node_t & DBG::get_node(node_idx_t node){
    return nodes.at(node);
}

uint32_t DBG::get_kmer_size() const {
    return kmer_size;
}

const vector<node_t> * DBG::get_nodes() {
    return &nodes;
}

void DBG::clear(){
    nodes.clear();
}
