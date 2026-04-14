//
// Created by enrico on 30/09/24.
//

#include <getopt.h>
#include <iostream>
#include "ColorGraph.h"

using namespace std;

struct params_t{
    string sequences_file_name;
    string colors_file_name;
    string output_file_name;
    int kmer_length = 0;
};

void print_help(){
    cout << "Usage: ./optimize-rle -k <kmer-length> -s <sequence-file> -c <color-file> -o <output-base-name>" << endl;
}

params_t parse_args(int argc, char **argv){
    params_t params;

    int choice;
    while((choice = getopt(argc, argv, "s:c:k:o:")) != -1){
        switch(choice){
            case 's':
                params.sequences_file_name = optarg;
                break;
            case 'c':
                params.colors_file_name = optarg;
                break;
            case 'k':
                params.kmer_length = stoi(optarg);
                break;
            case 'o':
                params.output_file_name = optarg;
                break;
            case '?':
                cerr << "Error parse_args(): Missing argument" << endl;
                print_help();
                exit(EXIT_FAILURE);
            default:
                cerr << "Error parse_args(): Unknown option " << choice << endl;
                print_help();
                exit(EXIT_FAILURE);
        }
    }
    // checking params
    if(params.output_file_name.empty()){
        params.output_file_name = "out";
        cerr << "WARN: using default output basename: " << params.output_file_name << endl;
    }
    if(params.kmer_length == 0){
        params.kmer_length = 31;
        cerr << "WARN: using default kmer-length: " << params.kmer_length << endl;
    }
    return params;
}

int main(int argc, char **argv){
    params_t params = parse_args(argc, argv);

    ColorGraph cg = ColorGraph(params.sequences_file_name, params.colors_file_name, params.kmer_length);

    //ColorGraph cg = ColorGraph(vector<string>{"AAAC", "AAG", "ATT", "CCCC", "TTTT"},
    //                           vector<vector<color_id_t>>{{0, 1}, {0}, {1}, {0, 1}, {2, 2}});

    cg.write_sequences(params.output_file_name + ".rleo.fa");
    cg.write_colors(params.output_file_name + ".rleo.colors");

    return EXIT_SUCCESS;
}