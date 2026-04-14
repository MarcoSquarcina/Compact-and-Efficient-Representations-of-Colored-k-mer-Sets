//
// Created by enrico on 30/08/24.
//

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char **argv){
    if(argc < 2){
        cerr << "Needed input file!" << endl;
        exit(EXIT_FAILURE);
    }
    string input_file_name = argv[1];
    string output_file_name = input_file_name + ".exploded";
    ifstream rle_file = ifstream(input_file_name);
    ofstream exploded_file = ofstream(output_file_name);
    string line;
    int value = 0;
    int count = 0;
    ulong pos = 0;

    cout << "Reading " << input_file_name << "..." << endl;

    while(rle_file >> line){
        // cout << line << endl;
        pos = line.find(':');
        if(pos == string::npos)
            count = 1;
        else {
            // cout << "found \":\" at " << pos << endl;
            count = stoi(line.substr(pos + 1, string::npos));
        }
        value = stoi(line.substr(0, pos));

        for(int i = 0; i < count; i++){
            exploded_file << value << "\n";
        }
    }

    exploded_file.close();
    rle_file.close();

    cout << output_file_name << " saved." << endl;

    return EXIT_SUCCESS;
}