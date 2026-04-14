//
// Created by enrico on 29/12/22.
// Modified by Marco Squarcina on 14/04/26.
//

#ifndef USTAR_CONSTS_H
#define USTAR_CONSTS_H

#include <map>

#define VERSION "3.0"
#define RLE_SEPARATOR ':'

using namespace std;

enum class seeding_method_t{
    FIRST,
    RANDOM,
    LOWER_MEDIAN_COLOR,
    LOWER_AVERAGE_COLOR,
    HIGHER_AVERAGE_COLOR,
    SIMILAR_COLORS,
    BIGGER_LENGTH,
    SMALLER_LENGTH,
    MORE_CONNECTED,
    LESS_CONNECTED,
    LESS_UNBALANCED,
    MORE_UNBALANCED
};

const map<string, seeding_method_t> seeding_method_names = {
        {"f", seeding_method_t::FIRST},
        {"r", seeding_method_t::RANDOM},
        {"-ma", seeding_method_t::LOWER_MEDIAN_COLOR},
        {"-aa", seeding_method_t::LOWER_AVERAGE_COLOR},
        {"+aa", seeding_method_t::HIGHER_AVERAGE_COLOR},
        {"=a", seeding_method_t::SIMILAR_COLORS},
        {"-l", seeding_method_t::SMALLER_LENGTH},
        {"+l", seeding_method_t::BIGGER_LENGTH},
        {"-c", seeding_method_t::LESS_CONNECTED},
        {"+c", seeding_method_t::MORE_CONNECTED},
        {"-u", seeding_method_t::LESS_UNBALANCED},
        {"+u", seeding_method_t::MORE_UNBALANCED}
};

enum class extending_method_t{
    FIRST,
    RANDOM,
    SIMILAR_COLOR,
    SIMILAR_MEDIAN_COLOR,
    LOWER_MEDIAN_COLOR,
    BIGGER_LENGTH,
    SMALLER_LENGTH,
    MORE_CONNECTED,
    LESS_CONNECTED,
    SIMILAR_COLOR1,
    SAME_BLOCK 
};

const map<string, extending_method_t> extending_method_names = {
        {"f", extending_method_t::FIRST},
        {"r", extending_method_t::RANDOM},
        {"=a", extending_method_t::SIMILAR_COLOR},
        {"=A", extending_method_t::SIMILAR_COLOR1},
        {"=ma", extending_method_t::SIMILAR_MEDIAN_COLOR},
        {"-l", extending_method_t::SMALLER_LENGTH},
        {"-ma", extending_method_t::LOWER_MEDIAN_COLOR},
        {"+l", extending_method_t::BIGGER_LENGTH},
        {"-c", extending_method_t::LESS_CONNECTED},
        {"+c", extending_method_t::MORE_CONNECTED},
        {"=b", extending_method_t::SAME_BLOCK} 
};

enum class postprocess_t{
    NONE,           // do nothing
    AVG,            // sort colors and sequences based on average colors
    AVG_FLIP,       // do AVG and flip colors and rc sequences if necessary
    AVG_FB,         // sort colors and sequences based on average of first and back of colors vector
    AVG_FB_FLIP     // do AVG_FB and flip colors and rc sequences if necessary
};

const map<string, postprocess_t> postprocess_method_names = {
        {"none", postprocess_t::NONE},
        {"avg", postprocess_t::AVG},
        {"avg_flip", postprocess_t::AVG_FLIP},
        {"avg_fb", postprocess_t::AVG_FB},
        {"avg_fb_flip", postprocess_t::AVG_FB_FLIP}
};

enum class encoding_t{
    PLAIN,
    RLE,
    AVG_RLE,
    FLIP,
    FLIP_RLE,
    AVG_FLIP_RLE,
    BINARY,
    BWT,
    OPT_RLE
};
// ----------------------------------------
const map<encoding_t, string> encoding_suffixes = {
        {encoding_t::PLAIN, ""},
        {encoding_t::RLE, ".rle"},
        {encoding_t::FLIP, ".flip"},
        {encoding_t::AVG_RLE, ".avg_rle"},
        {encoding_t::FLIP_RLE, ".flip_rle"},
        {encoding_t::AVG_FLIP_RLE, ".avg_flip_rle"},
        {encoding_t::BINARY, ".bin"},
        {encoding_t::BWT, ".bwt"},
        {encoding_t::OPT_RLE, ".opt_rle"}
};

const map<string, encoding_t> encoding_names = {
        {"plain", encoding_t::PLAIN},
        {"flip", encoding_t::FLIP},
        {"rle", encoding_t::RLE},
        {"avg_rle", encoding_t::AVG_RLE},
        {"flip_rle", encoding_t::FLIP_RLE},
        {"avg_flip_rle", encoding_t::AVG_FLIP_RLE},
        {"bin", encoding_t::BINARY},
        {"bwt", encoding_t::BWT},
        {"opt_rle", encoding_t::OPT_RLE}
};

// -----------------------------------------------------
#endif //USTAR_CONSTS_H