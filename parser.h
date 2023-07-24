#ifndef PARSER_H
#define PARSER_H

#include <vector>
#include <string>
#include <bits/stdc++.h>
#include <sys/resource.h>
#include <sdsl/bit_vectors.hpp>

using namespace sdsl;

namespace bio_fmi{

class Parser
{

    private:
        const int BLOCK_SIZE;
        std::string buffer_; 
        unsigned context_length;

        void processBlock(std::vector<char> buffer);
    public:
        std::string reference_string_;
        std::string string_of_changes_;

        bit_vector iloc_;
        bit_vector loc;

        Parser(unsigned context_length);
        ~Parser();

        int parseEDS(std::filesystem::path input_file_path);
};

}

#endif //PARSER_H