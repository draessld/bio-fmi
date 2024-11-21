#ifndef EDS_H
#define EDS_H

#include <iostream>
#include <filesystem>
#include <vector>
#include <random>
#include <algorithm>

class EDS{
    //  Accepts an elastic-degenerate string (EDS) in format {ACTAG,AGATA}{,ATCC}{ATTTTAA,AGGCGCC,AGAT}{,ATTT,CGCGC}... => every set is closed in curly brackets, strings in set divided by ',' an a size of empty string is supposed to be 0 => total N is 
    private:
        // std::string original_eds;
        // std::string l_eds;
        double calculateSize();
        int gen_pattern(std::ostream &os, unsigned int size);
    public:
        size_t original_input_size;

        std::vector<std::string> changes;
        std::vector<int> abs_change_number; // number of change in degenerate symbol
        std::vector<unsigned int> set_size;      
        std::vector<unsigned int> cum_set_size;      //  offset of positions in reference and non-reference sequence - for every change stored
        std::vector<bool> is_ref;  //  bit vector with one on every set with size 1

        std::vector<unsigned int> ref_position; //  for every string keeps its index in terms of reference sequence
        std::vector<unsigned int> base_position; //  for every string keeps its starting position (pointer into original text)
        std::vector<unsigned int> lengths;  //  for every string keeps its length

        unsigned int n;  //  number of sets
        unsigned int total_change_size;  // total number of characters in degenerate symbols  
        unsigned int n_common;  //  total number of characters in common parts 
        unsigned int l_common;  //  total length of common parts
        unsigned int N;  //  number of characters
        unsigned int m;  //  number of strings in all sets
        unsigned int n_empty_strings;    //  total number of empty string in EDS

        unsigned int min_l;  //  minimal context
        unsigned int max_l;  //  maximal context
        unsigned int avg_l;  //  average context length

        bool is_empty = true;

        //  TODO - phased EDS?

        EDS(std::istream &is);
        EDS(){is_empty = true;};
        ~EDS() = default;

        bool empty(){return is_empty;};
        int print();        //  print content of every block
        int stats();        //  print statistics about EDS
        int save(std::ostream &os); //  store EDS on given output stream
        void gen_patterns(std::ostream &os,unsigned int t,unsigned int size);
        void language(std::ostream &os);
        void extract(unsigned int position, std::vector<unsigned int> changes);
};

#endif //  EDS_H