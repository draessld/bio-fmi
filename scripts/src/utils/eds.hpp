#ifndef EDS_H
#define EDS_H

#include <iostream>
#include <filesystem>
#include <vector>

struct Change {
    std::string value;
    Change* next;
    
};

class EDS{
    //  Accepts an elastic-degenerate string (EDS) in format {ACTAG,AGATA}{,ATCC}{ATTTTAA,AGGCGCC,AGAT}{,ATTT,CGCGC}... => every set is closed in curly brackets, strings in set divided by ',' an a size of empty string is supposed to be 0 => total N is 
    private:
        // std::string original_eds;
        // std::string l_eds;

    public:
        size_t original_input_size;
        std::vector<std::string> changes;
        std::vector<int> base_position_; //  position where degenerate set starts
        std::vector<int> set_size_;      //  offset of positions in reference and non-reference sequence - for every change stored

        std::vector<unsigned int> lengths;

        unsigned int n;  //  number of sets
        unsigned int n_common;  //  number of common parts
        unsigned int l_common;  //  total length of common parts
        unsigned int N;  //  number of characters
        unsigned int m;  //  number of strings in all sets
        unsigned int n_empty_strings;    //  total number of empty string in EDS

        unsigned int min_l;  //  minimal context
        unsigned int max_l;  //  maximal context
        unsigned int avg_l;  //  average context length


        //  TODO - phased EDS?

        EDS(std::istream &is);
        ~EDS() = default;

        int print();        //  print content of every block
        int stats();        //  print statistics about EDS
        int save(std::ostream &os); //  store EDS on given output stream
};

#endif //  EDS_H