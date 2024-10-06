#ifndef EDS_H
#define EDS_H

#include <iostream>
#include <filesystem>

class EDS{
    //  Accepts an elastic-degenerate string (EDS) in format {ACTAG,AGATA}{,ATCC}{ATTTTAA,AGGCGCC,AGAT}{,ATTT,CGCGC}... => every set is closed in curly brackets, strings in set divided by ',' an a size of empty string is supposed to be 0 => total N is 
    private:
        std::string original_eds;
        std::string l_eds;

        size_t original_input_size;
        size_t modified_input_size = 0;

        unsigned int n;  //  number of sets
        unsigned int n_common;  //  number of common parts
        unsigned int l_common;  //  total length of common parts
        unsigned int N;  //  number of characters
        unsigned int m;  //  number of strings in all sets
        unsigned int n_empty_strings;    //  total number of empty string in EDS

        unsigned int min_l;  //  minimal context
        unsigned int max_l;  //  maximal context
        unsigned int avg_l;  //  maximal context

        //  TODO - phased EDS?

    public:
        EDS(std::istream &is);
        ~EDS();

        int make_l_cart(uint8_t l);    //  transform EDS using cartesian approach l-EDS
        int linearize();    //  returns two strings as result of linearization form of EDS
        int stats();        //  print statistics about EDS
        int save(std::ostream &os); //  store EDS on given output stream
};

#endif //  EDS_H