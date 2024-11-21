#ifndef BIO_FMI_H
#define BIO_FMI_H

#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <sys/resource.h>

#include <bits/stdc++.h>

#include <sdsl/suffix_arrays.hpp>
#include "eds.hpp"

using namespace sdsl;

namespace bio_fmi
{

    class Bio_FMi
    {
    public:

        typedef std::unordered_map<int, std::vector<std::pair<int,std::vector<int>>>> hash_type;
        typedef csa_wt<wt_huff<rrr_vector<127>>, 32, 64, text_order_sa_sampling<>> index_type; 
        
        void print(); //  print information about index structures

        //  statistics
        int context_length_;        //  input context length
        double total_index_size_;
        EDS eds;

        unsigned int n;                 //  length of EDS = number of nonempty sets
        unsigned int N;                 //  total length including the
        unsigned int m;                 //  number of strings in all sets

        //  methods
        Bio_FMi(EDS eds, int context_length); //  to create new index
        Bio_FMi(std::filesystem::path eds_file, int context_length); //  to create new index
        Bio_FMi(std::filesystem::path index_folder);                      //  to load index from the folder
        ~Bio_FMi();

        

        int build();               //  construct index - WT FM index structures
        int locate(const std::string& P); //    locate pattern in the index
        hash_type get_result();
        void print_result(const hash_type &hash_map); //  print information about index structures
        void print_stats();    

    private:
        hash_type old_hash_map_;
        hash_type new_hash_map_;
        void print_hash(const hash_type &hash_map);

        //  index structures
        index_type reference_index_; //  wavelet tree FM-index structure for reference string
        index_type changes_index_;   // wavelet tree FM-index structure for concatenation of changes

        std::filesystem::path eds_file_; //  input text file

        //  metadata
        std::filesystem::path reference_filepath_; //  metadata file
        std::filesystem::path changes_filepath_;   //  metadata file
        std::filesystem::path index_bed_;          //  path to store index


        select_support_mcl<> sloc_; //  select support structure for bit vector loc (1 on every change start)
        rank_support_v<> rloc_;     //  rank support structure for bit vector loc (1 on every change start)
        rank_support_v<> rtloc_;    //  rank support structure for bit vector loc (1 on every change start)
        rank_support_v<> riloc_;    //  rank support structure for bit vector iloc (1 on every sequence start hash)

        bit_vector iloc_; //  bit bector with one on every change start in concatenation of changes
        bit_vector tloc_; //
        bit_vector loc_;  //  bit bector with one on every sequence start in concatenation of changes


        std::vector<int> base_position_; //  position where degenerate set starts
        std::vector<int> set_size_;      //  offset of positions in reference and non-reference sequence - for every change stored
        std::vector<int> offset_;             //  offset of positions in reference and non-reference sequence - for every change stored

        //  methods
        int parse_eds(); //  load and parse data from eds
        int save();      //  save index files into folder
        int load();      //  load index files from folder
};

}

#endif //  BIO_FMI_H