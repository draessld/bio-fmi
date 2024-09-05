#ifndef BIO_FMI_H
#define BIO_FMI_H

#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <sys/resource.h>

#include <bits/stdc++.h>

#include <sdsl/suffix_arrays.hpp>

using namespace sdsl;

namespace bio_fmi
{

    class Bio_FMi
    {
    private:

        std::filesystem::path eds_file_; //  input text file

        //  metadata
        std::filesystem::path reference_filepath_; //  metadata file
        std::filesystem::path changes_filepath_;   //  metadata file
        std::filesystem::path index_bed_;          //  path to store index

        //  index structures
        csa_wt<wt_huff<rrr_vector<127>>, 32, 64, text_order_sa_sampling<>> reference_index_; //  wavelet tree FM-index structure for reference string
        csa_wt<wt_huff<rrr_vector<127>>, 32, 64, text_order_sa_sampling<>> changes_index_;   // wavelet tree FM-index structure for concatenation of changes

        select_support_mcl<> sloc_; //  select support structure for bit vector loc (1 on every change start)
        rank_support_v<> rloc_;     //  rank support structure for bit vector loc (1 on every change start)
        rank_support_v<> rtloc_;    //  rank support structure for bit vector loc (1 on every change start)
        rank_support_v<> riloc_;    //  rank support structure for bit vector iloc (1 on every sequence start hash)

        bit_vector iloc_; //  bit bector with one on every change start in concatenation of changes
        bit_vector tloc_; //
        bit_vector loc_;  //  bit bector with one on every sequence start in concatenation of changes

        typedef std::unordered_map<int, std::vector<std::pair<int,std::vector<int>>>> hash_type;
        hash_type old_hash_map_;
        hash_type new_hash_map_;

        std::vector<int> base_position_; //  position where degenerate set starts
        std::vector<int> set_size_;      //  offset of positions in reference and non-reference sequence - for every change stored
        std::vector<int> offset_;             //  offset of positions in reference and non-reference sequence - for every change stored

        //  methods
        int parse_eds(); //  load and parse data from eds
        int save();      //  save index files into folder
        int load();      //  load index files from folder
        void print_hash(const hash_type &hash_map);
        void print(); //  print information about index structures

    public:
        //  statistics
        int context_length_;        //  input context length
        double total_index_size_;

        size_t n;                 //  length of EDS = number of nonempty sets
        size_t N;                 //  total length including the
        size_t total_deg_strings; //  number of strings in the degenerate symbols

        //  methods
        Bio_FMi(std::filesystem::path eds_file, int context_length); //  to create new index
        Bio_FMi(std::filesystem::path index_folder);                      //  to load index from the folder
        ~Bio_FMi();

        int build();               //  construct index - WT FM index structures
        int locate(std::string P); //    locate pattern in the index
        hash_type get_result();
        void print_result(const hash_type &hash_map); //  print information about index structures
        void print_stats();                           //  just basic information about EDS extracted during the building
    };

}

#endif //  BIO_FMI_H