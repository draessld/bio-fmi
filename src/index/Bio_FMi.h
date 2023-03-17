#ifndef BIO_FMI_H
#define BIO_FMI_H

#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <sys/resource.h>

#include <bits/stdc++.h>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/wavelet_trees.hpp>

using namespace sdsl;

namespace bio_fmi
{

    class Bio_FMi
    {
    private:

    public:
        csa_wt<wt_huff<rrr_vector<127>>, 32, 64, text_order_sa_sampling<>> reference_index_;    //  wavelet tree FM-index structure for reference string
        csa_wt<wt_huff<rrr_vector<127>>, 32, 64, text_order_sa_sampling<>> changes_index_;  // wavelet tree FM-index structure for concatenation of changes
        
        select_support_mcl<> sloc_;     //  select support structure for bit vector loc (1 on every change start)
        rank_support_v<> rloc_;     //  rank support structure for bit vector loc (1 on every change start)
        rank_support_v<> riloc_;     //  rank support structure for bit vector iloc (1 on every sequence start hash)

        const unsigned minimal_acceptable_size_ = 6;    //  recommended border for pattern size and context length

        bit_vector iloc_;   //  bit bector with one on every change start in concatenation of changes
        bit_vector loc_;    //  bit bector with one on every sequence start in concatenation of changes
    
        std::string reference_string_;  //  reference string
        std::string string_of_changes_; //  concatenation of chnges

        double total_text_size_;
        double total_index_size_;
        bool is_empty_ = true;  //  check if index structure is empty
        bool original_text_change_ = false;

        unsigned number_of_segments_;   //  number of segments 
        unsigned number_of_changes_;    //  number of changes 
        
        std::string new_original_;      //  buffer for transformed EDS string

        unsigned context_length_;   //  input context length
        unsigned current_context_length_;   //  current context length - help improve time when the last chunk is smaller - universal distribution of two last chunk to similar size

        std::vector<unsigned> start_possitions_;    //  positions in concatenation of changes, where new sequence starts
        std::vector<unsigned> base_position_;   //  base postiion of change in reference string
        std::vector<unsigned> change_lengths_;  //  length of each change
        std::vector<int> offset_;   //  offset of positions in reference and non-reference sequence - for every change stored

        Bio_FMi(unsigned context_length);
        ~Bio_FMi();

        int build();    //  construct index - WT FM index structures
        int save(std::filesystem::path save_path);  //  save index files into folder
        int load(std::filesystem::path index_path); //  load index files from folder
        void print();   //  print information about index structures
        int read_patterns(std::filesystem::path patterns_input_file, std::vector<std::string> &patterns);   //  read pattern file in specific format
    };
    
}

#endif //  BIO_FMI_H