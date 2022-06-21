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
    public:
        csa_wt<wt_huff<rrr_vector<127>>, 32, 64, text_order_sa_sampling<>> reference_index_;
        csa_wt<wt_huff<rrr_vector<127>>, 32, 64, text_order_sa_sampling<>> changes_index_;
        
        select_support_mcl<> sloc_;
        rank_support_v<> rloc_;
        rank_support_v<> riloc_;

        const unsigned minimal_acceptable_size_ = 6;

        bit_vector iloc_;
        bit_vector loc_;
    
        std::string reference_string_;
        std::string string_of_changes_;

        double total_text_size_;
        double total_index_size_;
        bool is_empty_ = true;

        unsigned number_of_segments_;
        unsigned number_of_changes_;
        
        bool original_text_change_;
        std::string new_original_;

        unsigned context_length_;
        unsigned current_context_length_;

        std::vector<unsigned> start_possitions_;
        std::vector<unsigned> base_position_;
        std::vector<unsigned> change_lengths_;
        std::vector<int> offset_;

        Bio_FMi(unsigned context_length);
        ~Bio_FMi();

        int build();
        int read_patterns(std::filesystem::path patterns_input_file, std::vector<std::string> &patterns);
        int save(std::filesystem::path save_path);
        int load(std::filesystem::path index_path);
        void print();
    };
    
}

#endif //  BIO_FMI_H