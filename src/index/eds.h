#ifndef EDS_H
#define EDS_H

#include "Bio_FMi.h"

namespace bio_fmi{


class eds : public Bio_FMi
{
    using Bio_FMi::Bio_FMi;

private:
    size_t i_;  //  ivolve  per block
    size_t j_;  //  ivolve  per block

    std::vector<std::vector<unsigned>> old_set_;    //  result structure
    std::vector<std::vector<unsigned>> new_set_;    //  result structure - help

    int delete_file(std::filesystem::path what);
    void print(std::string s);
    template <typename T>
    void print_v(std::vector<T> v);
    std::vector<std::string> find_cartez(std::vector<std::string> base, std::vector<std::string> add);  //  create all combination of strings in vectors - for input reading 
    unsigned get_lcp(std::string ref,std::string change);   //  get longest common prefix - for input reading
    unsigned get_lcs(std::string ref,std::string change);   //  get longest common suffix - for input reading
    int get_change_possition(unsigned location,unsigned block_number,unsigned change_number,unsigned pre_hash_loc,unsigned pos_hash_loc);       //  get pointer to position to reference string from concatenation of changes (d)
    int parse_block(std::string buffer, std::filesystem::path &input_path);
    int flush(const char* start,size_t size, std::ofstream &out_file);
    int flush_change(unsigned number_of_chars_required, std::ofstream &changes_f, std::vector<unsigned> &tmp);
    int flush_segments(std::ofstream &reference_f, std::ofstream &changes_f);


    int first_context_length_;  //  length of first context - could be smaller than context length
    int last_context_length_;   //  length of last context - could be smaller than context length

    std::string buffer_; 
    char divider_ = '#';

public:

    int search(std::string pattern,bool silent);    //  locate pattern
    int read(std::filesystem::path input_path);    //  read input .eds text file
    int parse(std::filesystem::path input_path,const int BLOCK_SIZE);    //  read input .eds text file
    };
}

#endif //  EDS_H