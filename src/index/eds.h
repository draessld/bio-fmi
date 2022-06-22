#ifndef EDS_H
#define EDS_H

#include "Bio_FMi.h"

namespace bio_fmi{


class eds : public Bio_FMi
{
    using Bio_FMi::Bio_FMi;

private:

    std::vector<std::vector<unsigned>> old_set_;    //  result structure
    std::vector<std::vector<unsigned>> new_set_;    //  result structure - help

    std::vector<std::string> find_cartez(std::vector<std::string> base, std::vector<std::string> add);  //  create all combination of strings in vectors - for input reading 
    unsigned get_lcp(std::string ref,std::string change);   //  get longest common prefix - for input reading
    unsigned get_lcs(std::string ref,std::string change);   //  get longest common suffix - for input reading
    int get_change_possition(unsigned location,unsigned block_number,unsigned change_number,unsigned pre_hash_loc,unsigned pos_hash_loc);       //  get pointer to position to reference string from concatenation of changes (d)

    int first_context_length_;  //  length of first context - could be smaller than context length
    int last_context_length_;   //  length of last context - could be smaller than context length
public:

    int search(std::string pattern,bool silent);    //  locate pattern
    int read(std::filesystem::path text_input_file);    //  read input .eds text file
    };
}

#endif //  EDS_H