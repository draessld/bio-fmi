#ifndef EDS_H
#define EDS_H

#include "Bio_FMi.h"

namespace bio_fmi{


class eds : public Bio_FMi
{
    using Bio_FMi::Bio_FMi;

private:

    std::vector<std::vector<unsigned>> old_set_;
    std::vector<std::vector<unsigned>> new_set_;

    std::vector<std::string> find_cartez(std::vector<std::string> base, std::vector<std::string> add);
    int get_difference_from_hash(int difference);
    int get_change_possition(unsigned location,unsigned block_number,unsigned change_number,unsigned pre_hash_loc,unsigned pos_hash_loc);
    unsigned get_lcp(std::string ref,std::string change);
    unsigned get_lcs(std::string ref,std::string change);

    unsigned first_context_length_;
    unsigned last_context_length_;
public:

    int search(std::string pattern,bool silent);
    int read(std::filesystem::path text_input_file);
    };
}

#endif //  EDS_H