#ifndef ALN_H
#define ALN_H

#include "Bio_FMi.h"

namespace bio_fmi{

class aln : public Bio_FMi
{
    using Bio_FMi::Bio_FMi;
    
private:

    std::vector<std::unordered_set<unsigned>> old_set_;
    std::vector<std::unordered_set<unsigned>> new_set_;

    unsigned get_change_possition(unsigned location);
    int get_position_in_other_sequence(int location, unsigned sequence_number);
    unsigned store_l_context(unsigned r_index);
    void get_changes(std::string &ref, std::string &nonref);

public:
    
    int search(std::string pattern,bool silent);
    int read(std::filesystem::path text_input_file);

};
}

#endif //  ALN_H