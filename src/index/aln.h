#ifndef ALN_H
#define ALN_H

#include "Bio_FMi.h"

namespace bio_fmi{

class aln : public Bio_FMi
{
    using Bio_FMi::Bio_FMi;
    
private:

    std::vector<std::unordered_set<unsigned>> old_set_; //  result structure
    std::vector<std::unordered_set<unsigned>> new_set_; //  result structure - help

    unsigned get_change_possition(unsigned location);       //  get pointer to position to reference string from concatenation of changes (d)
    int get_position_in_other_sequence(int location, unsigned sequence_number);     //  look for non-changes fields in non-reference sequences on input position
    unsigned store_l_context(unsigned r_index);     //      store l-context of change ( while creating concatenation of changes)
    void get_changes(std::string &ref, std::string &nonref);    //  get all changes between reference and nonreference string in format of concatenation of changes  

public:
    
    int search(std::string pattern,bool silent);    //  locate pattern  
    int read(std::filesystem::path text_input_file);    //  read input .aln text file

};
}

#endif //  ALN_H