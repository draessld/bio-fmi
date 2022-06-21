#include "aln.h"

namespace bio_fmi
{
    unsigned aln::store_l_context(unsigned r_index)
    {
        std::string part;
        part = (r_index < context_length_) ? reference_string_.substr(0, r_index) : reference_string_.substr(r_index - (context_length_ - 1), context_length_ - 1);

        if (part.size() > 0)
        {
            string_of_changes_ += part;
            return part.size();
        }
        return 0;
    }

    void aln::get_changes(std::string &ref, std::string &nonref)
    {
        nonref += '$';
        unsigned r_index = 0;
        unsigned s_index = 0;
        int stored_r_context = 0;
        bool is_change_open = false;
        unsigned change_size = 0;

        for (size_t p = 0; p < ref.size(); p++)
        {
            /* INSERT   */
            if ((ref[p] == '-') && (nonref[p] != '-'))
            {
                stored_r_context = 0;
                if (is_change_open)
                {
                    /* l-context is stored */
                    string_of_changes_ += nonref[p];
                    change_size++;
                }
                else
                {
                    /* store l-context */
                    change_size += store_l_context(r_index);
                    string_of_changes_ += nonref[p];
                    change_size++;
                    is_change_open = true;
                    base_position_.push_back(0);
                }
                s_index++;
            }

            /*  DELETE  */
            if ((ref[p] != '-') && (nonref[p] == '-'))
            {
                stored_r_context = 0;
                if (!is_change_open)
                {
                    /* store l-context */
                    change_size += store_l_context(r_index);
                    is_change_open = true;
                    base_position_.push_back(0);
                }
                r_index++;
            }

            /*  UPDATE  */
            if ((ref[p] != '-') && (nonref[p] != '-') && (ref[p] != nonref[p]))
            {
                stored_r_context = 0;
                if (is_change_open)
                {
                    string_of_changes_ += nonref[p];
                    change_size++;
                }
                else
                {
                    /* store l-context */
                    change_size += store_l_context(r_index);
                    string_of_changes_ += nonref[p];
                    change_size++;
                    is_change_open = true;
                    base_position_.push_back(0);
                }
                r_index++;
                s_index++;
            }

            /*  store r-context ? */
            if ((ref[p] != '-') && (nonref[p] != '-') && (ref[p] == nonref[p]))
            {
                r_index++;
                s_index++;
                if (is_change_open)
                {
                    if (stored_r_context == 0)
                    {
                        base_position_[number_of_changes_] = r_index - 1;
                    }

                    string_of_changes_ += ref[p];
                    change_size++;
                    stored_r_context++;
                }
            }

            /* close change ? */
            if ((((stored_r_context - context_length_ + 1) == 0) || (ref[p] == '$')) && (is_change_open))
            {
                string_of_changes_ += '#';
                // base_position_.push_back(r_index - 1);
                offset_.push_back(r_index - s_index);
                change_lengths_.push_back(change_size-stored_r_context);
                is_change_open = false;
                number_of_changes_++;
                stored_r_context = 0;
                change_size = 0;
            }
        }
    }

    unsigned aln::get_change_possition(unsigned location)
    {
        /*  Get the number of change  ->    Rank_loc(i)*/
        unsigned change_number = rloc_(location);
        int difference = location - sloc_(change_number);
        std::cout << "location:" << location << "change number:" <<change_number << ", difference:" << difference<< std::endl;

        return base_position_[change_number - 1] - offset_[change_number - 1] - change_lengths_[change_number - 1] + difference - 1;
    }

    int aln::get_position_in_other_sequence(int location, unsigned sequence_number)
    {
        int offset = 0;
        int changes_start = rloc_(start_possitions_[sequence_number - 1]);
        int changes_end = rloc_(start_possitions_[sequence_number]) - 1;
        int precending_change_index = std::upper_bound(base_position_.begin() + changes_start, base_position_.begin() + changes_end, location) - base_position_.begin() - 1;

        std::cout <<", start,end" << changes_start << changes_end<< std::endl;
        std::cout << "pre index change " << precending_change_index << std::endl;

        // std::cout << "offset" << offset << "loc" << (location + current_context_length_) << "bool" << ((location + current_context_length_) < offset) << std::endl;

        /*  error when the location lays in half in change left context*/
        if ((location + current_context_length_) + change_lengths_[changes_start] < base_position_[changes_start])
        {
            /*  in first change  */
            std::cout << "first" << std::endl;
            return location;
        }

        if (location >= base_position_[changes_end])
        {
            /*  location is after all changes  */
            std::cout << "after" << std::endl;
            return location - offset_[changes_end];
        }

        // std::cout << "offset" << (location + current_context_length_) << "loc" << (base_position_[precending_change_index+1]) << "bool" << (location + current_context_length_ < base_position_[precending_change_index+1]) << std::endl;
        // offset = int(base_position_[precending_change_index+1]) - int(change_lengths_[precending_change_index+1]-context_length_+1);

        /*  is location between changes exclude their context?  */
        if ((location >= base_position_[precending_change_index]) && (location + current_context_length_ < base_position_[precending_change_index + 1]))
        {
            return location - offset_[precending_change_index];
        }

        if ((precending_change_index < changes_start) && (location + current_context_length_ < base_position_[precending_change_index + 1]))
        {
            return location;
        }

        // std::cout << "start_change_pos:" << start_change_position << ". offset:" << offset_[precending_change_index] << ", curr cl:" << current_context_length_ << std::endl;
        // if ((precending_change_index == -1) && ((location+current_context_length_) < start_change_position))
        // return location;
        // int res = ((location+current_context_length_ < start_change_position) || (location >= base_position_[precending_change_index]))?location - offset:-1;
        // int res = ((start_change_position + offset + current_context_length_-1) <= location) ? -1 : location - offset;
        // std::cout << "resL" << res << std::endl;
        return -1;
    }

    int aln::search(std::string pattern, bool silent)
    {
        try
        {
            old_set_ = std::vector<std::unordered_set<unsigned>>(number_of_segments_);
            new_set_ = std::vector<std::unordered_set<unsigned>>(number_of_segments_);

            unsigned sequence_number;
            int other_position;
            unsigned chunk_position_offset;
            unsigned chunk_index;
            std::string chunk;
            current_context_length_ = context_length_;
            unsigned pre_last_context_length = context_length_;
            unsigned last_context_length = context_length_;

            if (!silent)
            {

                if ((context_length_ <= minimal_acceptable_size_) || (pattern.size() <= minimal_acceptable_size_))
                    std::cout << "WARNING: minimal acceptable size: " << minimal_acceptable_size_ << " -> ineffective run" << std::endl;
            }

            /*  check pattern length and set the */
            if ((pattern.size() % context_length_) != 0)
            {
                pre_last_context_length = ceil((double)((pattern.size() % context_length_) + context_length_) / 2);
                last_context_length = floor((double)((pattern.size() % context_length_) + context_length_) / 2);
            }

            /*  for each chunk */
            for (chunk_index = 0; chunk_index <= std::round(pattern.size() / context_length_); chunk_index++)
            {

                if (pattern.size() <= context_length_){
                    current_context_length_ = pattern.size();
                    chunk_position_offset = chunk_index * context_length_;
                }
                else if (chunk_index == std::round(pattern.size() / context_length_) - 1)
                { //  pre-last chunk
                    current_context_length_ = pre_last_context_length;
                    chunk_position_offset = chunk_index * context_length_;
                }
                else if (chunk_index == std::round(pattern.size() / context_length_))
                { //  last chunk
                std::cout << "last" << std::endl;
                    current_context_length_ = last_context_length;
                    if (chunk_index > 1)
                        chunk_position_offset = (chunk_index - 1) * context_length_ + pre_last_context_length;
                    else if (chunk_index > 0)
                        chunk_position_offset = pre_last_context_length;
                    else
                        chunk_position_offset = 0;
                }
                else
                {
                    chunk_position_offset = chunk_index * context_length_;
                }
                // chunk_position_offset = chunk_index * current_context_length;

                chunk = pattern.substr(chunk_position_offset, current_context_length_);
                if (chunk.empty())
                    break;

                std::cout << context_length_ << ", " << current_context_length_ << std::endl;
                std::cout << chunk_index << ", " << chunk_position_offset << " -> " << chunk << std::endl;

                /*  SEARCH in changes */
                auto change_locations = sdsl::locate(changes_index_, chunk);
                std::cout << "in Id: " << change_locations.size() << std::endl;
                for (auto location : change_locations)
                {
                    /*  Get the number of sequence  ->  Rank_iloc(i)    */
                    sequence_number = riloc_(location);
                    other_position = get_change_possition(location);
                    std::cout << "location:" << location << ", sequence number:" << sequence_number << ", other position:" << other_position << ", chunk position offset" << chunk_position_offset << std::endl;

                    if (chunk_index == 0)
                    {
                        /*  SAVE in HASH table */
                        new_set_[sequence_number].insert(other_position);
                    }
                    else
                    {
                        auto it = old_set_[sequence_number].find(other_position - chunk_position_offset);
                        if (it != old_set_[sequence_number].end())
                            /*  SAVE in HASH table */
                            new_set_[sequence_number].insert(other_position - chunk_position_offset);
                    }
                }

                /*  SEARCH in reference */
                auto ref_locations = sdsl::locate(reference_index_, chunk);
                std::cout << "in I0: " << ref_locations.size() << std::endl;
                for (auto location : ref_locations)
                {
                    if (chunk_index == 0)
                    {
                        new_set_[0].insert(location);

                        for (sequence_number = 1; sequence_number < number_of_segments_; sequence_number++)
                        {

                            other_position = get_position_in_other_sequence(location, sequence_number);
                            std::cout << "location:" << location << ", sequence number:" << sequence_number << ", other position:" << other_position << ", chunk position offset" << chunk_position_offset << std::endl;
                            if (other_position != -1)
                                new_set_[sequence_number].insert(other_position);
                        }
                    }
                    else
                    {

                        auto it = old_set_[0].find(location - chunk_position_offset);
                        if (it != old_set_[0].end())
                            new_set_[0].insert(location - chunk_position_offset);

                        for (sequence_number = 1; sequence_number < number_of_segments_; sequence_number++)
                        {
                            other_position = get_position_in_other_sequence(location, sequence_number);
                            std::cout << "location:" << location << ", other position:" << other_position << ", chunk position offset" << chunk_position_offset << std::endl;

                            if (other_position != -1)
                            {
                                auto it = old_set_[sequence_number].find(other_position - chunk_position_offset);
                                if (it != old_set_[sequence_number].end())
                                    new_set_[sequence_number].insert(other_position - chunk_position_offset);
                            }
                        }
                    }
                }

                std::swap(old_set_, new_set_);
                new_set_ = std::vector<std::unordered_set<unsigned>>(number_of_segments_);

                //     for (size_t i = 0; i < number_of_segments_; i++)
                // {
                //     if (!silent)
                //     {
                //         /*  print results   */
                //         std::cout << "  {" << i << "}: ";
                //         for (auto loc : old_set_[i])
                //         {
                //             std::cout << loc << ',';
                //         }
                //         std::cout << std::endl;
                //     }
                // }
            }

            unsigned pattern_locations = 0;
            if (!silent)
                std::cout << pattern << "=> " << std::endl;
            for (size_t i = 0; i < number_of_segments_; i++)
            {
                if (!silent)
                {
                    /*  print results   */
                    std::cout << "  {" << i << "}: ";
                    for (auto loc : old_set_[i])
                    {
                        std::cout << loc << ',';
                    }
                    std::cout << std::endl;
                }
                pattern_locations += old_set_[i].size();
            }
            if (!silent)
                std::cout << "  total: " << pattern_locations << std::endl;

            return pattern_locations;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return -1;
        }

        return -1;
    }

    int aln::read(std::filesystem::path text_input_file)
    {

        number_of_segments_ = 0;
        number_of_changes_ = 0;

        /*  open file   */
        std::ifstream input_file(text_input_file);
        if (!input_file.is_open())
        {
            return -1;
        }

        std::vector<std::string> sequences_buffer;
        for (std::string line; std::getline(input_file, line); /**/)
        {
            if (!line.empty())
            {
                boost::trim(line);
                sequences_buffer.push_back(line);
                total_text_size_ += line.size();
            }
        }

        reference_string_ = sequences_buffer[0];
        sequences_buffer[0] += '$';
        string_of_changes_ += '#';
        number_of_segments_ = sequences_buffer.size();
        start_possitions_.push_back(0);

        boost::erase_all(reference_string_, "-");

        for (auto it = begin(sequences_buffer) + 1; it != end(sequences_buffer); ++it)
        {
            if (!(*it).empty())
            {
                get_changes(sequences_buffer[0], *it);
                start_possitions_.push_back(string_of_changes_.size() - 1);
            }
        }

        /* loc and iloc bitvectors */
        loc_ = bit_vector(string_of_changes_.size(), 0);
        iloc_ = bit_vector(string_of_changes_.size(), 0);

        for (size_t i = 0; i < string_of_changes_.size(); i++)
            if (string_of_changes_[i] == '#')
                loc_[i] = 1;

        for (auto i : start_possitions_)
            iloc_[i] = 1;

        // std::cout << "changes: " << string_of_changes_ <<  std::endl;
        return 0;
    }
}
