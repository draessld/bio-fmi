#include "eds.h"

namespace bio_fmi
{

    std::vector<std::string> eds::find_cartez(std::vector<std::string> base, std::vector<std::string> add)
    {
        std::vector<std::string> res;
        for (auto b : base)
            for (auto a : add)
                res.push_back(b + a);

        // for (auto r : res)
        // {
        //     std::cout << r << std::endl;
        // }

        return res;
    }

    unsigned eds::get_lcp(std::string ref, std::string change)
    {
        unsigned i = 0;
        unsigned a = ref.size();
        unsigned b = change.size();
        while (true)
        {
            if ((i < a) && (i < b) && (ref[i] == change[i]))
                i++;
            else
                break;
        }

        return i;
    }

    unsigned eds::get_lcs(std::string ref, std::string change)
    {
        unsigned i = 0;
        unsigned a = ref.size();
        unsigned b = change.size();
        while (true)
        {
            if ((i < a) && (i < b) && (ref[a - i - 1] == change[b - i - 1]))
                i++;
            else
                break;
        }

        return i;
    }

    int eds::get_change_possition(unsigned location, unsigned block_number, unsigned change_number, unsigned pre_hash_loc, unsigned pos_hash_loc)
    {
        // std::cout << block_number << ',' << first_context_length_ << ',' << last_context_length_ << std::endl;
        // std::cout << pre_hash_loc << ',' << pos_hash_loc << ',' << base_position_[change_number] << ',' << offset_[change_number] << std::endl;

        if (block_number == 0 && (pre_hash_loc + current_context_length_) <= (first_context_length_ + 1))
        {
            // std::cout << "pre" << std::endl;
            return -1;
        }
        else if ((block_number - 1) == number_of_segments_ && pos_hash_loc < (last_context_length_ + current_context_length_))
        {
            // std::cout << "pos" << std::endl;
            return -1;
        }
        else
            // if ((((pos_hash_loc+1-context_length_) < location) && (string_of_changes_[pos_hash_loc-1] != '$')) || (((pre_hash_loc+current_context_length_-1) < context_length_) && (block_number != 0))){
            //     std::cout << "other" << std::endl;
            //     return -1;
            // }

            return base_position_[change_number] + offset_[change_number] - pos_hash_loc;
    }

    int eds::search(std::string pattern, bool silent)
    {

        try
        {

            unsigned block_number = 0;
            unsigned change_number = 0;
            int other_position = -1;

            unsigned chunk_position_offset = 0;
            unsigned chunk_index = 0;
            std::string chunk = "";
            current_context_length_ = context_length_;
            unsigned pre_last_context_length = context_length_;
            unsigned last_context_length = context_length_;

            unsigned offset = 0;
            unsigned pre_hash_loc = 0;
            unsigned pos_hash_loc = 0;

            if ((context_length_ <= minimal_acceptable_size_) || (pattern.size() <= minimal_acceptable_size_))
                std::cout << "WARNING: minimal acceptable size: " << minimal_acceptable_size_ << " -> ineffective run" << std::endl;

            /*  check pattern length and set the */
            if ((pattern.size() % context_length_) != 0)
            {
                pre_last_context_length = ceil((double)((pattern.size() % context_length_) + context_length_) / 2);
                last_context_length = floor((double)((pattern.size() % context_length_) + context_length_) / 2);
            }

            /*  for each chunk */
            for (chunk_index = 0; chunk_index <= std::round(pattern.size() / context_length_); chunk_index++)
            {

                if (chunk_index == std::round(pattern.size() / context_length_) - 1)
                { //  pre-last chunk
                    current_context_length_ = pre_last_context_length;
                    chunk_position_offset = chunk_index * context_length_;
                }
                else if (chunk_index == std::round(pattern.size() / context_length_))
                { //  last chunk
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
                // chunk_position_offset = chunk_index * current_context_length_;

                chunk = pattern.substr(chunk_position_offset, current_context_length_);
                if (chunk.empty())
                    break;

                // std::cout << chunk_index << ", " << chunk_position_offset << " -> " << chunk << std::endl;

                // /*  SEARCH in changes */
                auto change_locations = sdsl::locate(changes_index_, chunk);
                // std::cout << "in Id: " << change_locations.size() << std::endl;

                /*  SEARCH in reference */
                auto ref_locations = sdsl::locate(reference_index_, chunk);
                // std::cout << "in I0: " << ref_locations.size() << std::endl;

                for (size_t l = 0; l < (ref_locations.size() + change_locations.size()); l++)
                {
                    if (l >= ref_locations.size())
                    {
                        // std::cout << change_locations[l-ref_locations.size()] << " : CHANGE" << std::endl;

                        other_position = change_locations[l - ref_locations.size()];
                        block_number = riloc_(other_position);
                        change_number = rloc_(other_position);
                        pre_hash_loc = other_position - sloc_(change_number + 1);
                        pos_hash_loc = sloc_(change_number + 2) - other_position;

                        // std::cout << pre_hash_loc << ',' << pos_hash_loc << ',' << base_position_[change_number] << ',' << offset_[change_number] << std::endl;
                        other_position = get_change_possition(other_position, block_number - 1, change_number - 1, other_position - pre_hash_loc, pos_hash_loc - other_position);
                        if (other_position == -1) // do not save change in contexts parts again - found in reference
                            continue;
                        offset = offset_[change_number - 1];
                    }
                    else
                    {
                        // std::cout << ref_locations[l] << " : REF" << std::endl;

                        block_number = 0;
                        change_number = 0;
                        other_position = ref_locations[l];
                        offset = 0;
                    }

                    // std::cout << "Location: " << other_position << ", block_number: " << block_number; // << ", change_number: " << change_number << ", offset: " << offset << ", chunk possition: " << chunk_position_offset << std::endl;

                    if (chunk_index == 0)
                    {
                        new_set_.push_back({other_position, offset, change_number});
                        //     // location_r->sequences.push_back(change_number);
                        //     // std::cout << "here" << std::endl;
                        //     // location_r->location = other_position;
                        //     // location_r->offset = offset;
                        //     // new_set_.push_back(location_r);
                    }
                    else
                    {
                        /*  VALIDATION  */
                        // std::cout << old_set_.size() << std::endl;
                        for (auto loc : old_set_)
                        {

                            // std::cout << loc->location << ", " << loc->offset << ", " << chunk_position_offset << ", " << loc->sequences.back()<< std::endl;
                            if (loc.back() == change_number || change_number == 0)
                            {
                                // std::cout << "a";
                                if ((other_position + loc[1] - offset - chunk_position_offset) == loc[0])
                                {
                                    // std::cout << " - yes" << std::endl;
                                    new_set_.push_back(loc);
                                }
                            }

                            if (loc.back() < change_number)
                            {
                                // std::cout << "b";
                                if ((other_position + loc[1] - chunk_position_offset) == loc[0])
                                {
                                    // std::cout << " - yes" << std::endl;
                                    loc[1] += offset;
                                    loc.push_back(change_number);

                                    new_set_.push_back(loc);
                                }
                            }
                        }
                    }
                }

                std::swap(old_set_, new_set_);
                new_set_.clear();

                /*  print results   */
                // for (auto loc : old_set_)
                // {
                //     std::cout << "      " << loc[0] << ',' << loc[1] << ": [";
                //     for (size_t i = 2; i < loc.size(); i++)
                //     {
                //         std::cout << loc[i] << ',';
                //     }
                //     std::cout << ']' << std::endl;
                // }

                // std::cout << "  total: " << old_set_.size() << std::endl;
            }

            if (!silent)
            {
                for (auto loc : old_set_)
                {
                    std::cout << "      " << loc[0] << ',' << loc[1] << ": [";
                    for (size_t i = 2; i < loc.size(); i++)
                    {
                        std::cout << loc[i] << ',';
                    }
                    std::cout << ']' << std::endl;
                }

                std::cout << "  total: " << old_set_.size() << std::endl;
            }

            return old_set_.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return -1;
        }

        return -1;
    }

    int eds::read(std::filesystem::path text_input_file)
    {
        std::cout << "reading eds" << std::endl;

        number_of_segments_ = 0;
        number_of_changes_ = 0;

        /*  open file   */
        std::ifstream input_file(text_input_file);
        if (!input_file.is_open())
        {
            return -1;
        }

        first_context_length_ = -1;
        bool open_block = false;
        bool open_change = false;

        unsigned pre_ref_length = 0;
        unsigned post_ref_length = 0;

        std::string r_context = "";
        std::string l_context = "";
        unsigned r_index = 0;
        string_of_changes_ = '#';
        unsigned s_index = 1;
        bool save_ref = false;

        unsigned block_end_possition = 0;

        std::vector<std::string> changes_base;
        std::vector<std::string> changes;
        std::string change = "";
        std::string ref_change = "";
        std::vector<unsigned> block_size;

        std::string line;
        getline(input_file, line);
        line += '$';

        for (auto character : line)
        {
            if (character == '{')
            {
                if (open_change)
                {
                    original_text_change_ = true;
                    std::string appending = reference_string_.substr(block_end_possition - 1, r_index - block_end_possition + 1);
                    // std::cout << "appending " << appending << ", starting on: " << block_end_possition << ", length: " << (r_index - block_end_possition) << std::endl;
                    for (auto c : changes)
                    {
                        changes_base.push_back(c + appending);
                        // std::cout << c + appending << std::endl;
                    }
                    // if (pre_ref_length != -1)
                    //     pre_ref_length += appending.size();

                    // post_ref_length = appending.size();
                    changes.clear();
                    r_context.clear();
                }
                else
                {
                    if (r_index >= (context_length_ - 1))
                    {
                        l_context = reference_string_.substr(r_index - context_length_ + 1, context_length_ - 1);
                    }
                    else
                    {
                        l_context = reference_string_.substr(0, r_index);
                    }
                    // std::cout << "setting r-context of " << reference_string_ << ", on "<< r_index << "=" << l_context<< std::endl;
                    start_possitions_.push_back(s_index - 1);
                }
                if (first_context_length_ == -1)
                {
                    first_context_length_ = l_context.size();
                }
                open_block = true;
                open_change = true;
                block_size.push_back(0);
                save_ref = true;
                continue;
            }

            if (character == ',')
            {
                if (!save_ref)
                {
                    changes.push_back(change);
                    change.clear();
                }
                else
                {
                    changes.push_back(change);
                    reference_string_ += change;
                    r_index += change.size();
                    // if (pre_ref_length == -1)
                    //     pre_ref_length = change.size();
                    change.clear();
                    save_ref = false;
                }
                block_size[number_of_segments_]++;
                continue;
            }

            if (character == '}')
            {
                if (block_size[number_of_segments_] != 0)
                {
                    changes.push_back(change);
                    change.clear();
                    if (!changes_base.empty())
                    {
                        // if (post_ref_length != -1)
                        //     post_ref_length += changes[0].size();
                        changes = find_cartez(changes_base, changes);
                    }
                }

                open_block = false;
                block_end_possition = r_index + 1;
                block_size[number_of_segments_]++;
                number_of_segments_++;
                continue;
            }

            if (open_block)
            {
                change += character;
                continue;
            }

            if (open_block && save_ref)
            {
                reference_string_ += character;
                r_index++;
                continue;
            }

            if (!open_block)
            {
                if (open_change)
                    r_context += character;
                else
                {
                    if (character != '$')
                        new_original_ += character;
                }

                reference_string_ += character;
                r_index++;
            }

            if ((open_change && (r_context.size() == (context_length_ - 1))) || character == '$')
            {
                // std::cout << "saving changes" << std::endl;
                if (changes.size() == 0)
                {
                    open_change = false;
                    break;
                }
                new_original_ += '{';
                // std::cout << "r_index: " << r_index << ", reference: " << reference_string_ << std::endl;

                new_original_ += changes[0];
                ref_change = changes[0];
                for (size_t i = 1; i < changes.size(); i++)
                {

                    new_original_ += ',';
                    new_original_ += changes[i];
                    pre_ref_length = get_lcp(ref_change, changes[i]);
                    post_ref_length = get_lcs(ref_change, changes[i]);
                    // std::cout << i << ": " << changes[i] << ", " << changes[0] << std::endl;
                    // std::cout << "ref length: " << pre_ref_length << ',' << post_ref_length << std::endl;
                    // std::cout << l_context << std::endl;
                    // std::cout << r_context << std::endl;

                    offset_.push_back(changes[i].size() - changes[0].size());

                    /*  cut l-context   */
                    if (pre_ref_length >= (context_length_ - 1))
                    {
                        string_of_changes_ += changes[i].substr(pre_ref_length - (context_length_ - 1), changes[i].size() - pre_ref_length + (context_length_ - 1));
                    }
                    else
                    {
                        if (l_context.size() + pre_ref_length < context_length_)
                        {
                            string_of_changes_ += l_context + changes[i];
                        }
                        else
                        {
                            string_of_changes_ += l_context.substr(l_context.size() - ((context_length_ - 1) - pre_ref_length), ((context_length_ - 1) - pre_ref_length)) + changes[i];
                        }
                    }
                    /* cut r-context by */
                    if (post_ref_length >= (context_length_ - 1))
                    {
                        // std::cout << 'd' << std::endl;
                        string_of_changes_ = string_of_changes_.substr(0, string_of_changes_.size() + context_length_ - post_ref_length - 2);
                        base_position_.push_back(r_index - post_ref_length - 1 + offset_.back());
                    }
                    else
                    {
                        if (r_context.size() + post_ref_length < context_length_)
                        {
                            // std::cout << 'e' << std::endl;
                            string_of_changes_ += r_context;
                            base_position_.push_back(r_index);
                        }
                        else
                        {
                            // std::cout << 'f' << std::endl;
                            string_of_changes_ += r_context.substr(0, r_context.size() - post_ref_length);
                            base_position_.push_back(r_index - post_ref_length);
                        }
                    }

                    string_of_changes_ += '#';
                    // std::cout << string_of_changes_ << std::endl;
                    s_index = string_of_changes_.size();
                    number_of_changes_++;
                    // base_position_.push_back(r_index-post_ref_length -1 + (changes[i].size() - changes[0].size()) );
                }
                new_original_ += '}';
                new_original_ += r_context;

                block_size[number_of_segments_] = changes.size();
                changes.clear();
                open_change = false;
                last_context_length_ = r_context.size();
                r_context.clear();
                changes_base.clear();
                pre_ref_length = -1;
            }
        }

        new_original_.pop_back();
        start_possitions_.push_back(s_index - 1);
        if (block_size.back() == 0)
            block_size.pop_back();

        /* loc and iloc bitvectors */
        loc_ = bit_vector(string_of_changes_.size(), 0);
        iloc_ = bit_vector(string_of_changes_.size(), 0);

        for (size_t i = 0; i < string_of_changes_.size(); i++)
            if (string_of_changes_[i] == '#')
                loc_[i] = 1;

        for (auto i : start_possitions_)
            iloc_[i] = 1;


        // std::cout << new_original_ << std::endl;
        return 0;
    }

}