#include "Bio_FMi.h"

namespace bio_fmi
{

    Bio_FMi::Bio_FMi(unsigned context_length)
        : context_length_(context_length){
            original_text_change_ = false;
            select_support_mcl<> sloc_;
            rank_support_v<> rloc_;
            rank_support_v<> riloc_;
        }

    Bio_FMi::~Bio_FMi() {
    }

    int Bio_FMi::build()
    {
        try
        {
            std::cout << "Context length: " << context_length_ << std::endl;
            std::cout << "Number of segments: " << number_of_segments_ << std::endl;
            std::cout << "Number of changes: " << number_of_changes_ << std::endl;
            std::cout << "Total text size: " << total_text_size_ << std::endl;
            std::cout << std::endl;

            std::cout << "(1/2) Building fm-index over reference string";
            construct_im(reference_index_, reference_string_, 1);
            std::cout << " ... done" << std::endl;

            std::cout << "(2/2) Building fm-index over string of changes";
            construct_im(changes_index_, string_of_changes_, 1);
            std::cout << " ... done" << std::endl;

            is_empty_ = false;

            std::cout << "Size of reference index: " << reference_index_.size() << std::endl;
            std::cout << "Size of change index: " << changes_index_.size() << std::endl;
            std::cout << std::endl;

        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return -1;
        }

        return 0;
    }

    int Bio_FMi::read_patterns(std::filesystem::path patterns_input_file, std::vector<std::string> &patterns)
    {

        std::cout << "reading pattern file" << std::endl;

        /*  open file   */
        std::ifstream input_file(patterns_input_file);
        if (!input_file.is_open())
        {
            // parameters.errorExitCode = -1;  //  TODO jaka bude chyba
            return -1;
        }

        unsigned number_of_patterns;
        unsigned pattern_length;
        std::string info;
        std::string line;

        std::getline(input_file, info);
        std::getline(input_file, line);

        sscanf(info.c_str(), "# number=%u length=%u", &number_of_patterns, &pattern_length);

        // std::cout << number_of_patterns << ',' << pattern_length << std::endl;
        int N = line.size();
        int j = 0;
        std::string res = "";
        while (j < N)
        {
            res += line[j];
            if (res.size() == pattern_length)
            {
                patterns.push_back(res);
                res = "";
            }
            j++;
        }
        return 0;
    }

    int Bio_FMi::save(std::filesystem::path save_path)
    {
        try
        {
            std::filesystem::path filename = save_path.filename();
            save_path = save_path.replace_extension();
            std::cout << "Saving index to " << save_path << std::endl;
            std::filesystem::create_directory(save_path);

            /*  save fm-indexes */
            store_to_file(reference_index_, save_path / filename.replace_extension(".ri")); // save I0
            store_to_file(changes_index_, save_path / filename.replace_extension(".ci"));   //    save Id
            /*  save bit_vectors    */
            store_to_file(loc_, save_path / filename.replace_extension(".loc"));   //  save bitvector loc_
            store_to_file(iloc_, save_path / filename.replace_extension(".iloc")); //  save bitvector iloc_
            /*  save context length*/
            store_to_file(context_length_, save_path / filename.replace_extension(".lc")); //  save lc
            /**/
            store_to_file(change_lengths_, save_path / filename.replace_extension(".cl"));
            store_to_file(start_possitions_, save_path / filename.replace_extension(".sp"));
            store_to_file(base_position_, save_path / filename.replace_extension(".abp"));
            store_to_file(offset_, save_path / filename.replace_extension(".aof"));

            if (original_text_change_)
            {
                std::cout << "saving new original text to: " << filename.replace_extension(".reposition.eds") << std::endl;
                std::ofstream outfile(save_path / filename);
                if (!outfile.is_open())
                {
                    std::cout << "cannot open output file" << std::endl;
                }
                outfile.write(new_original_.c_str(), new_original_.size());
                outfile.close();
            }   
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return -1;
        }
        return 0;
    }

    int Bio_FMi::load(std::filesystem::path index_path)
    {
        try
        {
            std::filesystem::path filename = index_path.filename();
            index_path = index_path.replace_extension();
            std::cout << "Loading index from " << index_path << std::endl;
            std::filesystem::create_directory(index_path);

            /*  load fm-indexes */
            load_from_file(reference_index_, index_path / filename.replace_extension(".ri")); // load I0
            load_from_file(changes_index_, index_path / filename.replace_extension(".ci"));   //    load Id
            /*  load bit_vectors    */
            load_from_file(loc_, index_path / filename.replace_extension(".loc"));   //  load bitvector loc
            load_from_file(iloc_, index_path / filename.replace_extension(".iloc")); //  load bitvector iloc
            /*  load context length*/
            load_from_file(context_length_, index_path / filename.replace_extension(".lc")); //  load lc
            /*  load additional inform vectors  */
            load_from_file(start_possitions_, index_path / filename.replace_extension(".sp")); //  load sample_start vector
            load_from_file(change_lengths_, index_path / filename.replace_extension(".cl"));  //  load aOp vector
            load_from_file(base_position_, index_path / filename.replace_extension(".abp")); //  load aBasepos vector
            load_from_file(offset_, index_path / filename.replace_extension(".aof"));  //  load aOffset vector

            riloc_ = rank_support_v<>(&iloc_);
            rloc_ = rank_support_v<>(&loc_);
            sloc_ = select_support_mcl<>(&loc_);

            number_of_changes_ = change_lengths_.size();
            number_of_segments_ = start_possitions_.size();

        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return -1;
        }
        return 0;
    }

    void Bio_FMi::print()
    {

        std::cout << "reference string:     " << reference_string_ << std::endl;
        std::cout << "string of changes:    " << string_of_changes_ << std::endl;
        std::cout << "loc_:                  " << loc_ << std::endl;
        std::cout << "iloc_:                 " << iloc_ << std::endl;
        std::cout << "number of segments: " << number_of_segments_ << std::endl;
        std::cout << "number of changes: " << number_of_changes_ << std::endl;

        std::cout << "start_possitions_: ";
        for (auto &i : start_possitions_)
            std::cout << i << ", ";
        std::cout << std::endl;

        std::cout << "aBasePos: ";
        for (auto &i : base_position_)
            std::cout << i << ", ";
        std::cout << std::endl;

        std::cout << "aOffset: ";
        for (auto &i : offset_)
            std::cout << i << ", ";
        std::cout << std::endl;

        std::cout << "change_lengths_: ";
        for (auto &i : change_lengths_)
            std::cout << i << ", ";
        std::cout << std::endl;
    }

}
