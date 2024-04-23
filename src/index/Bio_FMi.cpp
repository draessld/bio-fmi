#include "Bio_FMi.h"

namespace bio_fmi
{

    Bio_FMi::Bio_FMi(std::filesystem::path eds_file, int context_length)
        : context_length_(context_length), eds_file_(eds_file)
    {
        //  create folder for the index
        std::string filename = eds_file.filename();
        std::string ext = eds_file.extension().c_str();
        index_bed_ = eds_file;
        index_bed_ = index_bed_.replace_extension(ext + ".index");
        std::cout << "Index destination on " << index_bed_ << std::endl;
        if (!std::filesystem::exists(index_bed_))
            std::filesystem::create_directories(index_bed_);
        index_bed_.append(filename);

        //  metadata files - for reference string and the string of changes
        reference_filepath_ = index_bed_;
        reference_filepath_.replace_extension(ext + ".metadata.ref");

        changes_filepath_ = index_bed_;
        changes_filepath_.replace_extension(ext + ".metadata.chan");
    }

    Bio_FMi::Bio_FMi(std::filesystem::path index_folder)
    {
        //  find context length in folder name
        std::string filename = index_folder.filename().replace_extension("");
        index_bed_ = index_folder / filename;
        // eds_file_ = index_folder;
        // eds_file_.replace_extension("");

        std::cout << "Index will be loaded from " << index_bed_ << std::endl;

        // Extract the number from the filename

        size_t i = filename.size() - 5;
        for (; i > 0; i--)
        {
            if (filename[i] == '.')
                break;
        }
        std::string number_str = filename.substr(i + 1, filename.size() - i - 5);
        context_length_ = 0;
        // Convert the extracted string to an integer
        try
        {
            context_length_ = std::stoi(number_str);
            std::cout << "Context length was found as : " << context_length_ << std::endl;

            // Print the extracted number
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error: context length was not found. Please setup explicitly: " << std::endl;
            std::cin >> context_length_;
        }

        load();
        // print();
    }

    Bio_FMi::~Bio_FMi()
    {
    }

    int Bio_FMi::parse_eds()
    {
        //  open file
        std::ifstream in(eds_file_);
        if (!in.is_open())
        {
            std::cerr << "Error: Unable to open file " << eds_file_ << std::endl;
            return -1;
        }

        std::ofstream ref_file(reference_filepath_, std::ios::out);
        std::ofstream chan_file(changes_filepath_, std::ios::out);
        if (!ref_file.is_open() || !chan_file.is_open())
        {
            std::cerr << "Error: Unable to open metadata file " << eds_file_ << std::endl;
            return -3;
        }

        //  Data are validated for context length long enough
        std::string r_context, l_context; // keep context of length -1 (to search pattern of length 5 keep 4 characters of contexts)
        size_t cl_lessone = context_length_ - 1;

        set_size_.push_back(0);
        std::vector<std::string> changes;
        std::vector<int> start_positions_;
        size_t nchange = 0;

        bool context = 0; //  which context save 0=left,1=right
        bool change_open = false;

        chan_file << '#';

        size_t i = 0;
        char c;
        while (in.get(c))
        {
            // std::cout << c << std::endl;
            switch (c)
            {
            case '{':
                base_position_.push_back(i);
                n++;
                //  flush ref
                ref_file << l_context;

                if (!changes.empty())
                {
                    set_size_.push_back(set_size_.back());
                    set_size_.back() += changes.size();
                    ref_file << '#';
                }

                //  flush changes
                start_positions_.push_back((size_t)chan_file.tellp() - 1);
                for (size_t k = 0; k < changes.size(); k++)
                {
                    if (l_context.size() > cl_lessone)
                        chan_file << l_context.substr(l_context.size() - cl_lessone, cl_lessone);
                    else
                        chan_file << l_context;

                    chan_file << changes[k];
                    chan_file << r_context.substr(0, cl_lessone);
                    chan_file << '#';
                    offset_.push_back(changes[k].size());
                }
                l_context = r_context;
                r_context.clear();
                changes.clear();
                changes.push_back("");
                ;

                //  setup open degenerate symbol
                nchange = 0;
                change_open = true;
                break;
            case '}':
                context = context ^ 1;
                change_open = false;
                break;
            case ',':
                total_deg_strings++;
                changes.push_back("");
                nchange++;
                break;
            default:
                //  its character
                if (change_open)
                {
                    //  character belongs to string in degenerate set
                    changes[nchange].push_back(c);
                }
                else
                {
                    //  character belongs to common part
                    r_context.push_back(c);
                    i++;
                }
                break;
            }
        }

        N = i;
        //  flush ref
        ref_file << l_context;
        if (!changes.empty())
        {
            set_size_.push_back(set_size_.back());
            set_size_.back() += changes.size();
            ref_file << '#';
        }
        ref_file << r_context;

        //  flush changes
        start_positions_.push_back((size_t)chan_file.tellp() - 1);
        for (size_t k = 0; k < changes.size(); k++)
        {
            if (l_context.size() >= cl_lessone)
                chan_file << l_context.substr(l_context.size() - cl_lessone, cl_lessone);
            else
                chan_file << l_context;

            chan_file << changes[k];
            if (r_context.size() > cl_lessone)
                chan_file << r_context.substr(r_context.size() - context_length_, cl_lessone);
            else
                chan_file << r_context;
            chan_file << '#';
            offset_.push_back(changes[k].size());
        }

        in.close();
        ref_file.close();
        chan_file.close();

        n = (n * 2) + 1;

        // base_position_.erase(base_position_.begin());

        std::ifstream file(changes_filepath_, std::ios::binary);
        //  create inary vectors
        file.seekg(0, file.end);
        size_t pos = file.tellg();
        file.seekg(0, file.beg);
        loc_ = sdsl::bit_vector(pos, 0);
        iloc_ = sdsl::bit_vector(pos, 0);
        i = 0;
        while (file.get(c))
        {
            if (c == '#')
                loc_[i] = 1;
            i++;
        }

        for (i = 0; i < start_positions_.size(); i++)
        {
            iloc_[start_positions_[i]] = 1;
        }
        iloc_[pos - 1] = 1;

        file.close();

        std::ifstream file_ref(reference_filepath_, std::ios::binary);
        //  create inary vectors
        file_ref.seekg(0, file_ref.end);
        pos = file_ref.tellg();
        file_ref.seekg(0, file_ref.beg);
        tloc_ = sdsl::bit_vector(pos, 0);
        i = 0;
        while (file_ref.get(c))
        {
            if (c == '#')
                tloc_[i] = 1;
            i++;
        }

        file_ref.close();

        return 0;
    }

    int Bio_FMi::build()
    {
        try
        {
            std::cout << "  (0/3) Parsing EDS";
            if (Bio_FMi::parse_eds())
            {
                std::cout << "Error: Uncomplete EDS parsing" << std::endl;
                return -1;
            }
            std::cout << " ... done" << std::endl;

            std::cout << "  (1/3) Building fm-index over reference string";
            construct(reference_index_, reference_filepath_, 1);
            std::cout << " ... done" << std::endl;

            std::cout << "  (2/3) Building fm-index over string of changes";
            construct(changes_index_, changes_filepath_, 1);
            std::cout << " ... done" << std::endl;

            std::cout << "  (3/3) Building Rank and Select supports";
            riloc_ = rank_support_v<>(&iloc_);
            rloc_ = rank_support_v<>(&loc_);
            rtloc_ = rank_support_v<>(&tloc_);
            sloc_ = select_support_mcl<>(&loc_);
            std::cout << " ... done" << std::endl;

            total_index_size_ = size_in_mega_bytes(reference_index_) + size_in_mega_bytes(changes_index_) + size_in_mega_bytes(iloc_) + size_in_mega_bytes(loc_) + size_in_mega_bytes(base_position_) + size_in_mega_bytes(offset_) + size_in_mega_bytes(riloc_) + size_in_mega_bytes(rloc_) + size_in_mega_bytes(rtloc_) + size_in_mega_bytes(sloc_);

            // std::filesystem::remove(reference_filepath_);
            // std::filesystem::remove(changes_filepath_);

            save();
            return 0;
        }
        catch (const std::exception &e)
        {
            std::cerr << '\n'
                      << e.what() << '\n';
            return -2;
        }

        return 0;
    }

    void Bio_FMi::print_hash(const hash_type &map)
    {
        for (const auto &entry : map)
        {
            std::cout << "Key: (" << entry.first << ")\n";
            for (const auto &vec_entry : entry.second)
            {
                std::cout << "  Value Pair: {" << vec_entry.first << ", {";
                for (const auto &sub_vec_entry : vec_entry.second)
                {
                    std::cout << sub_vec_entry << " ";
                }
                std::cout << "}}\n";
            }
            std::cout << std::endl;
        }
    }

    void Bio_FMi::print_result(const hash_type &hash_map)
    {
        for (const auto &pair : hash_map)
        {
            for (auto occ : pair.second)
            {
                std::cout << occ.first << "[ ";
                for (int num : occ.second)
                {
                    std::cout << num << " ";
                }
                std::cout << "]" << std::endl;
            }
        }
    }

    int Bio_FMi::locate(std::string P)
    {
        new_hash_map_.clear();
        old_hash_map_.clear();

        size_t chunk_index;
        size_t chunk_start_position;
        std::string chunk;
        int block_number = 0;
        int change_number = 0;
        int offset = 0;
        // size_t tmp;
        // int max_n = offset_.size();
        // size_t position;
        int pre_hash_loc;
        // int pos_hash_loc;
        // auto it = new_hash_map_.find({0,0});
        auto it = new_hash_map_.find(0);
        // bool next_outside_change = false;
        bool previous_outside_change = false;
        // bool next_in_change = false;

        /*  check pattern length and set the */
        //  if the pattern size is not divisible by the context_length split tha last two chunk into evenly parts
        if ((P.size() % context_length_) != 0)
        {
            std::cout << "Unsupported pattern length = needs to be product of context_length" << std::endl;
            return -1;
        }

        /*  for each chunk */
        for (chunk_index = 0; chunk_index < (P.size() / context_length_); chunk_index++)
        {
            chunk_start_position = chunk_index * context_length_;
            chunk = P.substr(chunk_start_position, context_length_);
            // std::cout << "Searching for a chunk number: " << chunk_index << " starting on position in pattern P:" << chunk_start_position << "=" << chunk << std::endl;

            /*  SEARCH in reference */
            auto ref_locations = sdsl::locate(reference_index_, chunk);
            // std::cout << "in I0: " << ref_locations.size() << std::endl;

            // /*  SEARCH in changes */
            auto change_locations = sdsl::locate(changes_index_, chunk);
            // std::cout << "in Id: " << change_locations.size() << std::endl;

            for (auto loc : ref_locations)
            {
                //  validate with saved locations
                block_number = rtloc_(loc);

                //  is the next position in the change? == is the rank number on loc+context_length increased by 1?
                // next_in_change = (block_number != riloc_(loc + context_length_));
                // std::cout << "next position in change: " << next_in_change << std::endl;

                loc = loc - block_number + 1;
                // std::cout << "chunk found on position: " << loc << ", block number: " << block_number << std::endl;

                if (chunk_index == 0) //  first positions
                {
                    new_hash_map_[loc] = {{loc, {}}};
                }
                else
                {
                    //  VALIDATE

                    it = old_hash_map_.find(loc - context_length_);
                    if (it != old_hash_map_.end())
                    {
                        for (auto occ : it->second)
                        {
                            if (occ.second.empty())
                            {
                                //  case 1
                                // std::cout << "Case 1 found. Origin position " << occ.first << std::endl;
                                new_hash_map_[loc].push_back(occ);
                            }
                            else if (occ.second.back() <= set_size_[block_number])
                            {
                                //  case 3
                                // std::cout << "Case 3 found. Origin position " << occ.first << std::endl;
                                //  check if it across any change
                                if (occ.second.back() <= set_size_[block_number - 1])
                                {
                                    // it lies in more than previous change set => the last must contains empty change
                                    auto cp = occ;
                                    cp.second.push_back(set_size_[block_number - 1] + 1); //  we expect that changes are ordered either lexicographicaly or by size = empty element is always first
                                    new_hash_map_[loc].push_back(cp);
                                }
                                else
                                {
                                    new_hash_map_[loc].push_back(occ);
                                }
                            }
                        }
                    }
                }
            }

            for (auto loc : change_locations)
            {

                //  validate with saved locations
                //  it is location in change ( could be case 3 or 2)
                block_number = riloc_(loc) - 1;
                change_number = rloc_(loc);
                pre_hash_loc = sloc_(change_number);
                // pos_hash_loc = sloc_(change_number + 1);
                offset = loc - (pre_hash_loc + context_length_ - 1);


                //  there could be empty context when the eds starts with set of changes
                if (base_position_[block_number] < (context_length_ - 1))
                {
                    //  context length is not full
                    //  1+5 > 7-5+1
                    previous_outside_change = ((pre_hash_loc) >= (loc - context_length_));
                    loc -= pre_hash_loc;
                }
                else
                {
                    // previous_outside_change = ((pre_hash_loc + context_length_) > (loc - context_length_ + 1));
                    previous_outside_change = ((pre_hash_loc) >= (loc - context_length_));
                    loc = base_position_[block_number] + offset;
                }

                //  will be the next position outside of change? Is actual change long enough to catch another chunk?
                // std::cout << "next position outside of change:" << next_outside_change << std::endl;
                // std::cout << "previous position outside of change:" << previous_outside_change << std::endl;

                offset = offset_[change_number - 1];
                // if (offset > loc)
                // {
                //     offset = 0;
                // }
                // std::cout << "chunk found on position: " << loc << " while tracking change number " << change_number << std::endl;

                if (chunk_index == 0)
                {

                    int tmp = loc - offset;
                    // std::cout << "initial saving " << tmp << std::endl;
                        if (new_hash_map_.find(loc - offset) == new_hash_map_.end())
                            new_hash_map_[loc - offset] = {};
                        new_hash_map_[loc - offset].push_back({loc, {change_number}});
                }
                else
                {
                    //  VALIDATE
                    // std::cout << "Checking position" << loc - context_length_ << "in hash table" << std::endl;

                    if (!previous_outside_change)
                    {
                        // std::cout << "in same change" << std::endl;
                        it = old_hash_map_.find(loc - offset - context_length_);
                        if (it != old_hash_map_.end())
                        {
                            for (auto occ : it->second)
                            {
                                //  validate change numbers
                                if (occ.second.empty())
                                {
                                    //  case 3 in the beginning              ---[---]
                                    // std::cout << "Case 3 found. Origin position " << occ.first << std::endl;
                                    // if (new_hash_map_.find(loc - offset) == new_hash_map_.end())
                                    //     new_hash_map_[loc - offset] = {};

                                    // auto cp = occ;
                                    // cp.second.push_back(change_number);
                                    // new_hash_map_[loc - offset].push_back(cp);
                                }
                                else if (occ.second.back() == change_number)
                                {
                                    //  case 2
                                    // std::cout << "Case 2 found. Origin position: " << occ.first << std::endl;
                                    if (new_hash_map_.find(loc - offset) == new_hash_map_.end())
                                        new_hash_map_[loc - offset] = {};

                                    new_hash_map_[loc - offset].push_back(occ);
                                }
                            }
                        }
                    }
                    else
                    {
                        // std::cout << "in previus change" << std::endl;
                        it = old_hash_map_.find(loc - context_length_);
                        if (it != old_hash_map_.end())
                        {
                            for (auto occ : it->second)
                            {

                                if (occ.second.empty())
                                {
                                    //  case 3 in the beginning              ---[---]
                                    // std::cout << "Case 3 found. Origin position " << occ.first << std::endl;
                                    if (new_hash_map_.find(loc - offset) == new_hash_map_.end())
                                        new_hash_map_[loc - offset] = {};

                                    auto cp = occ;
                                    cp.second.push_back(change_number);
                                    new_hash_map_[loc - offset].push_back(cp);
                                }
                                else if (occ.second.back() <= set_size_[block_number])
                                {
                                    //  case 4
                                    // std::cout << "Case 4 found. Origin position: " << occ.first << std::endl;
                                    if (new_hash_map_.find(loc - offset) == new_hash_map_.end())
                                        new_hash_map_[loc - offset] = {};

                                    auto cp = occ;
                                    cp.second.push_back(change_number);
                                    new_hash_map_[loc - offset].push_back(cp);
                                }
                            }
                        }
                    }
                }
            }

            // std::cout << std::endl;

            if (new_hash_map_.empty()){
                old_hash_map_.clear();
                return -1;
            }

            // print_hash(new_hash_map_);

            std::swap(old_hash_map_, new_hash_map_);
            new_hash_map_.clear();
        }

        // print_hash(old_hash_map_);

        return 0;
    }

    int Bio_FMi::save()
    {
        try
        {
            store_to_file(reference_index_, index_bed_.replace_extension(".ri")); // save I0
            store_to_file(changes_index_, index_bed_.replace_extension(".ci"));   //    save Id
            // /*  save bit_vectors    */
            store_to_file(loc_, index_bed_.replace_extension(".loc"));   //  save bitvector loc
            store_to_file(iloc_, index_bed_.replace_extension(".iloc")); //  save bitvector iloc
            store_to_file(tloc_, index_bed_.replace_extension(".tloc")); //  save bitvector iloc
            store_to_file(set_size_, index_bed_.replace_extension(".ss"));
            store_to_file(base_position_, index_bed_.replace_extension(".abp"));
            store_to_file(offset_, index_bed_.replace_extension(".aof"));
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return -1;
        }
        return 0;
    }

    int Bio_FMi::load()
    {
        try
        {
            load_from_file(reference_index_, index_bed_.replace_extension(".ri")); // load I0
            load_from_file(changes_index_, index_bed_.replace_extension(".ci"));   //    load Id
            load_from_file(loc_, index_bed_.replace_extension(".loc"));            //  load bitvector loc
            load_from_file(iloc_, index_bed_.replace_extension(".iloc"));          //  load bitvector iloc
            load_from_file(tloc_, index_bed_.replace_extension(".tloc"));          //  load bitvector iloc
            load_from_file(base_position_, index_bed_.replace_extension(".abp"));  //  load aBasepos vector
            load_from_file(offset_, index_bed_.replace_extension(".aof"));         //  load aOffset vector
            load_from_file(set_size_, index_bed_.replace_extension(".ss"));

            riloc_ = rank_support_v<>(&iloc_);
            rloc_ = rank_support_v<>(&loc_);
            rtloc_ = rank_support_v<>(&tloc_);
            sloc_ = select_support_mcl<>(&loc_);
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
        std::cout << "loc:                  " << loc_ << std::endl;
        std::cout << "iloc:                 " << iloc_ << std::endl;
        std::cout << "tloc:                 " << tloc_ << std::endl;
        std::cout << "Context length: " << context_length_ << std::endl;
        std::cout << "Number of segments: " << n << std::endl;
        std::cout << "Number of changes: " << total_deg_strings << std::endl;
        std::cout << std::endl;

        std::cout << "aBasePos: ";
        for (auto &i : base_position_)
            std::cout << i << ", ";
        std::cout << std::endl;

        std::cout << "set size: ";
        for (auto &i : set_size_)
            std::cout << i << ", ";
        std::cout << std::endl;

        std::cout << "aOffset: ";
        for (auto &i : offset_)
            std::cout << i << ", ";
        std::cout << std::endl;
    }

    void Bio_FMi::print_stats()
    {
        std::cout << "Context length: " << context_length_ << std::endl;
        // std::cout << "n: " << n << std::endl;
        // std::cout << "N: " << N << std::endl;
        // std::cout << "# denegenerated sets: " << total_deg_sets << std::endl;
        // std::cout << "# strings in deg-sets: " << total_deg_strings << std::endl;
        std::cout << "Total index size: " << total_index_size_ << std::endl;
        std::cout << std::endl;
    }

    Bio_FMi::hash_type Bio_FMi::get_result()
    {
        return old_hash_map_;
    }
}
