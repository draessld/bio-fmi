#include "index.hpp"

double getFolderSize(const std::filesystem::path &folderPath)
{
    std::uintmax_t totalSize = 0;

    // Iterate through each file and subdirectory in the folder
    for (const auto &entry : std::filesystem::recursive_directory_iterator(folderPath))
    {
        if (std::filesystem::is_regular_file(entry.status()))
        { // Only count regular files
            totalSize += std::filesystem::file_size(entry.path());
        }
    }
    double totalSizeMB = static_cast<double>(totalSize) / (1024 * 1024); // Convert to MB
    return totalSizeMB;
}

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

        Bio_FMi::Bio_FMi(EDS eds, int context_length)
        : context_length_(context_length), eds(eds)
    {
        //  create folder for the index
        std::string name = "eds_index."+std::to_string(context_length);
        std::string ext = ".leds";
        index_bed_ = std::filesystem::current_path() / (name + ext + ".index");

        std::cout << "Index destination on " << index_bed_ << std::endl;
        if (!std::filesystem::exists(index_bed_))
            std::filesystem::create_directories(index_bed_);
        index_bed_.append(name);

        //  metadata files - for reference string and the string of changes
        reference_filepath_ = index_bed_;
        reference_filepath_.replace_extension(ext + ".metadata.ref");

        changes_filepath_ = index_bed_;
        changes_filepath_.replace_extension(ext + ".metadata.chan");
    }

    int Bio_FMi::parse_eds()
        {

            if (eds.empty())
            {
                std::ifstream ifs(eds_file_);
                if (!ifs.is_open())
                {
                    std::cerr << "Error: Unable to open eds file " << eds_file_ << std::endl;
                    return -3;
                }
                eds = EDS(ifs);
                ifs.close();
            }

            std::ofstream ref_file(reference_filepath_, std::ios::out);
            std::ofstream chan_file(changes_filepath_, std::ios::out);
            if (!ref_file.is_open() || !chan_file.is_open())
            {
                std::cerr << "Error: Unable to open metadata file " << eds_file_ << std::endl;
                return -3;
            }

            unsigned int cl = context_length_ -1;
            std::string context_r("");
            std::string context_l("");

            ref_file<<"#";
            chan_file<<"#";
            
            size_t chi = 0;
            bool ref = false;
            int basePos = 0;
            int setSize = 0;
            tloc_ = sdsl::bit_vector(eds.l_common + eds.n_common + 1, 0);
            loc_ = sdsl::bit_vector(eds.N + (eds.m * cl *2), 0);
            iloc_ = sdsl::bit_vector(eds.N + (eds.m * cl *2), 0);

            tloc_[0] = 1;
            loc_[0] = 1;
            iloc_[0] = 1;
            std::streampos pos;

            if (!eds.is_ref[0])
            {
                base_position_.push_back(0);
            }
            

            for (size_t x = 0; x < eds.n; x++){
                if (eds.is_ref[x])
                {
                    basePos += eds.changes[chi].size();
                    base_position_.push_back(basePos);
                    ref_file << eds.changes[chi];
                    tloc_[ref_file.tellp()] = 1;
                    ref_file << '#';
                    // std::cout << eds.changes[chi].size() << ',' << cl <<std::endl;
                    if(eds.changes[chi].size() < cl){
                        // std::cout << "if" << std::endl;
                        context_l = eds.changes[chi];
                    }else{
                        // std::cout << "else" << std::endl;
                        context_l = eds.changes[chi].substr(eds.changes[chi].size()-cl,cl);
                    }
                }else{
                    setSize += eds.set_size[x];
                    set_size_.push_back(setSize);
                    if ((chi + eds.set_size[x]) >= eds.m){
                        context_r = "";
                    }else{
                        // std::cout << eds.set_size[x]<<','<<eds.changes[chi + eds.set_size[x]].size()<< ',' << cl << std::endl;
                        if(cl > eds.changes[chi + eds.set_size[x]].size()){
                        // std::cout << "if2" << std::endl;
                            context_r = eds.changes[chi + eds.set_size[x]];
                        }else{
                        // std::cout << "else2" << std::endl;
                            context_r = eds.changes[chi + eds.set_size[x]].substr(0,cl);
                        }
                    }

                    for (size_t i = chi; i < (chi + eds.set_size[x]); i++)
                    {
                        offset_.push_back(eds.changes[i].size());
                        chan_file << context_l;
                        chan_file << eds.changes[i];
                        chan_file << context_r;
                        pos = chan_file.tellp(); 
                        loc_[pos] = 1;
                        chan_file << '#';
                    }
                    iloc_[pos] = 1;
                
                }

                chi += eds.set_size[x];
            }

            n = eds.n;
            m = eds.m;
            N = eds.N;

            // std::cout << T0<< std::endl;
            // std::cout << Td << std::endl;

            // std::ofstream ref_file(reference_filepath_, std::ios::out);
            // std::ofstream chan_file(changes_filepath_, std::ios::out);
            // if (!ref_file.is_open() || !chan_file.is_open())
            // {
            //     std::cerr << "Error: Unable to open metadata file " << eds_file_ << std::endl;
            //     return -3;
            // }

            // ref_file << T0;
            // chan_file << Td;

            ref_file.close();
            chan_file.close();

            return 0;
        }
/*
    int Bio_FMi::parse_eds()
    {

        std::ifstream ifs(eds_file_, std::ios::in | std::ios::binary);
        if (!ifs.is_open())
        {
            std::cerr << "Error: Unable to open eds file " << eds_file_ << std::endl;
            return -3;
        }

        EDS eds(ifs);
        ifs.clear();
        ifs.seekg(0, std::ios::beg);

        std::ofstream ref_file(reference_filepath_, std::ios::out);
        std::ofstream chan_file(changes_filepath_, std::ios::out);
        if (!ref_file.is_open() || !chan_file.is_open())
        {
            std::cerr << "Error: Unable to open metadata file " << eds_file_ << std::endl;
            return -3;
        }

        unsigned int cl = context_length_ - 1;
        std::string context_r("");
        std::string context_l("");

        std::string T0("#");
        std::string Td("#");
        ref_file << '#';
        chan_file << '#';

        unsigned int chi = 0;
        int basePos = 0;
        int setSize = 0;
        char c;
        // char* buffer;
        std::string buffer;

        for (size_t x = 0; x < eds.n; x++)
        {
            if (eds.is_ref[x])
            {
                base_position_.push_back(basePos);
                basePos += eds.changes[chi].size();

                buffer.resize(eds.lengths[chi], '\0');
                ifs.seekg(eds.base_position[chi]);
                ifs.read(&buffer[0], eds.lengths[chi]);
                ref_file << buffer;
                // T0+=eds.changes[chi];
                // T0+='#';
                ref_file << '#';

                //  TODO
                context_l = eds.changes[chi].substr(eds.changes[chi].size() - cl, cl);
            }
            else
            {
                set_size_.push_back(setSize);
                setSize += eds.set_size[x];
                if ((chi + eds.set_size[x]) >= eds.m)
                {
                //  TODO
                    context_r = "";
                }
                else
                {
                    //  TODO
                    // context_r = eds.changes[chi + eds.set_size[x]].substr(0, cl);


                    buffer.resize(cl, '\0');
                    ifs.seekg(eds.base_position[chi]);
                    ifs.read(&buffer[0], cl);
                }

                for (size_t i = chi; i < (chi + eds.set_size[x]); i++)
                {
                    offset_.push_back(eds.lengths[i]);

                    buffer.resize(eds.lengths[i], '\0');
                    ifs.seekg(eds.base_position[i]);
                    ifs.read(&buffer[0], eds.lengths[i]);

                    // Td += context_l;
                    // Td += eds.changes[i];
                    // Td += context_r;
                    // Td += '#';
                    chan_file << context_l;
                    chan_file << buffer;
                    chan_file << context_r;

                    chan_file << '#';
                }
            }

            chi += eds.set_size[x];
        }

        n = eds.n;
        m = eds.m;
        N = eds.N;

        // std::cout << T0<< std::endl;
        // std::cout << Td << std::endl;

        ifs.close();

        // ref_file << T0;
        // chan_file << Td;

        ref_file.close();
        chan_file.close();

        // std::ifstream file(changes_filepath_, std::ios::binary);
        // //  create inary vectors
        // file.seekg(0, file.end);
        // size_t pos = file.tellg();
        // file.seekg(0, file.beg);
        loc_ = sdsl::bit_vector(set_size_.size(), 0);
        iloc_ = sdsl::bit_vector(set_size_.size(), 0);
        // int i = 0;
        // char c;
        // chi = 0;
        // int si = 0;
        // while (file.get(c))
        // {
        //     if (c == '#'){
        //         chi++;
        //         loc_[i] = 1;
        //     }
        //     if (chi-1 == set_size_[si])
        //     {
        //         si++;
        //         iloc_[i] = 1;
        //     }
        //     i++;
        // }
        // iloc_[pos - 1] = 1;

        // file.close();

        // std::ifstream file_ref(reference_filepath_, std::ios::binary);
        // //  create inary vectors
        // file_ref.seekg(0, file_ref.end);
        // pos = file_ref.tellg();
        // file_ref.seekg(0, file_ref.beg);
        tloc_ = sdsl::bit_vector(n, 0);
        // i = 0;
        // while (file_ref.get(c))
        // {
        //     if (c == '#')
        //         tloc_[i] = 1;
        //     i++;
        // }

        // file_ref.close();
        return 0;
    }
    */

    Bio_FMi::Bio_FMi(std::filesystem::path index_folder)
    {
        //  find context length in folder name
        std::string filename = index_folder.filename().replace_extension("");
        index_bed_ = index_folder / filename;
        // eds_file_ = index_folder;
        // eds_file_.replace_extension("");

        std::cout << "Index will be loaded from " << index_bed_ << std::endl;

        // Extract the number from the filename

        // size_t i = filename.size() - 3;
        // for (; i > 0; i--)
        // {
        //     if (filename[i] == '.')
        //         break;
        // }
        // std::string number_str = filename.substr(i + 1, filename.size() - i - 5);
        // context_length_ = 0;
        // // Convert the extracted string to an integer
        // try
        // {
        //     context_length_ = std::stoi(number_str);
        //     std::cout << "Context length was found as : " << context_length_ << std::endl;

        //     // Print the extracted number
        // }
        // catch (const std::exception &e)
        // {
        //     std::cerr << "Error: context length was not found. Please setup explicitly: " << std::endl;
        //     std::cin >> context_length_;
        // }

        load();
        // print();
    }

    Bio_FMi::~Bio_FMi()
    {
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

    int Bio_FMi::locate(const std::string& P)
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

        //  check pattern length and set the 
        //  if the pattern size is not divisible by the context_length split tha last two chunk into evenly parts
        if ((P.size() % context_length_) != 0)
        {
            std::cout << "Unsupported pattern length = needs to be product of context_length" << std::endl;
            return -1;
        }

        //  for each chunk 
        for (chunk_index = 0; chunk_index < (P.size() / context_length_); chunk_index++)
        {
            chunk_start_position = chunk_index * context_length_;
            chunk = P.substr(chunk_start_position, context_length_);
            // std::cout << "Searching for a chunk number: " << chunk_index << " starting on position in pattern P:" << chunk_start_position << "=" << chunk << std::endl;

            //  SEARCH in reference 
            auto ref_locations = sdsl::locate(reference_index_, chunk);
            // std::cout << "in I0: " << ref_locations.size() << std::endl;

            // //  SEARCH in changes 
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
                // std::cout << "chunk found on position: " << loc << " "<< loc - context_length_<< ", block number: " << block_number-1 << std::endl;

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
                            
                            // std::cout << "set size: " << set_size_[block_number-1] << std::endl;
                            if (occ.second.empty())
                            {
                                //  case 1
                                // std::cout << "Case 1 found. Origin position " << occ.first << std::endl;
                                new_hash_map_[loc].push_back(occ);
                            }
                            else if (occ.second.back() <= set_size_[block_number-1])
                            {
                            // std::cout << occ.second.back() << std::endl;
                                //  case 3
                                // std::cout << "Case 3 found. Origin position " << occ.first << std::endl;
                                //  check if it across any change
                                if (occ.second.back() <= set_size_[block_number - 1])
                                {
                                    // it lies in more than previous change set => the last must contains empty change
                                    // auto cp = occ;
                                    // cp.second.push_back(set_size_[block_number - 1] + 1); //  we expect that changes are ordered either lexicographicaly or by size = empty element is always first
                                    new_hash_map_[loc].push_back(occ);
                                }
                                else
                                {
                                    new_hash_map_[loc].push_back(occ);
                                }
                            }else{
                                std::cout << occ.second.back() << std::endl;
                                std::cout << "Here " << occ.first << std::endl;
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

                // std::cout << loc << "hehe"<<block_number << " " << change_number << " " << pre_hash_loc << std::endl;

                // pos_hash_loc = sloc_(change_number + 1);
                offset = loc - (pre_hash_loc + context_length_ - 1);
                    previous_outside_change = ((pre_hash_loc) >= (loc - context_length_));

                //  there could be empty context when the eds starts with set of changes
                if (base_position_[block_number] < (context_length_ - 1))
                {
                    //  context length is not full
                    //  1+5 > 7-5+1
                    // std::cout << "here" << std::endl;
                    loc -= pre_hash_loc;
                }
                else
                {
                    // previous_outside_change = ((pre_hash_loc + context_length_) > (loc - context_length_ + 1));
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
                    // std::cout << "Checking position " << loc - context_length_ << " in hash table" << std::endl;

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
                                }else{
                                    //  case 3 
                                    // std::cout << "Case ???? found. Origin position " << occ.first << std::endl;
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

            if (new_hash_map_.empty())
            {
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
            // //  save bit_vectors    */
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
        // std::cout << "N: " << N << std::endl;
        // std::cout << "n: " << n << std::endl;
        // std::cout << "m: " << m << std::endl;
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
        std::cout << "n: " << n << std::endl;
        std::cout << "N: " << N << std::endl;
        std::cout << "m: " << m << std::endl;
        std::cout << "EDS file size: " << std::filesystem::file_size(eds_file_) << " B" << std::endl;
        std::cout << "Total index size: " << total_index_size_ << " MB" << std::endl;
        std::cout << "index folder size: " << getFolderSize(index_bed_.parent_path()) << " MB" << std::endl;
        std::cout << std::endl;
    }

    Bio_FMi::hash_type Bio_FMi::get_result()
    {
        return old_hash_map_;
    }
}
