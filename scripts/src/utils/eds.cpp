#include "eds.hpp"

EDS::EDS(std::istream &is){

    uint8_t local_m = 0;    //  number string in a set
    std::string line;

    //  get size of data
    std::streampos current_pos = is.tellg();
    is.seekg(0, std::ios::end);
    original_input_size = is.tellg();
    is.seekg(current_pos, std::ios::beg);

    min_l = original_input_size;
    max_l = 0;
    n = 0;
    N = 0;
    m = 0;
    n_common = 0;
    l_common = 0;
    n_empty_strings = 0;
    size_t ri = 0;
    
    char ch;
    while (is.get(ch))
    {
        // original_eds.push_back(ch);
        // changes[chi].push_back(ch);
        switch (ch)
        {
        case '{':
            if (!line.empty())
            {
                //  reference before
                changes.push_back(line);
                lengths.push_back(line.size());
                set_size_.push_back(1);
                n++;
                m++;

                if (line.size() > max_l)
                    max_l = line.size();
                if (line.size() < min_l)
                    min_l = line.size();
            }
            N+=line.size();
            local_m =1;
            n++;
            line.clear();
            break;
        case '}':
            changes.push_back(line);
            if (line.empty()){
                //  an empty string
                n_empty_strings++;
            }
            if (local_m == 1){
                n_common++;
                l_common+=line.size();
                //  Common part
                if (line.size() > max_l)
                    max_l = line.size();
                if (line.size() < min_l)
                    min_l = line.size();
            }
            set_size_.push_back(local_m);
            m+=local_m; //  add number of string in current set
            N+=line.size();
            line.clear();
            break;
        case ',':
            changes.push_back(line);
            lengths.push_back(line.size());
            if (line.empty()){
                //  an empty string
                n_empty_strings++;
            }
            local_m++;
            N+=line.size();
            line.clear();
            break;
        case '\n':
            break;
        default:
            line.push_back(ch);
            break;
        }
    }
    if (!line.empty()){
        //  reference in the end
        changes.push_back(line);
        lengths.push_back(line.size());
        set_size_.push_back(1);
        n++;
        N+=line.size();
        m++;

        if (line.size() > max_l)
            max_l = line.size();
        if (line.size() < min_l)
            min_l = line.size();

        line.clear();
    }

    if (n_common==0)
    {
        avg_l = 0;
    }else{
        avg_l = l_common / n_common;
    }
    
}

int EDS::stats(){
    std::cout << "n_common:" << int(n_common) << std::endl;
    std::cout << "l_common:" << int(l_common) << std::endl;
    std::cout << "n:" << int(n) << std::endl;
    std::cout << "N:" << int(N) << std::endl;
    std::cout << "m:" << int(m) << std::endl;
    std::cout << "min_l:" << int(min_l) << std::endl;
    std::cout << "max_l:" << int(max_l) << std::endl;
    std::cout << "avg_l:" << int(avg_l) << std::endl;
    std::cout << "empty string:" << int(n_empty_strings) << std::endl;
    std::cout << "changes:" << changes.size() << std::endl;
        for (auto change : changes)
    {
        std::cout << change << ", ";
    }
    std::cout <<std::endl;
    std::cout << "sets:" << set_size_.size() << std::endl;

        for (auto set : set_size_)
    {
        std::cout << set << ", ";
    }
    std::cout <<std::endl;
    return 0;
}