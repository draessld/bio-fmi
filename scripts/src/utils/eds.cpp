#include "eds.h"

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

    char ch;
    while (is.get(ch))
    {
        original_eds.push_back(ch);
        switch (ch)
        {
        case '{':
            local_m =1;
            n++;
            break;
        case '}':
            if (local_m == 1){
                n_common++;
                l_common+=line.size();
                //  Common part
                if (line.size() > max_l)
                    max_l = line.size();
                if (line.size() < min_l)
                    min_l = line.size();
                // std::cout << line.size() << ","<< int(min_l) << ","<< int(max_l) << std::endl;
            }
            m+=local_m; //  add number of string in current set
            // local_m = 0;
            N+=line.size();
            line.clear();
            break;
        case ',':
            if (line.empty()){
                //  an empty string
                n_empty_strings++;
            }
            local_m++;
            N+=line.size();
            line.clear();
            break;
        default:
            line.push_back(ch);
            break;
        }
    }
    avg_l = l_common / n_common;
}

EDS::~EDS() {
    // Destructor implementation (if needed)
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

    return 0;
}