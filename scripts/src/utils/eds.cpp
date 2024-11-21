#include "eds.hpp"

EDS::EDS(std::istream &is)
{

    unsigned int local_m = 0; //  number string in a set
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
    int x = 0;
    n_empty_strings = 0;
    size_t i = 0; // absolute iterator over file

    char ch;
    while (is.get(ch))
    {
        switch (ch)
        {
        case '{':
            if (!line.empty())
            {
                //  reference before
                n_common++;
                l_common += line.size();
                changes.push_back(line);
                lengths.push_back(line.size());
                base_position.push_back(i - line.size());
                abs_change_number.push_back(-1);
                ref_position.push_back(l_common);
                set_size.push_back(1);
                if (n > 0)
                    cum_set_size.push_back(cum_set_size.back() + 1);
                else
                    cum_set_size.push_back(1);

                is_ref.push_back(1);
                n++;
                m++;

                if (line.size() > max_l)
                    max_l = line.size();
                if (line.size() < min_l)
                    min_l = line.size();
            }
            N += line.size();
            local_m = 1;
            n++;
            line.clear();
            break;
        case '}':
            changes.push_back(line);
            lengths.push_back(line.size());
            base_position.push_back(i - line.size());
            ref_position.push_back(l_common);
            if (line.empty())
            {
                //  an empty string
                n_empty_strings++;
            }
            if (local_m == 1)
            {
                n_common++;
                l_common += line.size();
                abs_change_number.push_back(-1);
                //  Common part
                ref_position.push_back(l_common);
                if (line.size() > max_l)
                    max_l = line.size();
                if (line.size() < min_l)
                    min_l = line.size();
                is_ref.push_back(1);
            }
            else
            {
                abs_change_number.push_back(x++);
                total_change_size += line.size();
                is_ref.push_back(0);
            }

            set_size.push_back(local_m);
            if (n != 1)
                cum_set_size.push_back(cum_set_size.back() + local_m);
            else
                cum_set_size.push_back(local_m);
            m += local_m; //  add number of string in current set
            N += line.size();
            line.clear();
            break;
        case ',':
            changes.push_back(line);
            lengths.push_back(line.size());
            abs_change_number.push_back(x++);
            total_change_size += line.size();
            base_position.push_back(i - line.size());
            ref_position.push_back(l_common);
            if (line.empty())
            {
                //  an empty string
                n_empty_strings++;
            }
            local_m++;
            N += line.size();
            line.clear();
            break;
        case '\n':
            break;
        default:
            line.push_back(ch);
            break;
        }
        i++;
    }
    if (!line.empty())
    {
        //  reference in the end
        n_common++;
        l_common += line.size();
        changes.push_back(line);
        lengths.push_back(line.size());
        ref_position.push_back(l_common);
        base_position.push_back(i - line.size());
        abs_change_number.push_back(-1);
        set_size.push_back(1);
        if (n != 0)
            cum_set_size.push_back(cum_set_size.back() + 1);
        else
            cum_set_size.push_back(1);
        is_ref.push_back(1);
        n++;
        N += line.size();
        m++;

        if (line.size() > max_l)
            max_l = line.size();
        if (line.size() < min_l)
            min_l = line.size();

        line.clear();
    }

    if (n_common == 0)
    {
        avg_l = 0;
    }
    else
    {
        avg_l = l_common / n_common;
    }
    is_empty = false;
}

int EDS::stats()
{
    std::cout << "total_change_size:" << total_change_size << std::endl;
    std::cout << "n_common:" << n_common << std::endl;
    std::cout << "l_common:" << l_common << std::endl;
    std::cout << "n:" << n << std::endl;
    std::cout << "N:" << N << std::endl;
    std::cout << "m:" << m << std::endl;
    std::cout << "min_l:" << min_l << std::endl;
    std::cout << "max_l:" << max_l << std::endl;
    std::cout << "avg_l:" << avg_l << std::endl;
    std::cout << "empty string:" << n_empty_strings << std::endl;
    std::cout << "Total class size:" << calculateSize() << " MB" << std::endl;
    std::cout << "changes:" << changes.size() << std::endl;
    for (auto change : changes)
    {
        std::cout << change << ", ";
    }
    std::cout << std::endl;
    std::cout << "sets:" << set_size.size() << std::endl;

    std::cout << "cum_set_size:" << cum_set_size.size() << std::endl;
    for (auto set : cum_set_size)
    {
        std::cout << set << ", ";
    }
    std::cout << std::endl;

    std::cout << "set_size:" << set_size.size() << std::endl;
    for (auto set : set_size)
    {
        std::cout << set << ", ";
    }
    std::cout << std::endl;

    std::cout << "is_ref:" << is_ref.size() << std::endl;
    for (auto b : is_ref)
    {
        std::cout << b << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "base_position:" << base_position.size() << std::endl;
    for (auto bp : base_position)
    {
        std::cout << bp << ", ";
    }
    std::cout << std::endl;

    std::cout << "ref_position:" << ref_position.size() << std::endl;
    for (auto bp : ref_position)
    {
        std::cout << bp << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "lengths:" << lengths.size() << std::endl;
    for (auto l : lengths)
    {
        std::cout << l << ", ";
    }
    std::cout << std::endl;

    std::cout << "abs_change_number:" << abs_change_number.size() << std::endl;
    for (auto l : abs_change_number)
    {
        std::cout << l << ", ";
    }
    std::cout << std::endl;
    return 0;
}

template <typename T>
size_t calculateContainerSize(const T &container)
{
    size_t size = sizeof(container); // Base container size
    size += container.capacity() * sizeof(typename T::value_type);
    if constexpr (std::is_same<typename T::value_type, std::string>::value)
    {
        for (const auto &str : container)
        {
            size += str.capacity();
        }
    }
    return size;
}

double EDS::calculateSize()
{
    size_t totalSize = sizeof(this);

    // Include dynamic allocations
    totalSize += calculateContainerSize(changes);
    totalSize += calculateContainerSize(set_size);
    totalSize += calculateContainerSize(cum_set_size);
    totalSize += calculateContainerSize(is_ref);
    totalSize += calculateContainerSize(base_position);
    totalSize += calculateContainerSize(lengths);

    return static_cast<double>(totalSize) / (1024 * 1024);
}

int get_random_position(int n)
{
    if (n == 0)
    {
        return 0;
    }

    std::random_device rd;                         // Obtain a random number from hardware (if available)
    std::mt19937 gen(rd());                        // Initialize the generator with the random seed
    std::uniform_int_distribution<> dis(0, n - 1); // Define the range [0, n-1]

    return dis(gen); // Generate and return a random number in the range
}

int EDS::gen_pattern(std::ostream &os, unsigned int size)
{
    unsigned int rset;
    unsigned int rstr;
    unsigned int rpos;
    
    unsigned int spos;
    std::vector<int> uch;
    
    int attemps = 50;
    std::string res;
    std::string tmp;
    while (attemps)
    {
        res = "";
        uch.clear();
        rset = get_random_position(n); // get random starging set
        // std::cout << rset << ", "<< set_size[rset]<<std::endl;
        rstr = get_random_position(set_size[rset]); //  get random string in set
        // std::cout << rstr <<", " << changes[cum_set_size[rset] - rstr].size()<< std::endl;
        rpos = get_random_position(changes[cum_set_size[rset] -1- rstr].size()); // get random position in string
        // std::cout << rpos << std::endl;
        // std::cout << "startign position:" << ref_position[cum_set_size[rset] - rstr] + rpos << std::endl;
        spos = ref_position[cum_set_size[rset] -1 - rstr] + rpos;

        while (rset < n - 1)
        {
            if (abs_change_number[cum_set_size[rset]-1 - rstr] !=-1){
            if(uch.empty()){
                uch.push_back(abs_change_number[cum_set_size[rset] -1- rstr]);
            }else if(uch.back() != abs_change_number[cum_set_size[rset]-1 - rstr])
                uch.push_back(abs_change_number[cum_set_size[rset]-1 - rstr]);
            }
            
            tmp = changes[cum_set_size[rset] -1- rstr];
            // std::cout << tmp << ',' << rpos << ',' << (res.size()+tmp.size()-rpos) << std::endl;
            if ((res.size() + tmp.size() - rpos) < size)
            {
                // std::cout << tmp.substr(rpos, tmp.size() - rpos)<< std::endl;
                res.append(tmp.substr(rpos, tmp.size() - rpos));
                rset++; // move to next set
                if (rset >= n)
                {
                    res = "";
                    uch.clear();
                    continue;
                }
                rstr = get_random_position(set_size[rset]); //  get random string in set
                rpos = 0;
                // std::cout << "new set:" << rset << is_ref[rset]<< std::endl;
                // std::cout << "new position:" << rstr << std::endl;
            }
            else
            {
                // std::cout << tmp.substr(rpos, size - res.size())<< std::endl;
                res.append(tmp.substr(rpos, size - res.size()));
                // std::cout << "done" << std::endl;
                std::cout << spos+1 << "[";
                    for (auto ch: uch)
                    {
                        std::cout << ch+1 << ' ';
                    }
                std::cout << "]\n";
                os << res << '\n';
                return 0;
            }
        }
        // std::cout <<"reset" << std::endl; 
        attemps--;
    }
    std::cout << "unsucessful" << std::endl;
    return 1;
}

void EDS::gen_patterns(std::ostream &os, unsigned int t, unsigned int size)
{
    // os << "#t:" << t << "\t#l:" << size << '\n';
    for (size_t i = 0; i < t; i++)
    {
        gen_pattern(os, size);
    }
}


void EDS::extract(unsigned int position, std::vector<unsigned int> changes){
    unsigned int upp,low;
    upp = ref_position.size();
    low = 0;
    unsigned int pivot;
    while (upp>low)
    {
        pivot = (upp-low)/2;
        if (ref_position[pivot] < position)
        {
            low = pivot;
        }else{
            upp = pivot;
        }
    }
    std::cout << "Relative position: " << ref_position[low] << std::endl;
    position = base_position[low] + (position - ref_position[low]);
    std::cout << "File position: " << ref_position[low] << std::endl;
}


void EDS::language(std::ostream &os)
{
    std::cout << "Warning! This will take very much space" << std::endl;
    for (size_t i = 0; i < n; i++)
    {
        // for every set size
    }
}