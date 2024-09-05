#ifndef UTILS
#define UTILS

#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/resource.h>
#include <boost/program_options.hpp>
#include <sdsl/bit_vectors.hpp>

namespace po = boost::program_options;

/*  get peak memory RAM usage   */
inline long get_mem_usage()
{
    struct rusage myusage;

    getrusage(RUSAGE_SELF, &myusage);
    return myusage.ru_maxrss;
}

/*  get time usage   */
inline double get_time_usage()
{
    struct rusage myusage;

    getrusage(RUSAGE_SELF, &myusage);
    double system_time = (double)myusage.ru_stime.tv_sec + (double)myusage.ru_stime.tv_usec / 1000000.0;
    double user_time = (double)myusage.ru_utime.tv_sec + (double)myusage.ru_utime.tv_usec / 1000000.0;
    return system_time + user_time;
}

int read_patterns(std::filesystem::path pattern_file, std::vector<std::string> &patterns)
{

    std::cout << "Reading pattern file...";
    /*  open file   */
    std::ifstream in(pattern_file);
    if (!in.is_open())
        return -1;

    int npatterns = 0;
    std::string line;

    while (std::getline(in, line))
    {
        patterns.push_back(line);
        npatterns++;
    }

    std::cout << "  DONE. Total " << npatterns << std::endl;
    return 0;
}

/*  from given MSA explicitly create EDS under condition of common context of length at least l and save it into eds_file    */
int msa2leds(std::filesystem::path &msa_file, std::filesystem::path &eds_file, size_t l, bool save)
{
    uint64_t N = 0; //  total length given by the sum of all degenerate sets
    uint64_t n = 0; //  number of sets in EDS
    uint64_t m = 0; //  cardinality of EDS = total number of fractioned strings in the sets (sum of the set sizes icluding the "reference one")

    int min_l = std::filesystem::file_size(msa_file);
    int max_l = -1;

    if (eds_file.empty())
    {
        // save it on the same place as msa with differend extension
        eds_file = msa_file;
        eds_file.replace_extension("." + std::to_string(l) + ".s.leds");
    }

    //  read msa_file
    std::ifstream in(msa_file, std::ifstream::binary);
    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << msa_file << std::endl;
        return -1;
    }

    //  get the first "reference" sequence to compare with
    std::string line;
    uint64_t counter = 0;
    std::string first;
    sdsl::bit_vector B;
    size_t i = 0;
    int x = -1; //  to setup the width of fasta

    std::vector<int> start_positions;

    while (std::getline(in, line))
    {
        if (line.empty()) // Skip empty lines
            continue;

        if (line[0] == '>')
        {
            // reset
            if (counter == 1)
            {
                //  initialize bit vector
                B = sdsl::bit_vector(first.size() + 1, 1);
            }
            i = 0;
            counter++;
            start_positions.push_back(in.tellg());
        }
        else if (counter == 1)
        {
            //  first sequence
            first += line;
            if (x == -1)
                x = line.size();
        }
        else
        {
            // compare with the same position in first sequence and switch bit vector if not
            for (size_t j = 0; j < line.size(); j++)
            {
                if (line[j] != first[i] || line[j] == '-')
                {
                    B[i] = 0;
                }
                i++;
            }
        }
    }

    B[first.size()] = B[first.size() - 1] ^ 1;

    // std::cout << first << std::endl;
    // std::cout << B << std::endl;
    // for (size_t j = 0; j < start_positions.size(); j++)
    // {
    //     std::cout << start_positions[j] << ',';
    // }
    //     std::cout <<std::endl;

    sdsl::bit_vector::select_1_type selectOne(&B);
    sdsl::bit_vector::select_0_type selectZero(&B);

    size_t next_one = 0, next_zero = 0;
    size_t zeros = 0, ones = 0;
    size_t cl = 0;

    next_one = selectOne(1);
    next_zero = selectZero(1);

    sdsl::bit_vector H(B.size(), 0);

    i = 0;
    while (i < first.size())
    {
        // std::cout << i << ":" << B[i] << std::endl;

        if (B[i])
        {
            // match
            next_zero = selectZero(zeros + 1);
            // std::cout << "next zero: " << next_zero << ", total zeros: " << zeros << std::endl;

            cl = next_zero - i;
            // std::cout << "run of 1s length: " << cl << std::endl;
            if (cl >= l || i==0 || next_zero ==first.size())
            {
                // std::cout << "here is common part" << std::endl;
                for (size_t j = i; j < next_zero; j++)
                    H[j] = 1;
            }
            ones += cl;
            i = next_zero;
        }
        else
        {
            //  mismatch
            next_one = selectOne(ones + 1);
            // std::cout << "next one: " << next_one << ", total ones: " << ones << std::endl;
            cl = next_one - i;
            // std::cout << "run of 0s length: " << cl << std::endl;
            zeros += cl;
            i = next_one;
        }
    }

    H[first.size()] = H[first.size() - 1] ^ 1;
    // std::cout << H << std::endl;

    std::ofstream out;
    if (save)
    {
        out.open(eds_file, std::ios::out);

        // Check if the file was opened successfully
        if (!out.is_open())
        {
            std::cerr << "Error: Unable to create file." << std::endl;
            return 1;
        }
    }

    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << msa_file << std::endl;
        return -1;
    }

    selectOne = decltype(selectOne)(&H);
    selectZero = decltype(selectZero)(&H);

    next_one = 0;
    next_zero = 0;
    zeros = 0;
    ones = 0;
    cl = 0;
    int tmp;
    char *buffer = new char[first.size()+(first.size()/x)];

    in.clear();
    in.seekg(0);

    std::set<std::string> changes;
    std::string change;

    // std::cout << H << std::endl;

    i = 0;
    while (i < first.size())
    {
        if (H[i])
        {
            // std::cout << "match" << std::endl;
            // match - read cl characters from first sequence
            next_zero = selectZero(zeros + 1);

            cl = next_zero - i;

            ones += cl;
            // Move to specific position
            // std::cout << "moving to: " << start_positions[0] + i + (int)(i / x) << " position" << std::endl;
            tmp = (int)((i % x) + cl) / x;

            in.seekg(start_positions[0] + i + (int)(i / x));

            in.read(buffer, cl + tmp);
            buffer[cl + tmp] = '\0';

            n++;
            m++;
            counter = 0;          

            for (size_t k = 0; k < cl + tmp; k++)
                if (buffer[k] != '\n' && buffer[k] != '-')
                {
                    counter++;
                    if (save)
                        out << buffer[k];
                    N++;
                }
            if (counter > max_l)
                max_l = counter;

            if (counter < min_l && (ones != 0))
                min_l = counter;

            i = next_zero;

        }
        else
        {
            //  mismatch
            next_one = selectOne(ones + 1);

            cl = next_one - i;

            if (save)
                out << '{';

            for (size_t j = 0; j < start_positions.size(); j++)
            {
                // std::cout << i + cl << std::endl;
                in.seekg(start_positions[j] + i + static_cast<int>(i / x));
                tmp = static_cast<int>((i % x) + cl) / x;

                in.read(buffer, cl + tmp);
                buffer[cl + tmp] = '\0';

                // std::cout << first.size() << std::endl;
                // std::cout << cl+tmp << std::endl;

                for (size_t k = 0; buffer[k] != '\0'; ++k)
                // for (size_t k = 0; k < cl + tmp; k++)
                {
                    if (buffer[k] != '\n' && buffer[k] != '-')
                    {
                        // std::cout << buffer[k];
                        change.push_back(buffer[k]);
                        // if (save)
                        //     out << buffer[k];
                        // N++;
                    }
                }
                changes.insert(change);
                change.clear();
                // if (save)
                //     out << ',';
                // m++;
            }
            if (save){
                for(auto ch : changes){
                    out << ch;
                    N+=ch.size();
                    if (ch.empty())
                        N++;                    
                    m++;
                    out << ',';
                }
                changes.clear();
                tmp = out.tellp();
                out.seekp(tmp-1);
            }

            if (save)
                out << '}';

            n++;
            zeros += cl;
            i = next_one;
        }
    }

    // delete[] buffer;

    in.close();
    out.close();

    if (min_l == std::filesystem::file_size(msa_file))
        min_l = -1;

    std::cout << "N: " << N << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "m: " << m << std::endl;
    std::cout << "min_l: " << min_l << std::endl;
    std::cout << "max_l: " << max_l << std::endl;
    return 0;
}

std::set<std::string> cartesian(std::set<std::string> &a, std::set<std::string> &b)
{
    std::set<std::string> res;
    //  cartesian product
    for (auto ai : a)
        for (auto bi : b)
            res.insert(ai + bi);
    if (res.empty())
        return b;

    return res;
}

std::vector<size_t> cartesian(std::vector<size_t> &a, std::vector<size_t> &b)
{
    std::vector<size_t> res;
    //  cartesian product
    for (auto ai : a)
        for (auto bi : b)
            res.push_back(ai + bi);
    if (res.empty())
        return b;

    return res;
}

/*  from given MSA explicitly create EDS under condition of common context of length at least l using cartesian multiplication and save it into eds_file    */
inline int msa2cartleds(std::filesystem::path &msa_file, std::filesystem::path &eds_file, size_t l, bool save)
{
    int N = 0; //  total length given by the sum of all degenerate sets
    int n = 0; //  number of sets in EDS
    int m = 0; //  cardinality of EDS = total number of fractioned strings in the sets (sum of the set sizes icluding the "reference one")

    int min_l = std::filesystem::file_size(msa_file);
    int max_l = -1;

    if (eds_file.empty())
    {
        // save it on the same place as msa with differend extension
        eds_file = msa_file;
        eds_file.replace_extension("." + std::to_string(l) + "c.leds");
    }

    //  read msa_file
    std::ifstream in(msa_file, std::ifstream::binary);
    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << msa_file << std::endl;
        return -1;
    }

    //  get the first "reference" sequence to compare with
    std::string line;
    int counter = 0;
    std::string first;
    sdsl::bit_vector B;
    size_t i = 0;
    int x = -1; //  to setup the width of fasta

    std::vector<int> start_positions;

    while (std::getline(in, line))
    {

        if (line.empty()) // Skip empty lines
            continue;

        if (line[0] == '>')
        {
            // reset
            if (counter == 1)
            {
                //  initialize bit vector
                B = sdsl::bit_vector(first.size() + 1, 1);
            }
            i = 0;
            counter++;
            start_positions.push_back(in.tellg());
        }
        else if (counter == 1)
        {
            //  first sequence
            first += line;
            if (x == -1)
                x = line.size();
        }
        else
        {
            // compare with the same position in first sequence and switch bit vector if not
            for (size_t j = 0; j < line.size(); j++)
            {
                if (line[j] != first[i] || line[j] == '-')
                {
                    B[i] = 0;
                }
                i++;
            }
        }
    }

    B[first.size()] = B[first.size() - 1] ^ 1;

    // std::cout << B << std::endl; 
    sdsl::bit_vector::select_1_type selectOne(&B);
    sdsl::bit_vector::select_0_type selectZero(&B);

    size_t next_one = 0, next_zero = 0;
    size_t zeros = 0, ones = 0;
    size_t cl = 0;

    next_one = selectOne(1);
    next_zero = selectZero(1);

    std::ofstream out;
    if (save)
    {
        out.open(eds_file, std::ios::out);

        // Check if the file was opened successfully
        if (!out.is_open())
        {
            std::cerr << "Error: Unable to create file." << std::endl;
            return 1;
        }
    }

    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << msa_file << std::endl;
        return -1;
    }

    next_one = 0;
    next_zero = 0;
    zeros = 0;
    ones = 0;
    cl = 0;
    int tmp;
    std::set<std::string> changes_o;
    std::set<std::string> changes_n;
    char *buffer;
    buffer = new char[first.size()];
    line.clear();

    in.clear();

    i = 0;
    while (i < first.size())
    {
        if (B[i])
        {
            // match - read cl characters from first sequence
            next_zero = selectZero(zeros + 1);

            cl = next_zero - i;
            ones += cl;

            // Move to specific position
            tmp = (int)((i % x) + cl) / x;
            in.seekg(start_positions[0] + i + (int)(i / x));

            in.read(buffer, cl + tmp);
            buffer[cl + tmp] = '\0';

            for (size_t k = 0; k < cl + tmp; k++)
                if (buffer[k] != '\n' && buffer[k] != '-' && buffer[k] != '\0')
                {
                    line.push_back(buffer[k]);
                }

            if (next_zero == first.size())
            {
                break;
            }
            

            if (i==0 && !line.empty())
            {
                n++;
                m++;
                out << line;
                N+=line.size();

                if (cl > max_l)
                {
                    max_l = cl;
                }
                if (cl < min_l)
                {
                    min_l = cl;
                }
            }else if (cl < l)
            {   //  neproddo not output common part, but prolong changes
                for (auto ch : changes_o)
                {
                    changes_n.insert(ch + line);
                }
            }
            else
            {   //  print changes + common part
                if (cl > max_l)
                {
                    max_l = cl;
                }
                if (cl < min_l)
                {
                    min_l = cl;
                }
                
                n += 2;
                m += changes_o.size() +1; //   +1 is for the common part

                //  flush data into file
                if (save)
                {
                    out << '{';
                    for (const std::string &str : changes_o)
                    {
                        if (str.empty())
                            N++;
                        
                        N += str.size();
                        out << str; // Write each string followed by a newline
                        out << ',';
                    }
                    tmp = out.tellp();
                    out.seekp(tmp - 1);

                    if (tmp != 1)
                        out << '}';

                    out << line;
                    N+=line.size();

                }
                else
                {
                    for (const std::string &str : changes_o)
                        N += str.size();
                }
                changes_o.clear();
                line.clear();
            }
            i = next_zero;
            changes_o = changes_n;
            changes_n.clear();
            line.clear();
        }
        else
        {
            //  mismatch
            next_one = selectOne(ones + 1);
            cl = next_one - i;

            for (size_t j = 0; j < start_positions.size(); j++)
            {
                in.seekg(start_positions[j] + i + (int)(i / x));
                tmp = (int)((i % x) + cl) / x;

                in.read(buffer, cl + tmp);
                buffer[cl + tmp] = '\0';

                line.clear();
                for (size_t k = 0; k < cl + tmp; k++)
                    if (buffer[k] != '\n' && buffer[k] != '-' && buffer[k] != '\0')
                    {
                        line.push_back(buffer[k]);
                    }

                changes_n.insert(line);
            }

            changes_o = cartesian(changes_o, changes_n);
            changes_n.clear();

            zeros += cl;
            i = next_one;
            line.clear();
        }
    }


    if (!changes_o.empty())
    {
        //  flush data into file
        n += 1;
        m += changes_o.size(); //   +1 is for the common part
        if (save)
        {
            out << '{';
            for (const std::string &str : changes_o)
            {

                        if (str.empty())
                            N++;
                N += str.size();
                out << str; // Write each string followed by a newline
                out << ',';
            }
            tmp = out.tellp();
            out.seekp(tmp - 1);
            if (tmp != 1)
                out << '}';
            
            if (!line.empty()){
                out << line;
                N+=line.size();
                n++;
                m++;
            }
            if (line.size() > max_l)
                {
                    max_l = line.size();
                }
            if (line.size() < min_l)
                {
                    min_l = line.size();
                }
        }
        else
        {
            for (const std::string &str : changes_o)
                N += str.size();
        }
        changes_o.clear();
    }
    delete[] buffer;
    buffer = nullptr;

    in.close();
    out.close();

    if (min_l == std::filesystem::file_size(msa_file))
        min_l = -1;

    std::cout << "N: " << N << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "m: " << m << std::endl;
    std::cout << "min_l: " << min_l << std::endl;
    std::cout << "max_l: " << max_l << std::endl;

    return 0;
}

/*  from given EDS create EDS under condition of common context of length at least l using cartesian multiplication and save it into eds_file    */
inline int eds2leds(std::filesystem::path &original_file, std::filesystem::path &output_file, size_t l)
{
    uint64_t N = 0; //  total length given by the sum of all degenerate sets
    uint64_t n = 0; //  number of sets in EDS
    uint64_t m = 0; //  cardinality of EDS = total number of fractioned strings in the sets (sum of the set sizes icluding the "reference one")

    uint64_t min_l = std::filesystem::file_size(original_file);
    uint64_t max_l = 0;

    if (output_file.empty())
    {
        // save it on the same place as msa with differend extension
        output_file = original_file;
        output_file.replace_extension("." + std::to_string(l) + ".leds");
    }

    std::ifstream in(original_file, std::ifstream::binary);
    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << original_file << std::endl;
        return -1;
    }

    std::ofstream out;
    out.open(output_file, std::ios::out); // 3+3+1+2+3+0.5

    // Check if the file was opened successfully
    if (!out.is_open())
    {
        std::cerr << "Error: Unable to create file." << std::endl;
        return 1;
    }

    char ch;
    bool count = true;
    size_t cl = 0;
    bool first = true;
    std::string tmp = "";
    std::string line = "";
    size_t tmpi;
    std::set<std::string> changes_n;

    auto now = std::chrono::high_resolution_clock::now();
    auto duration = now.time_since_epoch();
    auto nanosecs = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
    std::stringstream ss;
    ss << nanosecs;
    ss.str();

    std::filesystem::path tmp_input = "./tmp_o." + ss.str() + ".tmp";
    std::filesystem::path tmp_output = "./tmp_n." + ss.str() + "tmp";
    bool merge = false;
    bool first_set = true;

    std::ofstream xy;
    std::ifstream tmp_in;
    std::ofstream tmp_out;

    while (in.get(ch))
    {
        // Process each character
        switch (ch)
        {
        case '{':
            //  reset counter and stop counting
            if (cl < l && !first)
            {
                //  append reference to changes
                if (cl != 0)
                {
                    tmp_in.open(tmp_input);
                    tmp_out.open(tmp_output);
                    while (std::getline(tmp_in, tmp))
                    {
                        tmp += line;
                        // std::cout << tmp << std::endl;
                        tmp_out << tmp << std::endl;
                    }
                    tmp_in.close();
                    tmp_out.close();
                    remove(tmp_input.c_str());
                    rename(tmp_output.c_str(), tmp_input.c_str());
                }
                merge = true;
            }
            else
            {
                //  flush changes
                out << '{';
                tmp_in.open(tmp_input);
                while (std::getline(tmp_in, tmp))
                {
                    out << tmp;
                    out << ',';
                    N += tmp.size();
                    m++;
                }
                tmp_in.close();
                tmpi = out.tellp();
                out.seekp(tmpi - 1);
                if (tmpi != 1)
                    out << '}';

                xy.open(tmp_input, std::ios::trunc);
                xy.close();
                n++;

                merge = false;
                out << line;
                N += cl;
                n++;
                m++;
            }

            if (cl > max_l && !merge)
                max_l = cl;

            if (cl < min_l && !first && !merge)
                min_l = cl;

            cl = 0;
            count = false;
            first = false;
            line.clear();
            break;
        case '}':
            changes_n.insert(line);

            if (merge)
            {
                tmp_in.open(tmp_input);
                tmp_out.open(tmp_output);
                while (std::getline(tmp_in, tmp))
                {
                    for (auto n : changes_n)
                    {
                        // std::cout << tmp << n << std::endl;
                        tmp_out << tmp;
                        tmp_out << n << std::endl;
                    }
                }
                tmp_in.close();
                tmp_out.close();
                remove(tmp_input.c_str());
                rename(tmp_output.c_str(), tmp_input.c_str());

                // changes_o = cartesian(changes_o, changes_n);
            }
            else // if (first_set)
            {
                tmp_out.open(tmp_output);
                for (auto n : changes_n)
                {
                    // std::cout << n << std::endl;
                    tmp_out << n << std::endl;
                }
                first_set = false;
                tmp_out.close();
                remove(tmp_input.c_str());
                rename(tmp_output.c_str(), tmp_input.c_str());
            }
            count = true;
            line.clear();
            changes_n.clear();
            break;
        case ',':
            changes_n.insert(line);
            line.clear();
            break;
        default:
            if (count)
                cl++;
            line.push_back(ch);
            break;
        }
    }
    if (ch == '}')
    {
        n--;
        m--;
    }

    //  flush data into file
    tmp_in.open(tmp_input);
    out << '{';
    while (std::getline(tmp_in, tmp))
    {
        out << tmp;
        out << ',';
        N += tmp.size();
        m++;
    }
    tmp_in.close();

    tmpi = out.tellp();
    out.seekp(tmpi - 1);
    if (tmpi != 1)
        out << '}';
    out << line;

    in.close();
    out.close();
    if (min_l == std::filesystem::file_size(original_file))
        min_l = -1;

    std::cout << "N: " << N << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "m: " << m << std::endl;
    std::cout << "min_l: " << min_l << std::endl;
    std::cout << "max_l: " << max_l << std::endl;

    remove(tmp_input.c_str());
    return 0;
}

/*  from given EDS create EDS under condition of common context of length at least l using  */
inline int eds2leds(std::filesystem::path &original_file, size_t l)
{
    uint64_t N = 0; //  total length given by the sum of all degenerate sets
    uint64_t n = 0; //  number of sets in EDS
    uint64_t m = 0; //  cardinality of EDS = total number of fractioned strings in the sets (sum of the set sizes icluding the "reference one")
    uint64_t min_l = std::filesystem::file_size(original_file);
    uint64_t max_l = 0;

    std::ifstream in(original_file, std::ifstream::binary);
    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << original_file << std::endl;
        return -1;
    }

    char ch;
    bool count = true;
    size_t cl = 0;
    bool first = true;
    std::string line = "";
    std::pair<uint64_t, uint64_t> Si = {0, 0};
    std::pair<uint64_t, uint64_t> Sj = {0, 0};

    bool merge = false;

    while (in.get(ch))
    {
        // Process each character
        switch (ch)
        {
        case '{':
            //  reset counter and stop counting
            if (cl < l && !first)
            {
                Si.second += (Si.first * line.size());
                merge = true;
            }
            else
            {
                //  ref
                if (cl != 0)
                {
                    N += cl;
                    n++;
                    m++;
                }
                merge = false;
            }

            if (cl > max_l && !merge)
                max_l = cl;

            if (cl < min_l && !first && !merge)
                min_l = cl;

            cl = 0;
            count = false;
            first = false;
            line.clear();
            break;
        case '}':

            Sj.first++;
            Sj.second += line.size();
            // std::cout << Sj.first << "," << Sj.second << std::endl;
            // start counting
            if (merge || Si.first == 0)
            {
                if (Si.first == 0)
                {
                    Si.first = Sj.first;
                    Si.second = Sj.second;
                }
                else
                {
                    Si.second = (Sj.first * Si.second) + (Si.first * Sj.second);
                    Si.first = Si.first * Sj.first;
                }
            }
            else
            {
                //  flush
                n++;
                N += Si.second;
                m += Si.first;
                Si.first = 0;
                Si.second = 0;
            }
            count = true;
            line.clear();
            Sj.first = 0;
            Sj.second = 0;

            // std::cout << Si.first << "," << Si.second << std::endl;
            break;
        case ',':
            Sj.first++;
            Sj.second += line.size();
            line.clear();
            break;
        default:
            if (count)
                cl++;
            line.push_back(ch);
            break;
        }
    }

    if (Si.first != 0)
    {
        //  flush data into file
        n++;
        N += Si.second;
        m += Si.first;
    }

    in.close();
    if (min_l == std::filesystem::file_size(original_file))
        min_l = -1;

    std::cout << "N: " << N << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "m: " << m << std::endl;
    std::cout << "min_l: " << min_l << std::endl;
    std::cout << "max_l: " << max_l << std::endl;

    return 0;
}

int count_eds(std::filesystem::path &eds_file)
{
    int N = 0; //  total length given by the sum of all degenerate sets
    int n = 0; //  number of sets in EDS
    int m = 0; //  cardinality of EDS = total number of fractioned strings in the sets (sum of the set sizes icluding the "reference one")
    int min_l = std::filesystem::file_size(eds_file);
    int max_l = -1;

    //  read msa_file
    std::ifstream in(eds_file, std::ifstream::binary);
    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << eds_file << std::endl;
        return -1;
    }

    char ch;
    bool count = true;
    size_t cl = 0;
    bool first = true;
    size_t line_size = 0;
    bool last = false;

    while (in.get(ch))
    {
        // Process each character
        // std::cout << "Character: " << ch << ", cl: " <<cl<< std::endl;
        switch (ch)
        {
        case '{':
            last = true;
            m++;
            n++;

            if (cl > max_l)
                max_l = cl;

            if (cl < min_l && !first)
                min_l = cl;

            if (first && N > 0)
            {
                n++;
                m++;
            }

            cl = 0;
            count = false;
            first = false;
            break;
        case '}':
            if (last){
                N++;
                last = false;
            }
            // start counting
            m++;
            n++;
            count = true;
            break;
        case '-':
            std::cout << "Invalid character: " << ch << std::endl;
            return -1;
            break;
        case ',':
            if (last){
                N++;
                last = false;
            }
            last = true;
            m++;
            break;
        default:
            last = false;
            N++;
            if (count)
                cl++;
            break;
        }
    }
    if (ch == '}')
    {
        n--;
        m--;
    }

    in.close();

    if (min_l == std::filesystem::file_size(eds_file))
        min_l = -1;

    std::cout << "N: " << N << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "m: " << m << std::endl;
    std::cout << "min_l: " << min_l << std::endl;
    std::cout << "max_l: " << max_l << std::endl;

    return 0;
}

// Function to read EDS from a file
std::string readEDSFromFile(const std::filesystem::path &filePath)
{
    std::ifstream file(filePath);
    if (!file)
    {
        std::cerr << "Error: Unable to open file for reading.\n";
        return "";
    }

    std::string eds, line;
    while (getline(file, line))
    {
        eds += line;
    }
    file.close();
    return eds;
}

// Function to extract a single string from the entire EDS with random choices
std::string generateRandomFullLengthEDSPattern(const std::string &eds)
{
    std::string result;
    size_t pos = 0;
    std::mt19937 rng(std::random_device{}()); // Create a random number generator

    while (pos < eds.size())
    {
        if (eds[pos] == '{')
        {
            size_t closePos = eds.find('}', pos);
            if (closePos == std::string::npos)
            {
                std::cerr << "Malformed EDS: No closing '}' found." << std::endl;
                break;
            }

            std::string options = eds.substr(pos + 1, closePos - pos - 1);
            std::vector<std::string> choices;
            std::string option;
            std::stringstream ss(options);
            while (getline(ss, option, ','))
            {
                choices.push_back(option);
            }

            if (!choices.empty())
            {
                std::uniform_int_distribution<size_t> dist(0, choices.size() - 1);
                result += choices[dist(rng)]; // Select a random option
            }

            pos = closePos + 1;
        }
        else
        {
            result += eds[pos];
            pos++;
        }
    }

    return result;
}

// Function to generate a random substring of a specified length from a given string
std::string getRandomSubstring(const std::string &inputString, int targetLength)
{
    if (inputString.length() < targetLength)
    {
        std::cerr << "Error: Target length is greater than the input string length." << std::endl;
        return "";
    }

    std::random_device rd;                                                              // Non-deterministic random number generator
    std::mt19937 rng(rd());                                                             // Seed the generator
    std::uniform_int_distribution<size_t> dist(0, inputString.length() - targetLength); // Ensure valid start range

    size_t startPos = dist(rng);
    return inputString.substr(startPos, targetLength);
}

void generate_patterns(std::string eds, int number, int length, std::vector<std::string> &patterns)
{
    std::string input;
    std::string p;

    // Initialize the random number generator
    std::mt19937 rng(std::time(nullptr));

    for (size_t i = 0; i < number; i++)
    {
        p = getRandomSubstring(generateRandomFullLengthEDSPattern(eds), length);
        if (!p.empty())
        {
            patterns.push_back(p);
        }
    }
}

std::string check_result(const std::string &eds, int pos, size_t l, const std::vector<int> &changeIndices)
{
    int current_change_i = 0;
    int j = 0;
    std::string result;
    bool ref = true;

    int change;

    for (int i = 0; i < eds.size(); i++)
    {
        // ACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGT
        // if (j == changeIndices.size() && i > start)
        // {
        // break;
        // }

        if (changeIndices.empty() || j >= changeIndices.size())
        {

            change = 0;
        }
        else
        {

            change = changeIndices[j];
        }

        switch (eds[i])
        {
        case '{':
            ref = false;
            current_change_i++;
            break;
        case '}':
            if (current_change_i == change)
            {
                j++;
            }
            ref = true;
            break;
        case ',':
            if (current_change_i == change)
            {
                j++;
            }
            current_change_i++;
            break;

        default:
            if (current_change_i == change)
            {
                result.push_back(eds[i]);
            }
            else if (ref)
            {
                result.push_back(eds[i]);
            }
            break;
        }
    }

    // std::cout << result << std::endl;
    return result.substr(pos - 1, l);
}

#endif // UTILS