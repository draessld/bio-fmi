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
inline int get_l_eds(std::filesystem::path &msa_file, std::filesystem::path &eds_file, size_t l)
{
    int N = 0; //  total length given by the sum of all degenerate sets
    int n = 0; //  number of sets in EDS
    int m = 0; //  cardinality of EDS = total number of fractioned strings in the sets (sum of the set sizes icluding the "reference one")

    if (eds_file.empty())
    {
        // save it on the same place as msa with differend extension
        eds_file = msa_file;
        eds_file.replace_extension("." + std::to_string(l) + ".eds");
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
            if (cl >= l)
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

    std::ofstream out(eds_file, std::ios::out);

    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << msa_file << std::endl;
        return -1;
    }

    // Check if the file was opened successfully
    if (!out.is_open())
    {
        std::cerr << "Error: Unable to create file." << std::endl;
        return 1;
    }

    selectOne = decltype(selectOne)(&H);
    selectZero = decltype(selectZero)(&H);

    next_one = 0;
    next_zero = 0;
    zeros = 0;
    ones = 0;
    cl = 0;
    int tmp;
    std::set<std::string> changes;
    char *buffer;
    buffer = new char[first.size()];

    in.clear();

    i = 0;
    while (i < first.size())
    {
        if (H[i])
        {
            // match - read cl characters from first sequence
            next_zero = selectZero(zeros + 1);

            cl = next_zero - i;
            ones += cl;

            // Move to specific position
            // std::cout << "moving to: " << start_positions[0]+i+(int)(i/x) << " position" << std::endl;
            tmp = (int)((i % x) + cl) / x;

            in.seekg(start_positions[0] + i + (int)(i / x));

            in.read(buffer, cl + tmp);
            buffer[cl + tmp] = '\0';

            n++;
            m++;
            for (size_t k = 0; k < cl + tmp; k++)
                if (buffer[k] != '\n' && buffer[k] != '-' && buffer[k] != '\0')
                {
                    out << buffer[k];
                    N++;
                }

            i = next_zero;
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

                changes.insert(line);
            }

            n++;
            m += changes.size();
            if (changes.size() > 1)
            {
                out << '{';

                for (const std::string &str : changes)
                {
                    N += str.size();
                    out << str; // Write each string followed by a newline
                    out << ',';
                }
                tmp = out.tellp();
                out.seekp(tmp - 1);

                out << '}';
            }
            zeros += cl;
            i = next_one;

            changes.clear();
        }
    }

    delete[] buffer;
    buffer = nullptr;

    in.close();
    out.close();

    std::cout << "N: " << N << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "m: " << m << std::endl;
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

/*  from given MSA explicitly create EDS under condition of common context of length at least l using cartesian multiplication and save it into eds_file    */
inline int get_cartl_eds(std::filesystem::path &msa_file, std::filesystem::path &eds_file, size_t l)
{
    int N = 0; //  total length given by the sum of all degenerate sets
    int n = 0; //  number of sets in EDS
    int m = 0; //  cardinality of EDS = total number of fractioned strings in the sets (sum of the set sizes icluding the "reference one")

    if (eds_file.empty())
    {
        // save it on the same place as msa with differend extension
        eds_file = msa_file;
        eds_file.replace_extension("." + std::to_string(l) + ".eds");
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

    sdsl::bit_vector::select_1_type selectOne(&B);
    sdsl::bit_vector::select_0_type selectZero(&B);

    size_t next_one = 0, next_zero = 0;
    size_t zeros = 0, ones = 0;
    size_t cl = 0;

    next_one = selectOne(1);
    next_zero = selectZero(1);

    std::ofstream out(eds_file, std::ios::out);

    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << msa_file << std::endl;
        return -1;
    }

    // Check if the file was opened successfully
    if (!out.is_open())
    {
        std::cerr << "Error: Unable to create file." << std::endl;
        return 1;
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
            line.clear();
            in.seekg(start_positions[0] + i + (int)(i / x));

            in.read(buffer, cl + tmp);
            buffer[cl + tmp] = '\0';

            for (size_t k = 0; k < cl + tmp; k++)
                if (buffer[k] != '\n' && buffer[k] != '-' && buffer[k] != '\0')
                {
                    line.push_back(buffer[k]);
                }

            if ((cl < l) && (i != 0))
            {
                for (auto ch : changes_o)
                {
                    changes_n.insert(ch + line);
                }
            }
            else
            {
                //  flush data into file
                out << '{';

                n += 2;
                m += changes_o.size() + 1; //   +1 is for the common part
                for (const std::string &str : changes_o)
                {
                    N += str.size();
                    out << str; // Write each string followed by a newline
                    out << ',';
                }
                tmp = out.tellp();
                out.seekp(tmp - 1);

                if (tmp != 1)
                    out << '}';

                out << line;
                changes_o.clear();
            }
            i = next_zero;
            changes_o = changes_n;
            changes_n.clear();
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
        }
    }

    if (!changes_o.empty())
    {
        //  flush data into file
        out << '{';
        n += 2;
        m += changes_o.size() + 1; //   +1 is for the common part
        for (const std::string &str : changes_o)
        {
            N += str.size();
            out << str; // Write each string followed by a newline
            out << ',';
        }
        tmp = out.tellp();
        out.seekp(tmp - 1);
        if (tmp != 1)
            out << '}';
        out << line;
        changes_o.clear();
    }
    

    delete[] buffer;
    buffer = nullptr;

    in.close();
    out.close();

    std::cout << "N: " << N << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "m: " << m << std::endl;
    
    return 0;
}

int validate_eds(std::filesystem::path &eds_file, int l)
{
    int N = 0; //  total length given by the sum of all degenerate sets
    int n = 0; //  number of sets in EDS
    int m = 0; //  cardinality of EDS = total number of fractioned strings in the sets (sum of the set sizes icluding the "reference one")

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

    while (in.get(ch))
    {
        // Process each character
        // std::cout << "Character: " << ch << ", cl: " <<cl<< std::endl;
        switch (ch)
        {
        case '{':
            m++;
            n++;
            //  reset counter and stop counting
            if (cl < l && !first){
                std::cout << "Context is not long enough!!" << std::endl;
                // return -1;
            }
            cl = 0;

            if (first && N>0)
            {
                n++;
                m++;
            }
            
            count = false;
            first = false;
            break;
        case '}':
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
            m++;
            break;
        default:
            N++;
            if (count)
                cl++;
            break;
        }
    }
    if (ch == '}'){
        n--;
        m--;
    }

    
    in.close();

    std::cout << "N: " << N << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "m: " << m << std::endl;
    
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