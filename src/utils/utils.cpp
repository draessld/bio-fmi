#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/resource.h>
#include <boost/program_options.hpp>
#include <sdsl/bit_vectors.hpp>

namespace po = boost::program_options;


/*  get peak memory RAM usage   */
long get_mem_usage()
{
    struct rusage myusage;

    getrusage(RUSAGE_SELF, &myusage);
    return myusage.ru_maxrss;
}

/*  get time usage   */
double get_time_usage()
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

    unsigned npatterns = 0;
    std::string line;

    while (std::getline(in,line)){
        patterns.push_back(line);
        npatterns++;
    }

    std::cout << "  DONE. Total " << npatterns << std::endl;
    return 0;
}


int transfer_msa_to_eds(std::filesystem::path &msa_file, std::filesystem::path &eds_file,size_t l)
{
    if (eds_file.empty())
    {
        // save it on the same place as msa with differend extension
        eds_file = msa_file;
        eds_file.replace_extension("." + std::to_string(l) + ".eds");
    }

    //  read msa_file
    std::ifstream in(msa_file,std::ifstream::binary);
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
            if (x==-1)
                x=line.size();
            
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

    sdsl::bit_vector H(B.size(),0);

    i=0;
    while( i < first.size())
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
    if (!out.is_open()) {
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
    char * buffer;

    in.clear();

    i=0;
    while( i < first.size())
    {
        // std::cout << i << ":" << H[i] << std::endl;

        if (H[i])
        {
            // match - read cl characters from first sequence
            next_zero = selectZero(zeros + 1);
            // std::cout << "next zero: " << next_zero << ", total zeros: " << zeros << std::endl;

            cl = next_zero - i;
            // std::cout << "reading " << cl << "characters" << std::endl;
            ones += cl;

            // Move to specific position
            // std::cout << "moving to: " << start_positions[0]+i+(int)(i/x) << " position" << std::endl;
            tmp = (int)((i%x)+cl)/x;
            
            in.seekg(start_positions[0]+i+(int)(i/x)); 
            // std::cout << "gaps: " <<tmp <<std::endl;
            
            buffer = new char [cl+tmp+1];
            in.read(buffer,cl+tmp);
            buffer[cl+tmp] = '\0';
            // std::cout << buffer << std::endl;

            for(size_t k=0; k<cl+tmp;k++)
                if (buffer[k] !='\n' && buffer[k] !='-') {
                     out << buffer[k];
                }

            i = next_zero;
        }
        else
        {
            //  mismatch
            next_one = selectOne(ones + 1);
            // std::cout << "next one: " << next_one << ", total ones: " << ones << std::endl;
            cl = next_one - i;

            for (size_t j = 0; j < start_positions.size(); j++)
            {
                in.seekg(start_positions[j]+i+(int)(i/x)); 
                tmp = (int)((i%x)+cl)/x;

                buffer = new char [cl+tmp+1];
                in.read(buffer,cl+tmp);
                buffer[cl+tmp] = '\0';
                // std::cout << buffer << std::endl;
                
                line.clear();
                for(size_t k=0; k<cl+tmp;k++)
                    if (buffer[k] !='\n' && buffer[k] !='-') {
                        // std::cout << buffer[k] << ',';
                        line.push_back(buffer[k]);
                    }

                changes.insert(line);
                
            }

            if (changes.size()>1)
            {
                out << '{';

                for (const std::string& str : changes) {
                    out << str; // Write each string followed by a newline
                    out << ',';
                }
                tmp = out.tellp();
                out.seekp(tmp-1);

                out << '}';
            }
            zeros += cl;
            i = next_one;

            changes.clear();
        }
    }

    in.close();
    out.close();
    return 0;
}

int validate_eds(std::filesystem::path &eds_file, unsigned l){
    //  read msa_file
    std::ifstream in(eds_file,std::ifstream::binary);
    if (!in.is_open())
    {
        std::cerr << "Error: Unable to open file " << eds_file << std::endl;
        return -1;
    }
    
    char ch;
    bool count=true;
    size_t cl = 0;
    bool first=true;

    while (in.get(ch)) {
        // Process each character
        // std::cout << "Character: " << ch << ", cl: " <<cl<< std::endl;
        if (ch == '{'){
            //  reset counter and stop counting
            if(cl < l && !first)
                return -1;
            cl = 0;
            count = false;
            first = false;
        }
        if(ch == '-'){
            std::cout << "Invalid character: " << ch << std::endl;
            return -1;
        }
        if (count)
            cl++;

        if (ch == '}')
            // start counting
            count = true;
    }
    in.close();
    return 0;
}
