#include "parser.h"

long get_mem_usage()
{
    struct rusage myusage;

    getrusage(RUSAGE_SELF, &myusage);
    return myusage.ru_maxrss;
}

double get_time_usage()
{
    struct rusage myusage;

    getrusage(RUSAGE_SELF, &myusage);
    double system_time = (double)myusage.ru_stime.tv_sec + (double)myusage.ru_stime.tv_usec / 1000000.0;
    double user_time = (double)myusage.ru_utime.tv_sec + (double)myusage.ru_utime.tv_usec / 1000000.0;
    return system_time + user_time;
}

// void print_v(std::vector<unsigned> v)
// {
//     for (const auto &i : v)
//     {
//         std::cout << i << " ";
//     }
//     std::cout << std::endl;
// }

namespace bio_fmi
{
    Parser::Parser(unsigned context_length) : BLOCK_SIZE(32), context_length(context_length){}

    Parser::~Parser() {}

    // void Parser::processBlock(std::vector<char> buffer)
    // {
    //     std::cout << buffer << std::endl;

    //     //  first cross
    //     std::vector<unsigned> start_b;
    //     std::vector<unsigned> end_b;
    //     std::vector<unsigned> number_of_changes;
    //     std::vector<unsigned> lenght_of_changes;
    //     unsigned len;

    //     // find starts of degenerate symbols
    //     for(auto i = 0; i != buffer.size(); i++){
    //         len++;
    //         switch (buffer[i])
    //         {
    //         case '{':
    //             start_b.push_back(i);
    //             number_of_changes.push_back(number_of_changes.empty()?0:number_of_changes.back());  //cumulative
    //             len = 0;
    //             break;
    //         case '}':
    //             end_b.push_back(i);
    //             lenght_of_changes.push_back(len);
    //             break;
    //         case ',':
    //             number_of_changes.back()++;
    //             lenght_of_changes.push_back(len);
    //             len = 0;
    //             break;
    //         default:
    //             break;
    //         }
            
    //     }

    //     unsigned start_len = 0;
    //     unsigned end_len = 0;
    //     if (!start_b.empty())
    //          start_len = start_b.front();
    //     if (!end_b.empty())
    //         end_len = buffer.size() - end_b.back();

    //     std::cout << "start_len: " << start_len << std::endl;
    //     std::cout << "end_len: " << end_len << std::endl;
    //     print_v(start_b);
    //     print_v(end_b);
    //     print_v(number_of_changes);
    //     print_v(lenght_of_changes);

    //     //  second crossing
    //     for(auto i = 0; i != buffer.size(); i++){
            
    //     }

    //     //flush
    // }

    int Parser::parseEDS(std::filesystem::path input_path)
    {
        //  vystupem chceme vsechny potrebne struktury a hlavni retezec
        //  nejdriv zkusime cely text najednou nacist do pameti

        // metadata zapisujeme na konci bloku, bloky chceme co nejvesti, abychom redukovali mnozstvi IO operaci

        //  open stream
        std::ifstream input_file(input_path, std::ios::binary);

        if (!input_file)
        {
            std::cerr << "Failed to open the input file." << std::endl;
            return 1;
        }

        buffer_.resize(BLOCK_SIZE);

        int counter = 0;
        // Read and process data in blocks
        while (input_file.read(buffer_.data(), BLOCK_SIZE))
        {
            std::size_t bytesRead = input_file.gcount(); // Number of bytes read
            counter++;
            // Process the block
            processBlock(std::vector<char>(buffer_.begin(), buffer_.begin() + bytesRead));
        }
        std::cout << "counter: " << counter << std::endl;

        // Process the remaining data (if any)
        std::size_t bytesRead = input_file.gcount(); // Number of bytes read after the last block
        if (bytesRead > 0)
        {
            processBlock(std::vector<char>(buffer_.begin(), buffer_.begin() + bytesRead));
        }

        input_file.close();
    }
}