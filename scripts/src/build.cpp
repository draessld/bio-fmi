#include <iostream>
#include <fstream>
#include <sstream>

#include "utils/eds.h"
#include "utils/index.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace bio_fmi;
namespace po = boost::program_options;

std::string usage = "";
std::string desc = "";
std::string method;
unsigned int l;
std::filesystem::path in_file;
std::filesystem::path out_file;


int handle_parameters(int argc, const char **argv)
{
    po::options_description desc("Allowed options");

    desc.add_options()("help", "produce help message")
    ("input,i", po::value<std::filesystem::path>(&in_file)->required(), "")
    ("output,o", po::value<std::filesystem::path>(&out_file), "")
    ("context_length,l", po::value<unsigned int>(&l), "length of chunk and stored context");
    
    po::positional_options_description posOptions;

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(posOptions).run(), vm);

        if (vm.count("help"))
        {
            std::cout << "Usage: " << argv[0] << " " << usage << std::endl
                      << std::endl;
            std::cout << desc << std::endl;

            return 1;
        }

        if ((vm.count("eds-file") == 0) && (vm.count("msa-file") == 0))
        {
            std::cout << "Error: No input given" << std::endl;

            std::cout << "Usage: " << argv[0] << " " << usage << std::endl
                      << std::endl;
            return -1;
        }


        po::notify(vm);
    }
    catch (const po::error &e)
    {
        std::cerr << "Usage: " << argv[0] << " " << usage << std::endl
                  << std::endl;
        std::cerr << desc << std::endl;

        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

int main(int argc, char const *argv[])
{
    /*  Parse input command line */
    int parameter_handle_result = handle_parameters(argc, argv);
    if (parameter_handle_result == -1)
    {
        std::cout << "Error while reading parameters\n"
                  << std::endl;
        return EXIT_FAILURE;
    }
    else if (parameter_handle_result == 1)
    {
        return EXIT_SUCCESS;
    }

    ifstream ifs(in_file);
    if (!ifs)  // If the file could not be opened, treat input as a sequence
    {
        std::cout << "Error while reading input_file\n"
                  << std::endl;
        return EXIT_FAILURE;
    }
    else  // If the file can be opened, treat input as a file
    {
        EDS eds(ifs);
        eds.stats();
    }
    
    return 0;
}