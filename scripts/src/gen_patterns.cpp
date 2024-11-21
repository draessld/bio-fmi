#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/resource.h>
#include <boost/program_options.hpp>
#include "utils/eds.hpp"

using namespace std;
namespace po = boost::program_options;

std::string usage = "";
std::string desc = "";
unsigned int l;
unsigned int t;
std::filesystem::path in_file = "";
std::filesystem::path out_file = "";

int handle_parameters(int argc, const char **argv)
{
    po::options_description desc("Allowed options");

    desc.add_options()("help", "produce help message")
    ("input,i", po::value<std::filesystem::path>(&in_file)->required(), "")
    ("output,o", po::value<std::filesystem::path>(&out_file), "")
    ("number,t", po::value<unsigned int>(&t), "number of patterns")
    ("size,l", po::value<unsigned int>(&l), "length of patterns");

    po::positional_options_description posOptions;
    posOptions.add("input", 1);

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

    //  setup out file
    if (out_file.empty()){
        out_file = in_file;
        out_file.replace_extension(std::to_string(t)+"_"+std::to_string(l)+".patterns");
    }
    cout << out_file <<endl;

    ifstream ifs(in_file);
    ofstream ofs(out_file);
    if (!ifs)  // If the file could not be opened, treat input as a sequence
    {
        cout << "Error: File was not found!" << endl;
        return 1;
    }
    
    EDS eds(ifs);
    // eds.stats();
    eds.gen_patterns(ofs,t,l);
    
    return 0;
}