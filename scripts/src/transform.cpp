#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/resource.h>
#include <boost/program_options.hpp>
#include "utils/utils.hpp"

using namespace std;
namespace po = boost::program_options;

std::string usage = "";
std::string desc = "";
std::string method;
unsigned int l;
std::filesystem::path in_file;
std::filesystem::path out_file;

int check_extension(std::string ext){
    if (ext == ".msa"){
        return 1;
    }else if(ext == ".vcf"){
        return 2;
    }else{
        return 0;
    }
}

int handle_parameters(int argc, const char **argv)
{
    po::options_description desc("Allowed options");

    desc.add_options()("help", "produce help message")
    ("input,i", po::value<std::filesystem::path>(&in_file)->required(), "")
    ("output,o", po::value<std::filesystem::path>(&out_file), "")
    ("context_length,l", po::value<unsigned int>(&l), "length of chunk and stored context, default 5")
    ("method", po::value<std::string>(&method)->default_value("linear"), "specify method of creating eds");

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

        if ((method != "cartesian") && (method != "linear")) {
            std::cerr << "Error: Invalid method. --method={cartesian, linear}\n";
        return 1;
    }
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
        if (l!=0)
        {
            out_file.replace_extension(std::to_string(l)+".eds");
        }else{
            out_file.replace_extension(".eds");
        }   
    }
    cout << out_file <<endl;

    ifstream ifs(in_file);
    ofstream ofs(out_file);
    if (!ifs)  // If the file could not be opened, treat input as a sequence
    {
        cout << "Error: File was not found!" << endl;
        return 1;
    }
    
    switch (check_extension(in_file.extension()))
    {
    case 1:
        /* MSA */
        cout << "msa" << endl;
        if (method == "linear")
        {
            msa2leds_linear(ifs,ofs,l);
        }
        if(method == "cartesian"){
            msa2leds_cartesian(ifs,ofs,l);
        }
        
        break;
    case 2:
        /*  VCF */
        cout << "vcf" << endl;
        TODO();
        break;
    default:
        cout << "Error: Unknow extension!" << endl;
        return 1;
        break;
    }
    
    return 0;
}