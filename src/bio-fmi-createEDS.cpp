#include <filesystem>
#include <string>
#include <iostream>
#include <sstream>

#include "utils/Informations_signs.h"
#include "utils/utils.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

struct Config
{
    int l = 6;
    std::filesystem::path file;
    std::filesystem::path eds_file;
    int method;  //  type of used method, 1=leds, 2=cart-leds
};

Informations_signs inf;
Config cfg;

int handle_parameters(int argc, const char **argv)
{
    po::options_description desc("Allowed options");

    desc.add_options()("help", "produce help message")
    ("help-verbose", "display verbose help message")
    ("version,v", "display version info")
    ("cart-l-eds,c", "run benchmark test for given file - create EDS with l-eds and with cart-l-eds method and compare printout results")
    ("l-eds,e", "run benchmark test for given file - create EDS with l-eds and with cart-l-eds method and compare printout results")
    ("stats,s", "calculate statistics about given eds")
    ("context_length,l", po::value<int>(&cfg.l)->required(), "length of chunk and stored context, default 5")
    ("file", po::value<std::filesystem::path>(&cfg.file), "");

    po::positional_options_description posOptions;

    posOptions.add("file", 1);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(posOptions).run(), vm);

        if (vm.count("help"))
        {
            std::cout << "Usage: " << argv[0] << " " << inf.createEDSInfoString << std::endl
                      << std::endl;
            std::cout << desc << std::endl;

            return 1;
        }
        if (vm.count("help-verbose"))
        {
            std::cout << inf.verboseInfoString << std::endl
                      << std::endl;

            std::cout << "Usage: " << argv[0] << " " << inf.createEDSInfoString << std::endl
                      << std::endl;

            std::cout << desc << std::endl;
            std::cout << inf.verboseparametersStringCreateEDS << std::endl;

            return 1;
        }
        if (vm.count("version"))
        {
            std::cout << inf.versionInfo << std::endl;
            return 1;
        }
        if (vm.count("cart-l-eds") == 0 && vm.count("l-eds") == 0){
            std::cout << "Usage: " << argv[0] << " " << inf.createEDSInfoString << std::endl;
        }
        if(vm.count("cart-l-eds")){
            cfg.method = 2;   
        }

        if(vm.count("l-eds")){
            cfg.method = 1;   
        }

        if(vm.count("stats")){
            cfg.method = 0;   
        }
        po::notify(vm);
    }
    catch (const po::error &e)
    {
        std::cerr << "Usage: " << argv[0] << " " << inf.createEDSInfoString << std::endl
                  << std::endl;
        std::cerr << desc << std::endl;

        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

int run(){
    auto base = std::chrono::high_resolution_clock::now();
    auto startTime = base;
    switch (cfg.method)
    {
    case 0:
        /* EDS under l-context condition method */
        std::cout << "Counting statistic on given EDS file" << std::endl;
        validate_eds(cfg.file, cfg.l);
        break;
    case 1:
        /* EDS under l-context condition method */
        std::cout << "EDS under l-context condition method" << std::endl;
        get_l_eds(cfg.file, cfg.eds_file, cfg.l);
        std::cout << "EDS saved on " << cfg.eds_file << std::endl;
        break;
    case 2:
        /* cartesian product under l-context condition method */
        std::cout << "EDS with cartesian product under l-context condition method" << std::endl;
        get_cartl_eds(cfg.file, cfg.eds_file, cfg.l);
        std::cout << "EDS saved on " << cfg.eds_file << std::endl;
        break;
    
    default:
        std::cout << "Unknown method, exiting ...";
        return -1;
        break;
    }
    std::cout << "MSA2EDS time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() << " um" << std::endl;
    return 0;
}

int main(int argc, char const *argv[])
{
    // ru_maxrss
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

    run();

    return 0;
}
