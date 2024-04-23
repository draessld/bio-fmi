#include <filesystem>
#include <string>
#include <iostream>
#include <sstream>

#include "utils/Informations_signs.h"
#include "index/Bio_FMi.h"
#include "utils/utils.hpp"

#include <boost/program_options.hpp>

using namespace bio_fmi;
namespace po = boost::program_options;

struct Config
{
    int l = 6;
    std::filesystem::path msa_file;
    std::filesystem::path eds_file;
    std::filesystem::path index_output_folder;
    std::filesystem::path pattern_file;
    bool implicit;
};

Informations_signs inf;
Config cfg;

int handle_parameters(int argc, const char **argv)
{
    po::options_description desc("Allowed options");

    desc.add_options()("help", "produce help message")
    ("help-verbose", "display verbose help message")
    ("version,v", "display version info")
    ("context_length,l", po::value<int>(&cfg.l)->required(), "length of chunk and stored context, default 5")
    ("msa-file,m", po::value<std::filesystem::path>(&cfg.msa_file), "")
    ("eds-file,e", po::value<std::filesystem::path>(&cfg.eds_file), "");

    po::positional_options_description posOptions;

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(posOptions).run(), vm);

        if (vm.count("help"))
        {
            std::cout << "Usage: " << argv[0] << " " << inf.buildUsageInfoString << std::endl
                      << std::endl;
            std::cout << desc << std::endl;

            return 1;
        }
        if (vm.count("help-verbose"))
        {
            std::cout << inf.verboseInfoString << std::endl
                      << std::endl;

            std::cout << "Usage: " << argv[0] << " " << inf.buildUsageInfoString << std::endl
                      << std::endl;

            std::cout << desc << std::endl;
            std::cout << inf.verboseparametersStringBuild << std::endl;

            return 1;
        }
        if (vm.count("version"))
        {
            std::cout << inf.versionInfo << std::endl;
            return 1;
        }

        if ((vm.count("eds-file") == 0) && (vm.count("msa-file") == 0))
        {
            std::cout << "Error: No input given" << std::endl;

            std::cout << "Usage: " << argv[0] << " " << inf.buildUsageInfoString << std::endl
                      << std::endl;
            return -1;
        }
       
        if (vm.count("msa-file") != 0)
        {
            cfg.implicit = true;
        }
        if (vm.count("eds-file") != 0)
        {
            cfg.implicit = false;
        }

        po::notify(vm);
    }
    catch (const po::error &e)
    {
        std::cerr << "Usage: " << argv[0] << " " << inf.buildUsageInfoString << std::endl
                  << std::endl;
        std::cerr << desc << std::endl;

        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}

int run_explicit()   // from explicitly given EDS create an index
{
    //  transform MSA to EDS based on given context length
    auto base = std::chrono::high_resolution_clock::now();
    auto startTime = base;

    // get_l_eds(cfg.msa_file, cfg.eds_file, cfg.l);
    
    std::cout << "MSA2EDS time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() << " um" << std::endl;

    long mem_baseline;
    Bio_FMi *index = new Bio_FMi(cfg.eds_file,cfg.l);

    mem_baseline = get_mem_usage();
    base = std::chrono::high_resolution_clock::now();
    startTime = base;

    if(index->build()){
        std::cout << "Building incomplete" << std::endl;
        return -1;
    }

    std::cout << "Build time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() << " um" << std::endl;
    std::cout << "Index size: " << index->total_index_size_ << " MB" << std::endl;
    std::cout << "Peak RAM usage: " << get_mem_usage() - mem_baseline << " kB" << std::endl;

    delete index;
    // std::cout  << std::endl;
    // std::cout << "Statistics about EDS: " << std::endl;
    // index->print_stats();
    return 0;
}

int run_implicit(){ // from given MSA create an index without creating explicit EDS file

    //  transform MSA to EDS based on given context length
    auto base = std::chrono::high_resolution_clock::now();
    auto startTime = base;

    // get_l_eds(cfg.msa_file, cfg.eds_file, cfg.l);   //  TODO 
    
    std::cout << "MSA2EDS time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() << " um" << std::endl;

    long mem_baseline;
    Bio_FMi *index = new Bio_FMi(cfg.msa_file,cfg.l);

    mem_baseline = get_mem_usage();
    base = std::chrono::high_resolution_clock::now();
    startTime = base;

    if(index->build()){
        std::cout << "Building incomplete" << std::endl;
        return -1;
    }

    std::cout << "Build time: " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - startTime).count() << " um" << std::endl;
    std::cout << "Index size: " << index->total_index_size_ << " MB" << std::endl;
    std::cout << "Peak RAM usage: " << get_mem_usage() - mem_baseline << " kB" << std::endl;

    delete index;
    // std::cout  << std::endl;
    // std::cout << "Statistics about EDS: " << std::endl;
    // index->print_stats();
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

    //  validate input data
    if (cfg.implicit)
    {
        //  run index
        run_implicit();
    
    }else{
        // check validity of context length and the shortest common part in EDS
        std::cout << "  Validating EDS ...";
        if(validate_eds(cfg.eds_file,cfg.l)!=0){

            std::cout << "Unable to work with context length shorter than is some common part in EDS\n" << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << "  done!" << std::endl;

        //  run index
        run_explicit();
    }


    return 0;
}
