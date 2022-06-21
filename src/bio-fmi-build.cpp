#include <filesystem>
#include <string>
#include <iostream>

#include "utils/Informations_signs.h"
#include "utils/Config.h"
#include "index/aln.h"
#include "index/eds.h"

#include <boost/program_options.hpp>

/*
 *
 *
 *
 *
 *
 */

using namespace bio_fmi;
namespace po = boost::program_options;

Informations_signs inf;
Config cfg;

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

int handle_parameters(int argc, const char **argv)
{
    po::options_description desc("Allowed options");

    desc.add_options()("help", "produce help message")("help-verbose", "display verbose help message")("version,v", "display version info")("silent,s", "silent mode")("repetition,r", po::value<unsigned>(&cfg.repetition), "number of repetition - for experiment needs")("context_length,l", po::value<unsigned>(&cfg.context_length), "length of chunk and stored context, default 5")("basefolder,o", po::value<std::filesystem::path>(&cfg.index_output_path), "use <basefolder> as prefix for all index files. Default: current folder is the specified input_file_name")("input-file,i", po::value<std::filesystem::path>(&cfg.index_input_path), "input file");

    po::positional_options_description posOptions;
    posOptions.add("input-file", 1);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(posOptions).run(), vm);

        if (vm.count("silent"))
        {
            cfg.silent = true;
        }
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
        if (vm.count("input-file") == 0)
        {
            std::cout << "Usage: " << argv[0] << " " << inf.buildUsageInfoString << std::endl
                      << std::endl;
            std::cout << desc << std::endl;

            return -1;
        }
        po::notify(vm);

        if (vm.count("basename") == 0)
        {
            cfg.index_output_path = cfg.index_input_path;
        }
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

template <typename T>
void run()
{

    double time_baseline;
    long mem_baseline;
    std::vector<double> times;
    std::vector<long> mems;

    T *index = new T(cfg.context_length);
    index->read(cfg.index_input_path);

    for (int i = 0; i < cfg.repetition; i++)
    {
        mem_baseline = get_mem_usage();
        time_baseline = get_time_usage();

        index->build();

        times.push_back((get_time_usage() - time_baseline));
        mems.push_back((get_mem_usage() - time_baseline));
    }

    std::cout << "Build time: " << std::reduce(times.begin(), times.end()) / cfg.repetition << std::endl;
    std::cout << "Peak RAM usage: " << std::reduce(mems.begin(), mems.end()) / cfg.repetition << " kB" << std::endl;

    if (!cfg.silent)
    {
        index->print();
    }
    index->save(cfg.index_output_path);
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

    if (cfg.index_input_path.extension() == ".aln")
    {
        cfg.type = inf.ALN;
        run<aln>();
    }
    else if (cfg.index_input_path.extension() == ".eds")
    {
        cfg.type = inf.EDS;
        run<eds>();
    }
    else
    {
        std::cout << "Unknown input format - use extension .aln or .eds" << std::endl;
        return EXIT_FAILURE;
    }

    /*  save configuration  */
    cfg.save();

    return 0;
}
