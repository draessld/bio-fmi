#include <filesystem>
#include <string>
#include <iostream>

#include "utils/Informations_signs.h"
#include "utils/Config.h"
#include "index/aln.h"
#include "index/eds.h"

#include <boost/program_options.hpp>

/*  
 *  This application can be used for pattern matching in aligned set of sequence or in eds format
 */

using namespace bio_fmi;
namespace po = boost::program_options;

Informations_signs inf;
Config cfg;

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

/*  parse parameters   */
int handle_parameters(int argc, const char **argv)
{
    po::options_description desc("Parameters");
    desc.add_options()("help,h", "display help message")
    ("help-verbose", "display verbose help message")
    ("version,v", "display version info")
    ("silent,s", "silent mode")
    ("pattern,p", "print occurences of every pattern")
    ("index-path,i", po::value<std::filesystem::path>(&cfg.index_input_path)->required(), "input text file path (positional arg 1)")
    ("pattern-file,I", po::value<std::filesystem::path>(&cfg.pattern_file)->required(), "input pattern file path (positional arg 2)");

    po::positional_options_description posOptions;

    posOptions.add("index-path", 1);
    posOptions.add("pattern-file", 1);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(posOptions).run(), vm);
        if (vm.count("silent"))
        {
            cfg.silent = true;
        }
        if (vm.count("pattern"))
        {
            cfg.pattern_mode = true;
        }
        if (vm.count("help"))
        {
            std::cout << "Usage: " << argv[0] << " " << inf.locateUsageInfoString << std::endl
                      << std::endl;
            std::cout << desc << std::endl;

            return 1;
        }
        if (vm.count("help-verbose"))
        {
            std::cout << inf.verboseInfoString << std::endl
                      << std::endl;
            std::cout << "Usage: " << argv[0] << " " << inf.locateUsageInfoString << std::endl
                      << std::endl;

            std::cout << desc << std::endl;
            std::cout << inf.verboseparametersStringSearch << std::endl;

            return 1;
        }
        if (vm.count("version"))
        {
            std::cout << inf.versionInfo << std::endl;
            return 1;
        }

        po::notify(vm);
    }
    catch (const po::error &e)
    {
        std::cerr << "Usage: " << argv[0] << " " << inf.locateUsageInfoString << std::endl
                  << std::endl;
        std::cerr << desc << std::endl;

        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    return 0;
}


/*  create index and run experiments   */
template<typename T>
void run(){
    double time_baseline;
    long mem_baseline;
    std::vector<double> times;
    std::vector<std::string> patterns;

    T *index = new T(cfg.context_length);
    
    /*  load index  */
    index->load(cfg.index_input_path);
    if (!cfg.silent)
    {
        index->print();
    }
    
    time_baseline = get_time_usage();

    /*  read patterns  */
    index->read_patterns(cfg.pattern_file, patterns);
    std::cout << "Patterns reading time: " << (get_time_usage()-time_baseline) << std::endl;
    std::cout << std::endl;

    mem_baseline = get_mem_usage();
    time_baseline = get_time_usage();

    double total_time = 0;
    int total_occurences = 0;
    int count = 0;

    std::cout << "Number of patterns n = " << patterns.size() << std::endl;
    std::cout << "Pattern length = " << patterns[0].size() << std::endl;
    
    /*  search patterns  */
    for (auto p : patterns)
    {
        count = index->search(p, cfg.silent);
        if(cfg.pattern_mode){
            std::cout << "Occurences for pattern " << p << ": " << count << std::endl;
        }
        total_occurences += count;
        total_time = get_time_usage()-time_baseline;
    }
    std::cout << "Peak RAM usage: " << double(get_mem_usage()-mem_baseline) / double(patterns.size()) << " kB" << std::endl;
    std::cout << total_occurences / patterns.size() << " average occurrences per pattern" << std::endl;
    std::cout << "Total number of occurrences  occ_t = " << total_occurences << std::endl;
    std::cout << "Total time : " << total_time << " seconds" << std::endl;
    std::cout << "Search time : " << total_time / patterns.size() << " seconds/pattern (total " << patterns.size() << " patterns)" << std::endl;
    std::cout << "Search time : " << total_time / total_occurences << " seconds/occurence (total " << total_occurences << " occurences)" << std::endl;

    delete index;
}


int main(int argc, char const *argv[])
{
    /*  Parse input text    */
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

    /*  load configuration  */
    cfg.load();

    /*  run  */
    if (cfg.type)
        run<eds>();
    else
        run<aln>();

    return EXIT_SUCCESS;
}
