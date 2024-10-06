#include <filesystem>
#include <string>
#include <iostream>

#include "utils/utils.hpp"
#include "utils/index.h"

#include <boost/program_options.hpp>

using namespace bio_fmi;
namespace po = boost::program_options;

#define BENCHMARK_REP 100

std::string usage = "";
std::string desc = "";
std::string pattern;
std::filesystem::path in_file;
std::filesystem::path out_file;
std::filesystem::path pattern_file;
bool benchmark = false;

/*  parse parameters   */
int handle_parameters(int argc, const char **argv)
{
    po::options_description desc("Parameters");
    desc.add_options()("help,h", "display help message")("input,i", po::value<std::filesystem::path>(&in_file)->required(), "index file path")("output,o", po::value<std::filesystem::path>(&out_file), "")("pattern,p", po::value<std::string>(&pattern), "print occurences of every pattern")("pattern-file,P", po::value<std::filesystem::path>(&pattern_file), "input pattern file path (positional arg 2)")("benchmark", "run benchmark testing");

    po::positional_options_description posOptions;

    posOptions.add("index-path", 1);

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

        if (vm.count("benchmark"))
        {
            benchmark = true;
        }

        if (!vm.count("pattern-file") && !vm.count("pattern"))
        {
            std::cerr << "Usage: " << argv[0] << " " << usage << std::endl
                      << std::endl;
            std::cerr << desc << std::endl;

            std::cerr << "Error: " << "Missing pattern" << std::endl;
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

/*  create index and run experiments   */
void run()
{
    // long mem_baseline;
    // auto time_baseline = std::chrono::high_resolution_clock::now();
    // std::vector<double> times;
    // std::vector<std::string> patterns;
    // // std::filesystem::path eds_file = cfg.index_bed;
    // // eds_file.replace_extension("");

    // // std::string eds = readEDSFromFile(eds_file);

    // Bio_FMi *index = new Bio_FMi(cfg.index_bed);

    // /*  read patterns  */
    // if (!cfg.pattern.empty())
    //     patterns.push_back(cfg.pattern);

    // if (!cfg.pattern_file.empty())
    //     if (read_patterns(cfg.pattern_file, patterns))
    //         std::cerr << "Unsucessfull pattern reading" << std::endl;

    // mem_baseline = get_mem_usage();

    // double total_time = 0;
    // int total_occurences = 0;
    // int p_occurences = 0;

    // std::cout << "Number of patterns n = " << patterns.size() << std::endl;

    // /*  search patterns  */
    // for (auto P : patterns)
    // {
    //     p_occurences = 0;

    //     time_baseline = std::chrono::high_resolution_clock::now();
    //     index->locate(P);
    //     auto time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - time_baseline);
    //     times.push_back(time.count());

    //     auto result = index->get_result();
    //     for (auto &pos : result)
    //     {
    //         p_occurences += pos.second.size();
    //         // for (auto &res: pos.second)
    //         // {
    //         // std::cout << "hop" << std::endl;
    //         // if(check_result(eds,res.first,P.size(),res.second) != P)
    //         //     std::cout << "Invalid result starting of pattern: " << P << " starting on position: " << res.first << " with changes: " << res.second << std::endl;            }
    //     }
    //     std::cout << ">" << P << '\t' << p_occurences << '\t' << times.back() << std::endl;
    //     index->print_result(result);
    //     total_occurences += p_occurences;
    // }
    // total_time = std::reduce(times.begin(), times.end());
    // std::cout << "Peak RAM usage: " << double(get_mem_usage() - mem_baseline) / double(patterns.size()) << " kB" << std::endl;
    // std::cout << "Average number of occurrences per pattern: " << total_occurences / patterns.size() << std::endl;
    // std::cout << "Total number of occurrences: " << total_occurences << std::endl;
    // std::cout << "Total time: " << total_time << " microseconds" << std::endl;
    // std::cout << "Average time per pattern: " << total_time / patterns.size() << " microseconds" << std::endl;
    // std::cout << "Average time per occurence: " << total_time / total_occurences << " microseconds" << std::endl;

    // delete index;
}

void run_benchmark()
{
    // long mem_baseline;
    // auto time_baseline = std::chrono::high_resolution_clock::now();
    // std::vector<double> times;
    // std::vector<std::string> patterns;
    // std::filesystem::path eds_file = cfg.index_bed;
    // eds_file.replace_extension("");

    // std::string eds = readEDSFromFile(eds_file);

    // Bio_FMi *index = new Bio_FMi(cfg.index_bed);

    // // for(size_t i = index->context_length_; i< (index->context_length_*5); i+=index->context_length_)
    // //     generate_patterns(eds,cfg.repetition,i,patterns);

    // generate_patterns(eds, cfg.repetition, index->context_length_ * 2, patterns);
    // generate_patterns(eds, cfg.repetition, index->context_length_ * 3, patterns);
    // generate_patterns(eds, cfg.repetition, index->context_length_ * 4, patterns);

    // mem_baseline = get_mem_usage();

    // double total_time = 0;
    // int total_occurences = 0;
    // int p_occurences = 0;

    // std::cout << "Number of patterns n = " << patterns.size() << std::endl;

    // /*  search patterns  */
    // for (auto P : patterns)
    // {
    //     p_occurences = 0;

    //     time_baseline = std::chrono::high_resolution_clock::now();
    //     index->locate(P);
    //     auto time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - time_baseline);
    //     times.push_back(time.count());

    //     auto result = index->get_result();
    //     for (auto &pos : result)
    //     {
    //         p_occurences += pos.second.size();
    //         for (auto &res : pos.second)
    //         {
    //             if (check_result(eds, res.first, P.size(), res.second) != P)
    //                 std::cout << "Invalid result starting of pattern: " << P << " starting on position: " << res.first << " with changes: " << res.second << std::endl;
    //         }
    //     }
    //     std::cout << ">" << P << '\t' << p_occurences << '\t' << times.back() << std::endl;
    //     index->print_result(result);
    //     total_occurences += p_occurences;
    // }
    // total_time = std::reduce(times.begin(), times.end());
    // std::cout << "Peak RAM usage: " << double(get_mem_usage() - mem_baseline) / double(patterns.size()) << " kB" << std::endl;
    // std::cout << "Average number of occurrences per pattern: " << total_occurences / patterns.size() << std::endl;
    // std::cout << "Total number of occurrences: " << total_occurences << std::endl;
    // std::cout << "Total time: " << total_time << " microseconds" << std::endl;
    // std::cout << "Average time per pattern: " << total_time / patterns.size() << " microseconds" << std::endl;
    // std::cout << "Average time per occurence: " << total_time / total_occurences << " microseconds" << std::endl;

    // delete index;
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

    // AAAACGTGTCATTAC
    if (benchmark)
        run_benchmark();
    else
        run();

    return EXIT_SUCCESS;
}
