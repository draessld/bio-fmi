#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <cassert>

#include "eds.hpp"
#include "index.hpp"
#include "utils.hpp"

using namespace bio_fmi;

void transform(std::filesystem::path path, unsigned context)
{
    std::filesystem::path out = path;
    out.replace_extension(std::to_string(context) + ".leds");
    if (!std::filesystem::exists(out))
    {
        //  transform the eds
        std::cout << "leds on " << out << " does not exists: transforming" << std::endl;
        std::ifstream ifs(path);
        std::ofstream ofs(out);
        eds2leds_cartesian(ifs, ofs, context);
        ifs.close();
        ofs.close();
    }
}

void build(std::filesystem::path path, unsigned context)
{
    std::filesystem::path out = path;
    out.replace_extension("leds.index");
    if (!std::filesystem::exists(out))
    {
        //  build the index
        std::cout << "index on " << out << " does not exists: building" << std::endl;
        Bio_FMi index = Bio_FMi(path, context);
        index.build();
    }
}

void load()
{
}

void locate_run(Bio_FMi &index, std::string p)
{
}

void locate_small1()
{
    //                     0             1              2
    //                     123456     78901234567     89012345678
    std::string eds_str = "TCACCG{A,T}GTGCATGTGCA{A,G}AGCATGACAGC";
    std::istringstream iss(eds_str);
    EDS eds(iss);
    // eds.stats();
    std::cout << "Parsing done" << std::endl;
    Bio_FMi index = Bio_FMi(eds, 5);
    index.build();
    // index.print();
    std::cout << "Building done" << std::endl;

    Bio_FMi::hash_type res;
    std::string p;
    {
        p = "TCACC"; //  1 []
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 1);
        assert(res.begin()->second.begin()->second.empty());
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "TGCATGTGCA"; //  8 []
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 8);
        assert(res.begin()->second.begin()->second.empty());
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "CCGAG"; //  4 [1]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 4);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 1);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "GTGTG"; //  6 [2]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 6);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 2);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "AGTGC"; //  7 [1]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 7);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 1);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "GAGCA"; //  18 [4]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 18);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 4);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "TGCAA"; //  14 [3]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 14);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 3);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "CCGAGTGCATGTGCAGAGCA"; //  4 [1,4]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 4);
        assert(res.begin()->second.begin()->second.size() == 2);
        assert(res.begin()->second.begin()->second[0] == 1);
        assert(res.begin()->second.begin()->second[1] == 4);
        std::cout << ">" << p << " -> OK" << std::endl;
    }
}

void locate_small2()
{

    //                          0        1
    //                          12345678901     234
    std::string eds_str = "{A,T}GTGCATGTGCA{A,G}AGC";
    std::istringstream iss(eds_str);
    EDS eds(iss);
    // eds.stats();
    std::cout << "Parsing done" << std::endl;
    Bio_FMi index = Bio_FMi(eds, 5);
    index.build();
    // index.print();
    std::cout << "Building done" << std::endl;

    Bio_FMi::hash_type res;
    std::string p;

    {
        p = "CATGT"; //  4 []
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 4);
        assert(res.begin()->second.begin()->second.empty());
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "GTGCATGTGC"; //  1 []
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 1);
        assert(res.begin()->second.begin()->second.empty());
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "AGTGC"; //  1 [1]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 1);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 1);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "AGAGC"; //  11 [4]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 11);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 4);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "TGCAA"; //  8 [3]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 8);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 3);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "TGTGCATGTGCAAAG"; //  1 [2,3]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 1);
        assert(res.begin()->second.begin()->second.size() == 2);
        assert(res.begin()->second.begin()->second[0] == 2);
        assert(res.begin()->second.begin()->second[1] == 3);
        std::cout << ">" << p << " -> OK" << std::endl;
    }
}

void locate_small3()
{
    //                     0             1             
    //                     123     45678901234     
    std::string eds_str = "CCG{A,T}GTGCATGTGCA{A,G}";
    std::istringstream iss(eds_str);
    EDS eds(iss);
    // eds.stats();
    std::cout << "Parsing done" << std::endl;
    Bio_FMi index = Bio_FMi(eds, 5);
    index.build();
    // index.print();
    std::cout << "Building done" << std::endl;

    Bio_FMi::hash_type res;
    std::string p;

    {
        p = "CATGT"; //  7 []
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 7);
        assert(res.begin()->second.begin()->second.empty());
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "GTGCATGTGC"; //  4 []
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 4);
        assert(res.begin()->second.begin()->second.empty());
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "CGTGT"; //  2 [2]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 2);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 2);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "TGCAA"; //  11 [3]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 11);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 3);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "TGCAA"; //  11 [3]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 11);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 3);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "CGAGTGCATGTGCAA"; //  2 [1,3]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->first == 2);
        assert(res.begin()->second.begin()->second.size() == 2);
        assert(res.begin()->second.begin()->second[0] == 1);
        assert(res.begin()->second.begin()->second[1] == 3);
        std::cout << ">" << p << " -> OK" << std::endl;
    }
}

void locate_small_nondeterministic()
{
    //                     0                 1             
    //                     123456         78901234     
    std::string eds_str = "TCACCG{AT,AA,T}GTGCATGTG";
    std::istringstream iss(eds_str);
    EDS eds(iss);
    // eds.stats();
    std::cout << "Parsing done" << std::endl;
    Bio_FMi index = Bio_FMi(eds, 5);
    index.build();
    // index.print();
    std::cout << "Building done" << std::endl;

    Bio_FMi::hash_type res;
    std::string p;

    {
        p = "ACCGA";    //  3[1], 3[2]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 2);
        assert(res.begin()->second.begin()->first == 3);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 2);

        assert((res.begin()->second.begin()+1)->first == 3);
        assert((res.begin()->second.begin()+1)->second.size() == 1);
        assert((res.begin()->second.begin()+1)->second[0] == 1);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

    {
        p = "TGTGC";    //  8[1], 7[3]
        index.locate(p);
        res = index.get_result();
        assert(res.begin()->second.size() == 2);
        assert(res.begin()->second.begin()->first == 7);
        assert(res.begin()->second.begin()->second.size() == 1);
        assert(res.begin()->second.begin()->second[0] == 3);

        assert((res.begin()->second.begin()+1)->first == 8);
        assert((res.begin()->second.begin()+1)->second.size() == 1);
        assert((res.begin()->second.begin()+1)->second[0] == 1);
        std::cout << ">" << p << " -> OK" << std::endl;
    }

}

/*  General EDS */
void locate3()
{
    std::filesystem::path test_file = "/home/draesdom/Projects/bio-fmi/scripts/test/test.eds";
    std::filesystem::path file_path = "/home/draesdom/Projects/bio-fmi/scripts/test/test.3.leds";
    std::filesystem::path index_path = "/home/draesdom/Projects/bio-fmi/scripts/test/test.3.leds.index";

    transform(test_file, 3);
    build(file_path, 3);

    std::string p = "TGATAAGGACAAGACCAACGTCAAGGCCGCCTGGATTGGTGG";
    Bio_FMi index = Bio_FMi(index_path);
    index.context_length_ = 3;
    index.print_stats();
    index.locate(p);

    Bio_FMi::hash_type results = index.get_result();
    for (const auto &pair : results)
    {
        for (auto occ : pair.second)
        {
            std::cout << occ.first << "[ ";
            for (int num : occ.second)
            {
                std::cout << num << " ";
            }
            std::cout << "]" << std::endl;
        }
    }
}

int main(int argc, char const *argv[])
{
    // locate3();
    locate_small1();
    locate_small2();
    locate_small3();
    locate_small_nondeterministic();
    return 0;
}