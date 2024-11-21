#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <cassert>

#include "eds.hpp"
#include "index.hpp"
#include "utils.hpp"

using namespace bio_fmi;

/*  General EDS */
void build_from_file() {
    std::filesystem::path path = "/home/draesdom/Projects/bio-fmi/scripts/test/test.eds";
    Bio_FMi index = Bio_FMi(path,1);
    index.build();

    index.print_stats();
}


/*  no {} around single changes*/
void build_from_eds() {
    std::cout << " General {}" << std::endl;
    std::string eds_s = "{,ATGAGT,GT}{CT}{C,G}{TCT}{,CCTGCCG,GATAAGG}{ACAA}{AAGCAACA,GACCAACG,GGCTGCTG}{T}{C,G}{AA}{A,G}{GC}{CCTA,CGCC,TGCC}{TGG}{,GGGAAG,GGTAAG}{ATTGGT}{G,T}{GCCATG}";
    std::istringstream iss(eds_s);
    EDS eds = EDS(iss);

    // Bio_FMi index = Bio_FMi(eds,1);
}

void build_ref_files() {
    std::vector<unsigned> test_contexts = {3,5,10};
    std::filesystem::path basepath = "./test.eds";
    std::filesystem::path path;

    for (auto it: test_contexts)
    {
        path = basepath.replace_extension(std::to_string(it)+".leds");
        std::cout << path;
        if(std::filesystem::exists(path)){
            //  build the index
            Bio_FMi index = Bio_FMi(path,it);
            index.build();
            // index.print_stats();
        }else{
            std::ifstream ifs(basepath);
            std::ofstream ofs(path);
            eds2leds_cartesian(ifs,ofs,it);
            Bio_FMi index = Bio_FMi(path,it);
            index.build();
            // index.print_stats();
        }
    }
}


int main(int argc, char const *argv[]) {

    build_from_file();
    build_from_eds();
    build_ref_files();
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}