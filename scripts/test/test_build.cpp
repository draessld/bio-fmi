#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <cassert>

#include "EDS.hpp"
#include "index.hpp"
#include "utils.hpp"

using namespace bio_fmi;

/*  General EDS */
void build_from_file() {
    std::cout << " Build from file" << std::endl;
    std::filesystem::path path = "../../test/test.eds";
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
    std::filesystem::path basepath = "../../test.eds";
    std::filesystem::path path;

    for (auto it: test_contexts)
    {
        path = basepath.replace_extension(std::to_string(it)+".leds");
        std::cout << " Build with context: "<<it << " for " << path << std::endl;
        if(std::filesystem::exists(path)){
            std::cout << "exists" << std::endl;
            
            //  build the index
            Bio_FMi index = Bio_FMi(path,it);
            index.build();
            // index.print_stats();
        }else{
            std::cout << "transforming" << std::endl;
            std::ifstream ifs(basepath);
            std::ofstream ofs(path);
            if (!ifs)  // If the file could not be opened, treat input as a sequence
            {
                std::cout << "Error: File was not found!" << std::endl;
            }
            eds2leds_cartesian(ifs,ofs,it);
            ifs.close();
            ofs.close();
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