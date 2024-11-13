#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <cassert>

#include "eds.hpp"
#include "index.hpp"

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

int main(int argc, char const *argv[]) {

    build_from_file();
    build_from_eds();
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}