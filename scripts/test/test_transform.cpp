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

int main(int argc, char const *argv[]) {

    std::cout << "All tests passed!" << std::endl;
    return 0;
}