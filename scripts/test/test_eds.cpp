#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <cassert>

#include "eds.hpp"


/*  General EDS */
void create_eds1() {
    std::cout << " General {}" << std::endl;
    std::string eds_s = "{,ATGAGT,GT}{CT}{C,G}{TCT}{,CCTGCCG,GATAAGG}{ACAA}{AAGCAACA,GACCAACG,GGCTGCTG}{T}{C,G}{AA}{A,G}{GC}{CCTA,CGCC,TGCC}{TGG}{,GGGAAG,GGTAAG}{ATTGGT}{G,T}{GCCATG}";
    std::istringstream iss(eds_s);
    EDS eds = EDS(iss);

    eds.stats();

    assert(eds.n == eds.set_size_.size());
    assert(eds.set_size_.size() == 18);
    
    int N = 0;
    for (auto change : eds.changes)
        N += change.size(); 

    assert(eds.N == 107);
    assert(eds.N == N);
    assert(eds.m == 32);
    assert(eds.changes.size() == 32);
    assert(eds.max_l == 6);
    assert(eds.min_l == 1);
    assert(eds.n_empty_strings == 3);
}


/*  no {} around single changes*/
void create_eds2() {
    std::cout << " NO reference {}" << std::endl;
    std::string eds_s = "{,ATGAGT,GT}CT{C,G}TCT{,CCTGCCG,GATAAGG}ACAA{AAGCAACA,GACCAACG,GGCTGCTG}T{C,G}AA{A,G}GC{CCTA,CGCC,TGCC}TGG{,GGGAAG,GGTAAG}ATTGGT{G,T}GCCATG";
    std::istringstream iss(eds_s);
    EDS eds = EDS(iss);

    eds.stats();

    assert(eds.n == eds.set_size_.size());
    assert(eds.set_size_.size() == 18);
    
    int N = 0;
    for (auto change : eds.changes)
        N += change.size(); 

    assert(eds.N == 107);
    assert(eds.N == N);
    assert(eds.m == 32);
    assert(eds.changes.size() == 32);
    assert(eds.max_l == 6);
    assert(eds.min_l == 1);
    assert(eds.n_empty_strings == 3);
}

/*  no common parts*/
void create_eds3() {
    std::cout << " Degenerate symbols next to each other " << std::endl;
    std::string eds_s = "{,ATGAGT,GT}{CT,AA}{C,G}{TCT,G}{,CCTGCCG,GATAAGG}{,ACAA}{AAGCAACA,GACCAACG,GGCTGCTG}{A,T,G}{C,G}{AA}{A,G}{GC}{CCTA,CGCC,TGCC}{TGG}{,GGGAAG,GGTAAG}{A,ATTGGT}{G,T}{GCCATG,}";
    std::istringstream iss(eds_s);
    EDS eds = EDS(iss);
    eds.stats();

    assert(eds.n == eds.set_size_.size());
    assert(eds.set_size_.size() == 18);
    
    int N = 0;
    for (auto change : eds.changes)
        N += change.size(); 

    assert(eds.N == 113);
    assert(eds.N == N);
    assert(eds.m == 39);
    assert(eds.changes.size() == 39);
    assert(eds.max_l == 3);
    assert(eds.min_l == 2);
    assert(eds.n_empty_strings == 5);

}

/*  no common parts*/
void create_eds4() {
    std::cout << " Degenerate symbols next to each other, without {} for |DS| = 1 " << std::endl;
    std::string eds_s = "{,ATGAGT,GT}{CT,AA}{C,G}{TCT,G}{,CCTGCCG,GATAAGG}{,ACAA}{AAGCAACA,GACCAACG,GGCTGCTG}{A,T,G}{C,G}AA{A,G}GC{CCTA,CGCC,TGCC}TGG{,GGGAAG,GGTAAG}{A,ATTGGT}{G,T}{GCCATG,}";
    std::istringstream iss(eds_s);
    EDS eds = EDS(iss);
    eds.stats();

    assert(eds.n == eds.set_size_.size());
    assert(eds.set_size_.size() == 18);
    
    int N = 0;
    for (auto change : eds.changes)
        N += change.size(); 

    assert(eds.N == 113);
    assert(eds.N == N);
    assert(eds.m == 39);
    assert(eds.changes.size() == 39);
    assert(eds.max_l == 3);
    assert(eds.min_l == 2);
    assert(eds.n_empty_strings == 5);
}


int main(int argc, char const *argv[]) {
    create_eds1(); 
    create_eds2(); 
    create_eds3(); 
    create_eds4(); 

    std::cout << "All tests passed!" << std::endl;
    return 0;
}