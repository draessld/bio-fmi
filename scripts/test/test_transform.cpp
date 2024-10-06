#include <iostream>
#include <fstream>
#include <filesystem>
#include <cassert>

#include "utils/utils.hpp"


void create_test_file(const std::filesystem::path& path, const std::string& content) {
    std::ofstream out(path);
    if (!out.is_open()) {
        std::cerr << "Error: Could not create test file " << path << std::endl;
        exit(1);
    }
    out << content;
    out.close();
}

void test_msa2leds_linear() {
    std::string msa = ">test1\nATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG------ATTGGTGGC\nCATG\n>test2\n------CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGATTGGTTGC\nCATG\n>test3\n----GTCTCTCT-------ACAAAAGCAACATCAAGGCTGCCTGGGGGAAGATTGGTGGC\nCATG";
    std::string leds0 = "{,ATGAGT,GT}{CT}{C,G}{TCT}{,CCTGCCG,GATAAGG}{ACAA}{AAGCAACA,GACCAACG,GGCTGCTG}{T}{C,G}{AA}{A,G}{GC}{CCTA,CGCC,TGCC}{TGG}{,GGGAAG,GGTAAG}{ATTGGT}{G,T}{GCCATG}";
    std::string leds3 = "{ATGAGTCTC,CTG,GTCTC}{TCT}{,CCTGCCG,GATAAGG}{ACAA}{AAGCAACATCAAGGCTGCC,GACCAACGTCAAGGCCGCC,GGCTGCTGTGAAAGCCCTA}{TGG}{,GGGAAG,GGTAAG}{ATTGGT}{G,T}{GCCATG}";
    std::string leds5 = "{ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG,CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG,GTCTCTCTACAAAAGCAACATCAAGGCTGCCTGGGGGAAG}{ATTGGT}{G,T}{GCCATG}";
    std::istringstream iss(msa);
    std::ostringstream oss;

    std::cout << msa <<std::endl;

    msa2leds_linear(iss,oss,0);
    std::cout << 0 << ": " << oss.str() <<std::endl;
    assert(leds0 == oss.str());

    iss.seekg(0);
    oss.str("");
    oss.clear();
    msa2leds_linear(iss,oss,3);
    std::cout << 3 << ": " << oss.str() <<std::endl;
    assert(leds3 == oss.str());


    iss.seekg(0);
    oss.str("");
    oss.clear();
    msa2leds_linear(iss,oss,5);
    std::cout << 5 << ": " << oss.str() <<std::endl;
    assert(leds5 == oss.str());
}

void test_msa2leds_cartesian() {
    std::string msa = ">test1\nATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG------ATTGGTGGC\nCATG\n>test2\n------CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGATTGGTTGC\nCATG\n>test3\n----GTCTCTCT-------ACAAAAGCAACATCAAGGCTGCCTGGGGGAAGATTGGTGGC\nCATG";
    std::string leds0 = "{,ATGAGT,GT}{CT}{C,G}{TCT}{,CCTGCCG,GATAAGG}{ACAA}{AAGCAACA,GACCAACG,GGCTGCTG}{T}{C,G}{AA}{A,G}{GC}{CCTA,CGCC,TGCC}{TGG}{,GGGAAG,GGTAAG}{ATTGGT}{G,T}{GCCATG}";
    std::string leds3 = "{ATGAGTCTC,ATGAGTCTG,CTC,CTG,GTCTC,GTCTG}{TCT}{,CCTGCCG,GATAAGG}{ACAA}{AAGCAACATCAAAGCCCTA,AAGCAACATCAAAGCCGCC,AAGCAACATCAAAGCTGCC,AAGCAACATCAAGGCCCTA,AAGCAACATCAAGGCCGCC,AAGCAACATCAAGGCTGCC,AAGCAACATGAAAGCCCTA,AAGCAACATGAAAGCCGCC,AAGCAACATGAAAGCTGCC,AAGCAACATGAAGGCCCTA,AAGCAACATGAAGGCCGCC,AAGCAACATGAAGGCTGCC,GACCAACGTCAAAGCCCTA,GACCAACGTCAAAGCCGCC,GACCAACGTCAAAGCTGCC,GACCAACGTCAAGGCCCTA,GACCAACGTCAAGGCCGCC,GACCAACGTCAAGGCTGCC,GACCAACGTGAAAGCCCTA,GACCAACGTGAAAGCCGCC,GACCAACGTGAAAGCTGCC,GACCAACGTGAAGGCCCTA,GACCAACGTGAAGGCCGCC,GACCAACGTGAAGGCTGCC,GGCTGCTGTCAAAGCCCTA,GGCTGCTGTCAAAGCCGCC,GGCTGCTGTCAAAGCTGCC,GGCTGCTGTCAAGGCCCTA,GGCTGCTGTCAAGGCCGCC,GGCTGCTGTCAAGGCTGCC,GGCTGCTGTGAAAGCCCTA,GGCTGCTGTGAAAGCCGCC,GGCTGCTGTGAAAGCTGCC,GGCTGCTGTGAAGGCCCTA,GGCTGCTGTGAAGGCCGCC,GGCTGCTGTGAAGGCTGCC}{TGG}{,GGGAAG,GGTAAG}{ATTGGT}{G,T}{GCCATG}";
    std::string leds5 = "{ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG,CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG,GTCTCTCTACAAAAGCAACATCAAGGCTGCCTGGGGGAAG}{ATTGGT}{G,T}{GCCATG}";
    std::istringstream iss(msa);
    std::ostringstream oss;

    std::cout << msa <<std::endl;

    msa2leds_cartesian(iss,oss,0);
    std::cout << 0 << ": " << oss.str() <<std::endl;
    assert(leds0 == oss.str());

    iss.seekg(0);
    oss.str("");
    oss.clear();
    msa2leds_cartesian(iss,oss,3);
    std::cout << 3 << ": " << oss.str() <<std::endl;
    assert(leds3 == oss.str());


    iss.seekg(0);
    oss.str("");
    oss.clear();
    msa2leds_cartesian(iss,oss,5);
    std::cout << 5 << ": " << oss.str() <<std::endl;
    assert(leds5 == oss.str());
}

int main(int argc, char const *argv[]) {
    test_msa2leds_linear();  // Test with valid data
    test_msa2leds_cartesian();  // Test with empty blocks

    std::cout << "All tests passed!" << std::endl;
    return 0;
}