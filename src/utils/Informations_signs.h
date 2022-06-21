#ifndef INFORMATIONS_SIGNS_H
#define INFORMATIONS_SIGNS_H

#include <string>
#include <filesystem>
#include <iostream>

struct Informations_signs
{
    const std::string versionInfo = "bio-fmi v1.0";
    const std::string buildUsageInfoString = "[options] <input_file_name>";
    const std::string locateUsageInfoString = "[options] <index_basename> <pattern_file>";
    const std::string usageInfoString = "[options] <index_basename> <pattern_file>";
    
    const std::string verboseInfoString =
        "This software is called BIO-FMI (Index for set of genomes and pangenomes). It can be used for pattern matching in elastic-degenerate (ED text). Authors: Petr Prochazka, Jan Holub.\n"
        "Input format can be saved in .aln (ALN) file - alignment of sequences, every sequence in new line; or in .eds file - {A,C,}GAAT{,A,AT}ATT, where the strings in bracket are ordered ascending."
        "BIO-FMI return start positions of pattern occurences.";

    const std::string verboseparametersStringBuild =
        "Input text file (positional parameter 1 or named parameter -i or --in-text-file) should contain the file with extension .aln in the format of sequence -AC--GT-CGTA and every sequence in new line or .eds in the format {A,C,}GAAT{,A,AT}ATT.\n"
        "Context length (optional named parameter -l or --context_length with default value of 6) divide input pattern into chunks of this length and during the construction save same length right and left context of text input.";

    const std::string verboseparametersStringSearch =
        "Input text file (positional parameter 1 or named parameter -i or --in-text-file) should contain the file with extension .aln in the format of sequence -AC--GT-CGTA and every sequence in new line or .eds in the format {A,C,}GAAT{,A,AT}ATT.\n"
        "Input pattern file (positional parameter 2 or named parameter -I or --in-pattern-file) should contain the information about number an length of patterns followed by one line of concatenated patterns.\n";

    void build_information_sign(std::filesystem::path input, std::filesystem::path output)
    {
        std::cout << "Building bio-fm index of input_file " << input << std::endl;
        std::cout << "Index will be saved to " << output << std::endl;
    }
    const std::string search_information_sign = "Searching pattern";

    const unsigned ALN = 0;
    const unsigned EDS = 1;
};

#endif // INFORMATIONS_SIGNS_H