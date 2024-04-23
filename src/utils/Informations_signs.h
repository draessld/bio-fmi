#ifndef INFORMATIONS_SIGNS_H
#define INFORMATIONS_SIGNS_H

#include <string>
#include <filesystem>
#include <iostream>

struct Informations_signs
{
    const std::string versionInfo = "bio-fmi v1.1";
    const std::string createEDSInfoString = "bio-fmi-createEDS [-c|-e] <msa_file>";
    const std::string buildUsageInfoString = "bio-fmi-build [-m|-e]<input_text_file>";
    const std::string locateUsageInfoString = "bio-fmi-locate <index_basename> -P<pattern_file> -p<pattern>";
    
    const std::string verboseInfoString =
        "";

    const std::string verboseparametersStringBuild =
        "";

    const std::string verboseparametersStringCreateEDS =
        "";

    const std::string verboseparametersStringSearch =
        "Input text file (positional parameter 1 or named parameter -i or --in-text-file) should contain the file with extension .aln in the format of sequence -AC--GT-CGTA and every sequence in new line or .eds in the format {A,C,}GAAT{,A,AT}ATT.\n"
        "Input pattern file (positional parameter 2 or named parameter -I or --in-pattern-file) should contain the information about number an length of patterns followed by one line of concatenated patterns.\n";

};

#endif // INFORMATIONS_SIGNS_H