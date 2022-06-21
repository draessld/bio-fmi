#ifndef CONFIG_H
#define CONFIG_H

#include <fstream>
#include <iostream>
#include <filesystem>

class Config
{
public:

    unsigned type;
    unsigned context_length = 6;
    std::filesystem::path index_input_path;
    std::filesystem::path index_output_path;
    std::filesystem::path pattern_file;
    bool silent;
    bool pattern_mode;
    unsigned repetition = 1;

    Config(/* args */);
    ~Config();
    int save();
    int load();
};

Config::Config(/* args */){}

Config::~Config(){}

int Config::save(){
    // std::cout << "saving to:"<<index_output_path.replace_extension() / "config.cfg" << std::endl; 
    std::ofstream file(index_output_path.replace_extension() / "config.cfg");
    if (!file.is_open())
    {
        std::cerr << "Error while reading file" << std::endl;
        return -1;
    }

    file << type << ' ' << context_length;
    file.close();
    return 0;
}

int Config::load(){
    std::ifstream file(index_input_path.replace_extension()/ "config.cfg");
    if (!file.is_open())
    {
        std::cerr << "Error while reading file" << std::endl;
        return -1;
    }
    
    file >> type >> context_length;

    return 0;
}

#endif // CONFIG_H
