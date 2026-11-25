// #include <filesystem>
// #include <fstream>
// #include <iostream>

// class SourceManager{
//     private:
//     public:
//     std::ifstream ss;
//     std::filesystem::path tmpfile;
//     const int nsamples;
//     SourceManager(std::filesystem::path sourceFile){
//         ss.open(sourceFile);
//         if (!ss.is_open())
//             std::cerr << "Error: cannot open source file!" << std::endl; 
    
    
//         //  TODO: read source number = first int
//     }
    
//         ~SourceManager() {
//             if (ss.is_open())
//                 ss.close();
//         }
    
//         int transfer(std::istream from, std::ostream to, int which){
//             ss.seekg(which*nsamples - which);
//             std::getline(ss,line);
//             BitReader::get(ss,nsamples);
//         }

//         std::bitset<MAXBITS> read(std::istream from, int which);
// }