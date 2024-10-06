#include <iostream>
#include <fstream>
#include <sstream>

#include "utils/eds.h"

using namespace std;

int main(int argc, char const *argv[])
{
    // Check if an input argument is provided
    if (argc < 2) {
        cout << "Error: No input provided!" << endl;
        return 1;
    }

    string input = argv[1];

    ifstream ifs(input);
    if (!ifs)  // If the file could not be opened, treat input as a sequence
    {
        stringstream iss(input);
        EDS eds(iss);  // Pass the stringstream as istream
        eds.stats();
    }
    else  // If the file can be opened, treat input as a file
    {
        EDS eds(ifs);  // Pass the ifstream as istream
        eds.stats();
    }
    
    return 0;
}