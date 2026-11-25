#ifndef UTILS_H
#define UTILS_H
#include <string>
#include <iostream>
#include <fstream>
#include <sys/resource.h>
#include <sdsl/bit_vectors.hpp>
#include <vector>
#include <type_traits>
#include <bitset>
// #include "SourceManager.cpp"

#define MAXBITS 1024

void TODO()
{
    std::cout << "Error: Not implemented!" << std::endl;
    exit(1);
}

void printch(std::set<std::string> &a)
{
    for (auto &c : a)
        std::cout << c << ", ";
    std::cout << std::endl;
}

void printch(std::vector<std::string> &a)
{
    for (auto &c : a)
        std::cout << c << ", ";
    std::cout << std::endl;
}

std::set<std::string> cartesian(std::set<std::string> &a, std::set<std::string> &b)
{
    std::set<std::string> res;
    //  cartesian product
    for (auto ai : a)
        for (auto bi : b)
            res.insert(ai + bi);
    if (res.empty())
        return b;

    return res;
}

std::vector<std::string> match_variants(std::vector<std::string> &a, std::vector<std::string> &b, std::istream &in_source, int var_counter, int nsources)
{
    //  check if tmp file exists - change var_count and read from it
    const std::filesystem::path tmpfile("./tmp.out");
    std::fstream tmp(tmpfile,std::ios::in | std::ios::out | std::ios::binary);
    tmp.seekg(0, std::ios::end);
    std::streampos endPos = tmp.tellg();
    std::streampos pos = endPos;
    if(!tmp.is_open()){
        std::cerr<<"Error: open tmp files was not possible" << std::endl;
    }
    std::vector<std::string> result;
    std::vector<std::string> outps;
    std::string line;
    for (int x = 0; x < b.size(); x++)
    {
        // std::cout << "b:" << x << "=" << b[x] << std::endl;
        // std::cout << (var_counter - (b.size()-x)) * nsources + (var_counter-(b.size()-x)) << std::endl;
        in_source.seekg((var_counter - (b.size()-x)) * nsources + (var_counter-(b.size()-x)), std::ios::beg);
        std::getline(in_source,line);
        std::bitset<MAXBITS> bitsb(line);
        for (int y = 0; y < a.size(); y++)
        {
            // std::cout << "a:" << y << "=" << a[y] << std::endl;
            // std::cout << (var_counter - b.size()-(a.size()-y)) * nsources + (var_counter - b.size()-(a.size()-y)) << std::endl;
            tmp.seekg((var_counter - b.size()-(a.size()-y)) * nsources + (var_counter - b.size()-(a.size()-y)), std::ios::beg); //    jump onto line with sources
            std::getline(tmp, line);
            std::bitset<MAXBITS> bitsa(line);
            auto res = bitsa & bitsb; //  which sources they share
            if (!res.none())
            {
                //  add variant merge
                result.push_back(a[y] + b[x]);
                outps.push_back(res.to_string().substr(MAXBITS-nsources,nsources));
            }
        }
    }

    tmp.seekp((var_counter - b.size()-a.size()) * nsources + (var_counter - b.size()-a.size()), std::ios::beg);    //  TODO move to the position in the file
    for (auto source: outps)
    {
        tmp << source << std::endl; //  rewrite old sources
    }
    tmp.close();

    return result;
}

/*  Transformation fucntion */
int msa2eds()
{
    TODO();
    return 0;
}

int msa2leds_linear(std::istream &in, std::ostream &out, unsigned int l)
{
    //  get the first "reference" sequence to compare with
    std::string line;
    uint64_t counter = 0;
    std::string first;
    sdsl::bit_vector B;
    size_t i = 0;
    int x = -1; //  to setup the width of fasta

    std::vector<int> start_positions;

    while (std::getline(in, line))
    {
        if (line.empty()) // Skip empty lines
            continue;

        if (line[0] == '>')
        {
            // reset
            if (counter == 1)
            {
                //  initialize bit vector
                B = sdsl::bit_vector(first.size() + 1, 1);
            }
            i = 0;
            counter++;
            start_positions.push_back(in.tellg());
        }
        else if (counter == 1)
        {
            //  first sequence
            first += line;
            if (x == -1)
                x = line.size();
        }
        else
        {
            // compare with the same position in first sequence and switch bit vector if not
            for (size_t j = 0; j < line.size(); j++)
            {
                if (line[j] != first[i] || line[j] == '-')
                {
                    B[i] = 0;
                }
                i++;
            }
        }
    }

    B[first.size()] = B[first.size() - 1] ^ 1;

    // std::cout << first << std::endl;
    // std::cout << B << std::endl;
    // for (size_t j = 0; j < start_positions.size(); j++)
    // {
    //     std::cout << start_positions[j] << ',';
    // }
    //     std::cout <<std::endl;

    sdsl::bit_vector::select_1_type selectOne(&B);
    sdsl::bit_vector::select_0_type selectZero(&B);

    size_t next_one = 0, next_zero = 0;
    size_t zeros = 0, ones = 0;
    size_t cl = 0;

    next_one = selectOne(1);
    next_zero = selectZero(1);

    sdsl::bit_vector H(B.size(), 0);

    i = 0;
    while (i < first.size())
    {
        // std::cout << i << ":" << B[i] << std::endl;

        if (B[i])
        {
            // match
            next_zero = selectZero(zeros + 1);
            // std::cout << "next zero: " << next_zero << ", total zeros: " << zeros << std::endl;

            cl = next_zero - i;
            // std::cout << "run of 1s length: " << cl << std::endl;
            if (cl >= l || i == 0 || next_zero == first.size())
            {
                // std::cout << "here is common part" << std::endl;
                for (size_t j = i; j < next_zero; j++)
                    H[j] = 1;
            }
            ones += cl;
            i = next_zero;
        }
        else
        {
            //  mismatch
            next_one = selectOne(ones + 1);
            // std::cout << "next one: " << next_one << ", total ones: " << ones << std::endl;
            cl = next_one - i;
            // std::cout << "run of 0s length: " << cl << std::endl;
            zeros += cl;
            i = next_one;
        }
    }

    H[first.size()] = H[first.size() - 1] ^ 1;
    // std::cout << H << std::endl;

    selectOne = decltype(selectOne)(&H);
    selectZero = decltype(selectZero)(&H);

    next_one = 0;
    next_zero = 0;
    zeros = 0;
    ones = 0;
    cl = 0;
    int tmp;
    char *buffer = new char[first.size() + (first.size() / x)];

    in.clear();
    in.seekg(0);

    std::set<std::string> changes;
    std::string change;

    // std::cout << H << std::endl;

    i = 0;
    while (i < first.size())
    {
        if (H[i])
        {
            // std::cout << "match" << std::endl;
            // match - read cl characters from first sequence
            next_zero = selectZero(zeros + 1);

            cl = next_zero - i;

            ones += cl;
            // Move to specific position
            // std::cout << "moving to: " << start_positions[0] + i + (int)(i / x) << " position" << std::endl;
            tmp = (int)((i % x) + cl) / x;

            in.seekg(start_positions[0] + i + (int)(i / x));

            in.read(buffer, cl + tmp);
            buffer[cl + tmp] = '\0';

            counter = 0;

            out << "{";
            for (size_t k = 0; k < cl + tmp; k++)
                if (buffer[k] != '\n' && buffer[k] != '-')
                {
                    counter++;
                    out << buffer[k];
                }
            out << "}";
            i = next_zero;
        }
        else
        {
            //  mismatch
            next_one = selectOne(ones + 1);

            cl = next_one - i;

            out << '{';

            for (size_t j = 0; j < start_positions.size(); j++)
            {
                // std::cout << i + cl << std::endl;
                in.seekg(start_positions[j] + i + static_cast<int>(i / x));
                tmp = static_cast<int>((i % x) + cl) / x;

                in.read(buffer, cl + tmp);
                buffer[cl + tmp] = '\0';

                // std::cout << first.size() << std::endl;
                // std::cout << cl+tmp << std::endl;

                for (size_t k = 0; buffer[k] != '\0'; ++k)
                // for (size_t k = 0; k < cl + tmp; k++)
                {
                    if (buffer[k] != '\n' && buffer[k] != '-')
                    {
                        // std::cout << buffer[k];
                        change.push_back(buffer[k]);
                        // if (save)
                        //     out << buffer[k];
                        // N++;
                    }
                }
                changes.insert(change);
                change.clear();
                // if (save)
                //     out << ',';
                // m++;
            }
            for (auto ch : changes)
            {
                out << ch;
                out << ',';
            }
            changes.clear();
            tmp = out.tellp();
            out.seekp(tmp - 1);

            out << '}';

            zeros += cl;
            i = next_one;
        }
    }

    // delete[] buffer;

    return 0;
}

int msa2leds_cartesian(std::istream &in, std::ostream &out, unsigned int l)
{
    //  get the first "reference" sequence to compare with
    std::string line;
    int counter = 0;
    std::string first;
    sdsl::bit_vector B;
    size_t i = 0;
    int x = -1; //  to setup the width of fasta

    std::vector<int> start_positions;

    while (std::getline(in, line))
    {

        if (line.empty()) // Skip empty lines
            continue;

        if (line[0] == '>')
        {
            // reset
            if (counter == 1)
            {
                //  initialize bit vector
                B = sdsl::bit_vector(first.size() + 1, 1);
            }
            i = 0;
            counter++;
            start_positions.push_back(in.tellg());
        }
        else if (counter == 1)
        {
            //  first sequence
            first += line;
            if (x == -1)
                x = line.size();
        }
        else
        {
            // compare with the same position in first sequence and switch bit vector if not
            for (size_t j = 0; j < line.size(); j++)
            {
                if (line[j] != first[i] || line[j] == '-')
                {
                    B[i] = 0;
                }
                i++;
            }
        }
    }

    B[first.size()] = B[first.size() - 1] ^ 1;

    // std::cout << B << std::endl;
    sdsl::bit_vector::select_1_type selectOne(&B);
    sdsl::bit_vector::select_0_type selectZero(&B);

    size_t next_one = 0, next_zero = 0;
    size_t zeros = 0, ones = 0;
    size_t cl = 0;

    next_one = selectOne(1);
    next_zero = selectZero(1);

    next_one = 0;
    next_zero = 0;
    zeros = 0;
    ones = 0;
    cl = 0;
    int tmp;
    std::set<std::string> changes_o;
    std::set<std::string> changes_n;
    char *buffer;
    buffer = new char[first.size()];
    line.clear();

    in.clear();

    i = 0;
    while (i < first.size())
    {
        if (B[i])
        {
            // match - read cl characters from first sequence
            next_zero = selectZero(zeros + 1);

            cl = next_zero - i;
            ones += cl;

            // Move to specific position
            tmp = (int)((i % x) + cl) / x;
            in.seekg(start_positions[0] + i + (int)(i / x));

            in.read(buffer, cl + tmp);
            buffer[cl + tmp] = '\0';

            for (size_t k = 0; k < cl + tmp; k++)
                if (buffer[k] != '\n' && buffer[k] != '-' && buffer[k] != '\0')
                {
                    line.push_back(buffer[k]);
                }

            if (next_zero == first.size())
            {
                break;
            }

            if (i == 0 && !line.empty())
            {
                out << '{';
                out << line;
                out << '}';
            }
            else if (cl < l)
            { //  neproddo not output common part, but prolong changes
                for (auto ch : changes_o)
                {
                    changes_n.insert(ch + line);
                }
            }
            else
            { //  print changes + common part

                //  flush data into file

                out << '{';
                for (const std::string &str : changes_o)
                {
                    out << str; // Write each string followed by a newline
                    out << ',';
                }
                tmp = out.tellp();
                out.seekp(tmp - 1);

                if (tmp != 1)
                    out << '}';

                out << '{';
                out << line;
                out << '}';

                changes_o.clear();
                line.clear();
            }
            i = next_zero;
            changes_o = changes_n;
            changes_n.clear();
            line.clear();
        }
        else
        {
            //  mismatch
            next_one = selectOne(ones + 1);
            cl = next_one - i;

            for (size_t j = 0; j < start_positions.size(); j++)
            {
                in.seekg(start_positions[j] + i + (int)(i / x));
                tmp = (int)((i % x) + cl) / x;

                in.read(buffer, cl + tmp);
                buffer[cl + tmp] = '\0';

                line.clear();
                for (size_t k = 0; k < cl + tmp; k++)
                    if (buffer[k] != '\n' && buffer[k] != '-' && buffer[k] != '\0')
                    {
                        line.push_back(buffer[k]);
                    }

                changes_n.insert(line);
            }

            changes_o = cartesian(changes_o, changes_n);
            changes_n.clear();

            zeros += cl;
            i = next_one;
            line.clear();
        }
    }

    if (!changes_o.empty())
    {
        //  flush data into file
        out << '{';
        for (const std::string &str : changes_o)
        {

            out << str; // Write each string followed by a newline
            out << ',';
        }
        tmp = out.tellp();
        out.seekp(tmp - 1);
        if (tmp != 1)
            out << '}';

        if (!line.empty())
        {
            out << '{';
            out << line;
            out << '}';
        }
        changes_o.clear();
    }
    delete[] buffer;
    buffer = nullptr;

    return 0;
}

int vcf2eds(std::istream &fastas, std::istream &vcfs, std::ostream &out, unsigned int l)
{

    return 0;
}

int vcf2leds()
{
    TODO();
    return 0;
}

int eds2leds_cartesian(std::istream &in, std::ostream &out, unsigned int l)
{
    size_t cl = l - 1;
    bool first = true;
    bool merge = false;
    std::string line = "";
    std::set<std::string> a;
    std::set<std::string> b;
    int counter = 0;

    char ch;
    while (in.get(ch))
    {
        switch (ch)
        {
        case '{':
            // a set is open - may be degenerate or reference
            if (!line.empty())
            {
                // reference part without {}
                if (line.size() < cl && !first)
                {
                    //  need to merges
                    merge = true;
                    for (auto change : a)
                    {
                        change.append(line);
                    }
                }
                else
                {
                    //  flush reference part
                    if (a.size() > 1)
                    {
                        out << '{';
                        for (auto change : a)
                        {
                            out << change << ',';
                        }
                        out.seekp(-1, std::ios::cur);
                        out << '}';
                    }
                    else if (a.size() == 1)
                    {
                        out << *(a.begin());
                    }
                    out << line;
                }
                first = false;
                line.clear();
            }
            counter = 1;
            break;
        case '}':
            b.insert(line);
            if (counter == 1)
            {
                //  reference
                if (line.size() < cl && !first)
                {
                    //  need to merges
                    merge = true;
                    a = cartesian(a, b);
                }
                else
                {
                    //  flush reference part
                    if (a.size() > 1)
                    {
                        out << '{';
                        for (auto change : a)
                        {
                            out << change << ',';
                        }
                        out.seekp(-1, std::ios::cur);
                        out << '}';
                    }
                    else if (a.size() == 1)
                    {
                        out << *(a.begin());
                    }
                    out << line;
                }
                first = false;
            }
            else
            {
                //  degenerate symbol
                if (merge)
                {
                    a = cartesian(a, b);
                    merge = false;
                }
                else
                {
                    if (first){
                        //  store a into tmp
                        // std::ofstream tmp("./tmp.out");
                        // for (int i=0;i<a.size();i++){
                        //     in_source.seekg(i*nsources + i);
                        //     std::getline(in_source,line);
                        //     tmp << line;
                        // }
                        // tmp.close();
                        first = false;
                    }
                    a = b;
                }
            }
            b.clear();
            line.clear();
            break;
        case ',':
            b.insert(line);
            line.clear();
            counter++;
            break;
        case '\n':
            break;
        default:
            line.push_back(ch);
            break;
        }
    }
    if (a.size() > 1)
    {
        out << '{';
        for (auto change : a)
        {
            out << change << ',';
        }
        out.seekp(-1, std::ios::cur);
        out << '}';
    }

    if (b.size() > 1)
    {
        out << '{';
        for (auto change : b)
        {
            out << change << ',';
        }
        out.seekp(-1, std::ios::cur);
        out << '}';
    }
    else if (b.size() == 1)
    {
        out << *(b.begin());
    }

    if (!line.empty())
    {
        // reference part without {}
        out << line;
    }

    return 0;
}

int eds2leds_linear(std::istream &in, std::ostream &out, unsigned int l, std::istream &in_source, std::ostream &out_source)
{
    std::cout << "linear EDS2lEDS" << std::endl;
    const std::filesystem::path tmpfile("./tmp.out");
    std::ofstream tmp;
    size_t cl = l - 1;
    bool first = true;
    bool merge = false;
    std::string line = "";
    std::string slinea, slineb;
    std::vector<std::string> a;
    std::vector<std::string> b;
    int counter = 0;
    int var_counter = 0; //  to track number of

    std::filesystem::remove(tmpfile);

    std::getline(in_source, slinea); //  first variant is always reference
    int nsources = slinea.size();

    char ch;
    while (in.get(ch))
    {
        switch (ch)
        {
        case '{':
            // a set is open - may be degenerate or reference
            if (!line.empty())
            {
                // reference part without {}
                if (line.size() < cl && !first)
                {
                    //  need to merges
                    merge = true;
                    for (auto change : a)
                        change.append(line);
                }
                else
                {
                    //  flush reference part
                    if (a.size() > 1)
                    {
                        out << '{';
                        for (auto change : a)
                        {
                            out << change << ',';
                        }
                        out.seekp(-1, std::ios::cur);
                        out << '}';
                    }
                    else if (a.size() == 1)
                    {
                        out << *(a.begin());
                    }
                    out << line;
                }
                first = false;
                line.clear();
            }

            counter = 1;
            break;
        case '}':
            b.push_back(line);
            if (counter == 1)
            {
                //  reference with {}
                if (line.size() < cl && !first)
                {
                    //  need to merges
                    merge = true;
                    //append reference to the variant in a
                    for (auto &change: a)
                        change.append(line);     
                }
                else
                {
                    //  flush reference part
                    if (a.size() > 1)
                    {

                        // std::cout << "Flushing previous set: " << std::endl;
                        // printch(a);

                        out << '{';
                        for (auto change : a)
                        {
                            out << change << ',';
                        }
                        out.seekp(-1, std::ios::cur);
                        out << '}';

                    }
                    else if (a.size() == 1)
                    {
                        //????  reference stored in change set
                        out << *(a.begin());
                    }
                    out << line;
                }
                first = false;
            }
            else
            {
                //  degenerate symbol
                if (merge)
                {                    
                    var_counter -= a.size();
                    a = match_variants(a, b, in_source, var_counter+a.size()+b.size(), nsources);
                    var_counter += a.size();
                    merge = false;
                }
                else
                {
                    var_counter+=counter;
                    //  store a into tmp
                    tmp.open(tmpfile,std::ios::app);
                    if (first){
                        tmp.seekp(0,std::ios::beg);
                        for (int i=0;i<b.size();i++){
                            in_source.seekg(i*nsources + i);
                            std::getline(in_source,line);
                            tmp << line << std::endl;
                        }
                        first = false;
                    }else{
                        tmp.seekp(0,std::ios::end); //  start on the end of file
                        for (int i=b.size();i>0;i--){
                            in_source.seekg((var_counter-i)*nsources + (var_counter-i));
                            std::getline(in_source,line);
                            tmp << line << std::endl;
                        }
                    }
                    tmp.close();
                    a = b;
                }
            }
            b.clear();
            line.clear();
            break;
        case ',':
            // var_counter++;
            b.push_back(line);
            line.clear();
            counter++;
            break;
        case '\n':
            break;
        default:
            line.push_back(ch);
            break;
        }
    }
    if (a.size() > 1)
    {
        out << '{';
        for (auto change : a)
        {
            out << change << ',';
        }
        out.seekp(-1, std::ios::cur);
        out << '}';
    }

    if (b.size() > 1)
    {
        out << '{';
        for (auto change : b)
        {
            out << change << ',';
        }
        out.seekp(-1, std::ios::cur);
        out << '}';
    }
    else if (b.size() == 1)
    {
        out << *(b.begin());
    }

    if (!line.empty())
    {
        // reference part without {}
        out << line;
    }

    return 0;
}

#endif //   UTILS_H