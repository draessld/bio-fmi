#ifndef TRANSFORM_H
#define TRANSFORM_H
#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/resource.h>
#include <sdsl/bit_vectors.hpp>

void TODO()
{
    std::cout << "Error: Not implemented!" << std::endl;
    exit(1);
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

int vcf2eds()
{
    TODO();
    return 0;
}

int vcf2leds()
{
    TODO();
    return 0;
}

int eds2leds(std::istream &in, std::ostream &out, unsigned int l)
{
    char ch;
    bool count = true;
    size_t cl = 0;
    bool first = true;
    std::string line = "";
    std::pair<uint64_t, uint64_t> Si = {0, 0};
    std::pair<uint64_t, uint64_t> Sj = {0, 0};

    bool merge = false;

    while (in.get(ch))
    {
        // Process each character
        switch (ch)
        {
        case '{':
            //  reset counter and stop counting
            if (cl < l && !first)
            {
                Si.second += (Si.first * line.size());
                merge = true;
            }
            else
            {
                //  ref
                merge = false;
            }

            cl = 0;
            count = false;
            first = false;
            line.clear();
            break;
        case '}':

            Sj.first++;
            Sj.second += line.size();
            // std::cout << Sj.first << "," << Sj.second << std::endl;
            // start counting
            if (merge || Si.first == 0)
            {
                if (Si.first == 0)
                {
                    Si.first = Sj.first;
                    Si.second = Sj.second;
                }
                else
                {
                    Si.second = (Sj.first * Si.second) + (Si.first * Sj.second);
                    Si.first = Si.first * Sj.first;
                }
            }
            else
            {
                //  flush
                Si.first = 0;
                Si.second = 0;
            }
            count = true;
            line.clear();
            Sj.first = 0;
            Sj.second = 0;

            // std::cout << Si.first << "," << Si.second << std::endl;
            break;
        case ',':
            Sj.first++;
            Sj.second += line.size();
            line.clear();
            break;
        default:
            if (count)
                cl++;
            line.push_back(ch);
            break;
        }
    }

    return 0;
}

#endif //   TRANSFORM_H