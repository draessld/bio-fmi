#include <iostream>
#include <bitset>
#include <stdexcept>

#define MAXBITS 1024

class BitReader
{
private:
    std::istream &input;
    unsigned char buffer;
    int bitsLeft;

public:
    BitReader(std::istream &in) : input(in), buffer(0), bitsLeft(0) {}

    int getBit()
    {
        if (bitsLeft == 0)
        {
            input.read(reinterpret_cast<char *>(&buffer), 1);
            if (input.eof())
            {
                throw std::runtime_error("End of stream reached");
            }
            bitsLeft = 8;
        }
        int bit = (buffer >> (bitsLeft - 1)) & 1;
        bitsLeft--;
        return bit;
    }

    std::bitset<MAXBITS> get(std::istream &in, int nbits)
    {
        std::bitset<MAXBITS> result;
        for (int i = 0; i < nbits; ++i)
        {
            result[i] = getBit();
        }
        return result;
    }
};

class BitWriter
{
private:
    std::ostream &output;
    unsigned char buffer;
    int bitsFilled;

public:
    BitWriter(std::ostream &out) : output(out), buffer(0), bitsFilled(0) {}

    void writeBit(int bit)
    {
        buffer = (buffer << 1) | (bit & 1);
        bitsFilled++;
        if (bitsFilled == 8)
        {
            flush();
        }
    }

    void write(std::ostream &out, const std::bitset<MAXBITS> &bits, int nbits)
    {
        for (int i = 0; i < nbits; ++i)
        {
            writeBit(bits[i]);
        }
    }

    void flush()
    {
        if (bitsFilled > 0)
        {
            buffer <<= (8 - bitsFilled);
            output.write(reinterpret_cast<char *>(&buffer), 1);
            buffer = 0;
            bitsFilled = 0;
        }
    }

    ~BitWriter()
    {
        flush();
    }
};

