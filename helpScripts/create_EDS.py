"""
Grossi et al., On-line pattern matching on similar texts, CPM'17

page 11:
Synthetic data. Synthetic ED texts were created randomly (uniform distribution over
the DNA alphabet) with n ranging from 100,000 to 1,600,000; and the percentage of
degenerate positions was set to 10%. For each degenerate position within the synthetic ED
texts, the number of strings was chosen randomly, with an upper bound set to 10. The
length of each string of a degenerate position was chosen randomly, with an upper bound
again set to 10. Every non-degenerate position within the synthetic ED texts contained a
single letter. Four different patterns of length m = 8, 16, 32, or 64 were given as input to all
three programs, along with the aforementioned synthetic ED texts, resulting in four sets of
output.
"""

import random
import sys,os
from functools import reduce

sizeMB,nDegeneratePositions = sys.argv[1:-1]
sizeMB = float(sizeMB)
nDegeneratePositions = int(nDegeneratePositions)

byte_size = int(sizeMB * 1e6) # Total number of segments: 100, 500, 1000, 1600 thousands segments.
#   z velikosti souboru o pravdepodobnosti na vstupu
alphabet = "ACGTN" # Alphabet for character sampling.
psts = [0.24,0.24,0.24,0.24,0.04]   #   Probability of each character in alphabet


# Maximum number of variants (a), the number of variants for each degenerate segment will be sampled from the interval [2, a].
nMaxSegmentVariants = 10

# Maximum length of each segment variant (b), the length for each variant will be sampled from the interval [0, b] (segments might contain empty words).
nMaxVariantLength = 10

# Number of segments (must be smaller than or equal to nSegments) which are degenerate (indeterminate), i.e. contain multiple variants.
# 10% of the text as in Grossi et al.

outFile = sys.argv[3] # Output file path.

def main():

    os.remove(outFile)
        
    # nDegeneratePositions = int(prob * byte_size/10)
    degenLengths = []

    for i in range(nDegeneratePositions):
        # degenLengths.append(gen_avg(nMaxVariantLength,random.randint(2, nMaxSegmentVariants),0,20))
        degenLengths.append(random.sample(range(nMaxVariantLength), random.randint(2, nMaxSegmentVariants) )) # Randomly number of seq in degenerate symbol.
    
    symbols_size = sum([sum(i) for i in degenLengths])
    ref_size = abs(int(byte_size - symbols_size))

    print("Started, alph = \"{0}\", text size = {1}m".format(alphabet, sizeMB))

    text = randomString(alphabet, ref_size)

    degenPosList = sorted(random.sample(range(ref_size), nDegeneratePositions )) # Randomly drawn degenerate positions.
    degenStrings = {} # Dictionary: position in text -> list of a few strings

    print("Generating degenerate strings for #positions = {0}k".format(nDegeneratePositions / 1000.0))

    textPos = 0

    for curPos,curLen in zip(degenPosList,degenLengths):
        curSet = set()

        for l in curLen:
            curStr = randomString(alphabet, l)
            curSet.add(curStr)

        dumpToFile(text[textPos:curPos-1], curSet)
        textPos = curPos-1

    dumpToFile(text[textPos:],list())

def gen_avg(expected_avg,n,a,b):
    while True:
        l = [random.randint(a,b) for i in range(n)]
        avg = reduce(lambda x,y: x+y,l)/len(l)

        if avg == expected_avg:
            return l

def randomString(alph, size):
    sigma = len(alph)
    # return "".join(alph[random.randint(0, sigma - 1)] for _ in range(size))
    return "".join(random.choices(alph,psts)[0] for _ in range(size))

def dumpToFile(text, degenStrings):
    # print("Generating output elastic-degenerate text...")
    outStr = ""
    with open(outFile, "a") as f:
        f.write(text)

        if(not degenStrings == list()):
            degenStrings = sorted(list(degenStrings),key=lambda x: len(x))

            f.write('{')
            for iD in range(len(degenStrings)):
                f.write(degenStrings[iD])
                if iD != len(degenStrings)-1:
                    f.write(',')
            f.write('}')

if __name__ == "__main__":
    main()
