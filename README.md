# BIO-FMI

Author: Dominika Draesslerová (dominikadraesslerova@gmail.com)

From the paper: "Draesslerová, Dominika. Bioinformatics index tool for elastic degenerate string matching. Master’s thesis. Czech Technical University in Prague, Faculty of Information Technology, 2022. Also available from: ⟨https://github.com/draessld/bio-fmi⟩"

##  Brief description

This software is called BIO-FMI (Index for set of genomes and pangenomes). It can be used for pattern matching in elastic-degenerate (ED) text. Authors: Petr Prochazka, Jan Holub.
Input format can be saved in .aln (ALN) file - alignment of sequences, every sequence in new line; or in .eds file - {A,C,}GAAT{,A,AT}ATT, where the strings in bracket are ordered ascending. BIO-FMI return start positions of pattern occurences with sequence index; or positions of occurences with absolute number sequence in ED symbol.

##  Download

To clone the repository, use following option:
```
$   git clone https://github.com/draessld/bio-fmi
```
##  Compile
You need the SDSL library installed on your system (https://github.com/simongog/sdsl-lite).

bio-fmi uses cmake to generate the Makefile. Create a build folder in the main bio-fmi folder:
```
$   mkdir build
```
execute cmake:
```
$   cd build; cmake ..
```
and compile:
```
$   make
```
##  Run
After compiling, run

```
$   bio-fmi-build input.txt
```

This command will create the BIO-FMI index of the text file "input.aln" or "input.eds" and will store it using as filename as prefix. Use option -o to specify a different output folder for the index files.

Run

```
$   bio-fmi-locate index_basename patterns_file
```

to locate a set of patterns in followinf format. Generated using script *genpattern* from pizza&chilli (http://pizzachili.dcc.uchile.cl/experiments.html) testbed.
```
# number=3 length=5 file=pattern_file forbidden=\n
GTGCNAGCGAAGCTA
```

##  TODO
-   ALN repair mismatches
-   EDS text reading time booster
-   Experiments with other index over EDS
-   Experiments with real datasets
-   Chromosome selection


