njhseq
======
Version 3.0.2-dev

C++ library for dealing with dna sequence data and various other biological data created by Nicholas John Hathaway (njh)  


# Dependencies  

The majority of dependencies is downloaded by ./setup.py script but several of the libraries install depend on cmake and should be present for the install of njhseq to work. Also a C++17 compiler is needed, by default g++ on Linux (gcc 13 on Ubuntu 24.04) and clang++ on macOS.   

```
./configure.py
./setup.py --compfile compfile.mk --outMakefile makefile-common.mk
make

```
