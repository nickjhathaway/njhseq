njhseq
======
Version 2.6.2-dev

C++ library for dealing with dna sequence data and various other biological data created by Nicholas John Hathaway (njh)  


# Dependencies  

The majority of dependencies is downloaded by ./setup.py script but several of the libraries install depend on cmake and should be present for the install of njhseq to work. Also a modern day c++ (clang++-3.8 or g++-6) compilier is needed.   

```
./configure.py 
./setup.py --compfile compfile.mk --outMakefile makefile-common.mk 
make 

```

