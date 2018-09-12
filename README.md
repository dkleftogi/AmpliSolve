# AmpliSolve

AmpliSolve, is a new tool for precise error estimation and variant calling in targeted deep sequencing data. AmpliSolve has been designed specifically for amplicon-based libraries sequenced with the Ion AmpliSeq technology, and it can in principle be applied to other platforms. AmpliSolve uses a population of normal samples as a training set to infer position specific, nucleotide specific and strand specific background noise levels, and deploys a Poisson-based statistical model to identify single nucleotide variants. We have tested AmpliSolve on circulating tumor DNA (ctDNA) samples sequenced using a custom AmpliSeq panel. 


## Dependencies and System Requirements

AmpliSolve depends on the following:

```
1. SAM tools downloaded from http://samtools.sourceforge.net/
``` 

Make sure that is installed and included in the PATH of your system. If not you may need a command like: 

```
PATH=$PATH:/your/path/to/Samtools
```

During the program development we used Samtools v1.3.1 and bedtools v2.17.9 (both provided with the original licences)
   
```
2. ASEQ software downloaded from https://demichelislab.unitn.it/doku.php?id=public:aseq
```
AmpliSolve runs without ASEQ but we use ASEQ to generate read count files ending with .PILEUP.ASEQ, from the original BAM files. 

AmpliSolve does not work on the actual BAM files, and thus it requires to pre-process the input BAMs and generate the input in the appropriate format. Examples are provided at Execution_examples.md  

```
3. Boost libraries from http://www.boost.org/
```

Download the libraries, untar using tar -xvf and simply store the folder in your system.
   

## Compilation

AmpliSolve is a C++ program, thus you need the C++ compiler to be installed in your computer. The program has been developed in a Mac OS computer with El Capitan version 10.11.5. The C++ compiler used for development is the clang C, C++ and Objective-C compiler.The program works for Unix-like systems. The program might work on Windows under cygwin emulator, however it has not been tested.

### To install AmpliSolve type one of the following:

Linux i386 with kernel >= 2.6.15  
```
make -fMakefile.linux i386 BOOST_FLAG=/your/dir/to/boost/libraries
``` 

MacOSX with XCode
```
make -fMakefile.macosx BOOST_FLAG=/your/dir/to/boost/libraries
```


You will see two commands similar to the following:

```
cc -I /Users/dkleftog/boost_1_61_0 -o AmpliSolvePreProVC AmpliSolvePreProVC.cpp -std=c++0x -lstdc++
cc -I /Users/dkleftog/boost_1_61_0 -o AmpliSolveVC AmpliSolveVC.cpp -std=c++0x -lstdc++
```

and if the compilation is successful you will find the following binaries:
AmpliSolveErrorEstimation and AmpliSolveVariantCalling

To make the program accessible from any location in the system you need to include all binaries in your PATH. This can be done easily with a command like 

```
PATH=$PATH:/write/your/dir/to/AmpliSolve/source_codes
```

## Licence

This project is licenced under the the Educational Community License, Version 2.0. You may not use this file except in compliance with the License. You may obtain a copy of the License at https://opensource.org/licenses/ECL-2.0

## Contact

If you face any compilation problems, issues with dependencies, versions and so on, or you found bugs please contact me at dimitrios DOT kleftogiannis AT icr DOT ac DOT uk 

We would also appreciate hearing about how you used this code, improvements that you have made to it. 




