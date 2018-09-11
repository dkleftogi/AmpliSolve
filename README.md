# AmpliSolve

AmpliSolve, is a new bioinfrmatics tool for precise error estimation and variant calling in targeted deep sequencing data. AmpliSolve has been designed specifically for amplicon-based libraries sequenced with the Ion AmpliSeq technology, and it can be in principle be applied to other platforms (but not tested yet). AmpliSolve uses a population of normal samples as a training set to infer a position specific, nucleotide specific and strand specific background noise levels, and deploys a Poisson-based statistical model to identify single nucleotide variants. We have tested AmpliSolve on circulating tumor DNA (ctDNA) samples sequenced using a custom AmpliSeq panel. 


## Dependencies and System Requirements

AmpliSolve dependes on the following:

```
1. SAM tools http://samtools.sourceforge.net/
``` 
```
2. Bed tools http://bedtools.readthedocs.io/en/latest/
```

Make sure that both are installed and included in the PATH of your system. If not you may need a command like: 

```
PATH=$PATH:/your/path/to/Samtools
```

During the program development we used Samtools v1.3.1 and bedtools v2.17.9 (both provided with the original licences)
   
```
3. ASEQ software downloaded from https://demichelislab.unitn.it/doku.php?id=public:aseq
```
We use ASEQ to generate read count files ending with .PILEUP.ASEQ, from the original BAM files. AmpliSolve does not work on the actual BAM files, and thus it requires to pre-process the input and generate the input in the appropriate format. Examples are provided below.  

```
4. Boost libraries from http://www.boost.org/
```

Download the libraries, untar using tar -xvf and simply store the folder in your system.
   

Compilation
====================================================================================================

This is a C++ program, thus you need the C++ compiler to be installed in your computer.

The AmpliSolve program has been developed in a Mac OS computer with El Capitan version 10.11.5

The C++ compiler used for development is the clang C, C++ and Objective-C compiler

The program works for Unix-like systems and has been tested but not extensively. 

Samtools library `libbam.a' has been generated for GNU/Linux and MacOSX systems.
If libraries are not working we suggest to download them again.

The program might work on Windows under cygwin emulator, however it has not been tested.

Linux x64:
If your kernel is >= 2.6.15 type `make -fMakefile.linux BOOST_FLAG=/your/dir/to/boost/libraries' to compile all programs.

Linux i386:
If your kernel is >= 2.6.15 type `make -fMakefile.linux i386 BOOST_FLAG=/your/dir/to/boost/libraries' to compile all programs.

MacOSX:
Requires XCode.
Type `make -fMakefile.macosx BOOST_FLAG=/your/dir/to/boost/libraries' to compile all programs.

You will see the following:

gcc -g -w -Wall -O2 -I./include/ main.c -o computeCounts -L./lib/macosx/ -lbam -lm -lz -lpthread
cc -I /Users/dkleftog/boost_1_61_0 -o AmpliSolvePreProVC AmpliSolvePreProVC.cpp -std=c++0x -lstdc++
cc -I /Users/dkleftog/boost_1_61_0 -o AmpliSolveVC AmpliSolveVC.cpp -std=c++0x -lstdc++
cc -o AmpliSolveCN AmpliSolveCN.cpp -std=c++0x -lstdc++
cc -o AmpliSolvePreProCN AmpliSolvePreProCN.cpp -std=c++0x -lstdc++

and if the compilation is successful you will find the following binaries:
computeCounts, AmpliSolvePreProVC, AmpliSolveVC, AmpliSolveCN and AmpliSolvePreProCN

To make the program accessible from any location in the system you need to include all binaries created in your PATH.

This can be done easily with a command like PATH=$PATH:/write/your/dir/to/AmpliSolve/source_codes

If you face any problem with the dependencies, versions or other compilation issues 
please contact the developer at dimitrios DOT kleftogiannis AT icr DOT ac DOT uk 



Pre-compiled Binaries
====================================================================================================

We provide pre-compiled binaries of all AmpliSolve programs generated for MacOSX (El Capitan 10.11.5) and Centos7 (7.1.1503).

Please contact directly Dimitrios at dimitrios DOT kleftogiannis AT icr DOT ac DOT uk

