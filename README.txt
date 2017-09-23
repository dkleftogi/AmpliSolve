AmpliSolve: detecting somatic single nucleotide variants and copy number aberrations from Ion Torrent, amplicon-based sequencing data
========================================================================================================================================================

AmpliSolve, is a new tool to identify somatic single nucleotide variants (SNVs) and copy number (CN) aberrations in matched tumour-normal samples from targeted deep sequencing experiments. AmpliSolve has been designed specifically for amplicon-based libraries sequenced with the Ion Torrent platform but it can be equally applied to other platforms (not tested yet). It relies on a population of normal samples to infer a position specific sequencing error (noise) and to estimate the expected coverage at copy neutral regions. We tested AmpliSolve on circulating tumor DNA (ctDNA) samples sequenced using a custom amplicon (AmpliSeq) panel. Our results show that AmpliSolve discriminates effectively SNVs from background noise even at low allele frequency and estimates accurately CN aberrations.       


Dependencies and System Requirements
========================================================================================================================================================

AmpliSolve dependes on the following:

1. SAM tools http://samtools.sourceforge.net/
    
2. Bed tools http://bedtools.readthedocs.io/en/latest/

Make sure that both are intalled and also included in the PATH of your system. If not you may need a command like: PATH=$PATH:/your/path/to/Samtools

During the program development and testing Samtools v1.3.1 and bedtools v2.17.9 were used.
   
3. ASEQ https://demichelislab.unitn.it/doku.php?id=public:aseq

4. Boost libraries.
    
Download boost libraries from http://www.boost.org/ untar using tar -xvf and store the folder in your system.
   
More details about AmpliSolve, and toy execution examples will be availalbe soon. 



Compilation
========================================================================================================================================================

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

Pre-compiled Binaries
========================================================================================================================================================

We provide pre-compiled binaries of all programs generated on MacOSX (El Capitan 10.11.5) and Centos7 (7.1.1503) 


