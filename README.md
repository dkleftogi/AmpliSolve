# AmpliSolve

Targeted deep sequencing is a highly effective technology to identify known and novel single nucleotide variants (SNVs) with many applications in translational medicine, disease monitoring and cancer profiling. However, identification of SNVs using deep sequencing data is a challenging computational problem as different sequencing artifacts limit the analytical sensitivity of SNV detection, especially at low variant allele frequencies (VAFs). To address the problem of relatively high noise levels in amplicon-based deep sequencing data (e.g. with the Ion AmpliSeq technology) in the context of SNV calling, we have developed a new bioinformatics tool called AmpliSolve. AmpliSolve uses a set of normal samples to model position-specific, strand-specific and nucleotide-specific background artifacts (noise), and deploys a Poisson model-based statistical framework for SNV detection. Our tests on both synthetic and real data indicate that AmpliSolve achieves a good trade-off between precision and sensitivity, even at VAF below 5% and as low as 1%. We further validate AmpliSolve by applying it to the detection of SNVs in 96 circulating tumor DNA samples at three clinically relevant genomic positions and compare the results to digital droplet PCR experiments. AmpliSolve is a new tool for in-silico estimation of background noise and for detection of low frequency SNVs in targeted deep sequencing data. Although AmpliSolve has been specifically designed for and tested on amplicon-based libraries sequenced with the Ion Torrent platform it can, in principle, be applied to other sequencing platforms as well. 


## Dependencies and System Requirements

AmpliSolve depends on the following:

```
1. SAM tools downloaded from http://samtools.sourceforge.net/
``` 

Make sure that is installed and included in the PATH of your system. If not you may need a command like: 

```
PATH=$PATH:/your/path/to/Samtools
```

During the program development we used Samtools v1.3.1 (provided with the original licences)
   
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

This project is licenced under the GNU GPLv3 Licence.

Copyright (C) 2018 The Institute of Cancer Research (ICR) -- Dimitrios Kleftogiannis

You may not use this file except in compliance with the License. A copy of the licence is availle with the source codes.

## Contact

If you face any compilation problems, issues with dependencies, versions and so on, or you found bugs please contact me at dimitrios DOT kleftogiannis AT icr DOT ac DOT uk 

We would also appreciate hearing about how you used this code, improvements that you have made to it. 




