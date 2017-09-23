System Requirements
===================

Amplisolve dependes on boost libraries, samtools and some other libraries: zlib (v1.2.3+), math, pthread.

Download boost libraries from http://www.boost.org/ untar using tar -xvf and store the folder in your system.

You need to give this path as BOOST_FLAG when you make the software (see below).

This Amplisolve distribution, generate the following binaries: 

1. computeCounts is a simple version of ASEQ software that can be used as a separate program for estimating read counts per position

2. AmpliSolvePreProVC is the pre-processing program required for the variant calling

3. AmpliSolveVC is the actual variant caller

Please make sure that always computeCounts is available and in the same directory with the other two programs.

Make also sure that you have installed samtools and bedtools in your system to run AmpliSolveCN (NOT YET AVAILABLE -> coming soon).

Compilation
===========

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


Binaries
========

We provide pre-compiled binaries of all programs generated on MacOSX (El Capitan 10.11.5) and Centos7 (7.1.1503) 


