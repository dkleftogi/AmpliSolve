## How to run AmpliSolve ?? 

First make sure that you have the correct environment, resolved all program’s dependencies and generated the binaries: 

```
AmpliSolveErrorEstimation 
```
and 

```
AmpliSolveVariantCalling
```

### Prepare the input files for AmpliSolve

AmpliSolve requires a pre-processing step that utilizes ASEQ. To make the execution faster, this is a separate step that can be executed in parallel using a bash script for example in a cluster. ASEQ software is available at https://demichelislab.unitn.it/doku.php?id=public:aseq , but we also provide a simplified version of ASEQ named computeCounts (only the binary, please contact us for more info) in our Versions_used_during_development DIR. This program deploys the PILEUP mode of ASEQ which is the one we need for this pre-processing step. 

We recommend to place all available normal samples (germline) in a separete DIR say NORMAL_BAM_DIR , and all tumour samples in another DIR named TUMOUR_BAM_DIR 

To run ASEQ (or equally computeCounts) type: 

```
./ASEQ [vcf=dummyVCF.txt] [bam=myBAM.bam] [threads=int] [mbq=int] [mrq=int] [mdc=int] [out=Out_DIR]
```

where: 
1. dummyVCF.txt is a VCF-like file that contains all positions in the gene panel. Assuming that we have a small panel with one amplicon chr8:23437012-23437020
the dummyVCF.txt will look like:

```
chr8	23437012	.	.	.	.	.	.
chr8	23437013	.	.	.	.	.	.
chr8	23437014	.	.	.	.	.	.
chr8	23437015	.	.	.	.	.	.
chr8	23437016	.	.	.	.	.	.
chr8	23437017	.	.	.	.	.	.
chr8	23437018	.	.	.	.	.	.
chr8	23437019	.	.	.	.	.	.
chr8	23437020	.	.	.	.	.	.
```
2. myBam.bam is the bam files of interest. Remember that we need to run this command for all files in NORMAL_BAM_DIR and TUMOUR_BAM_DIR 

3. threads is the number of threads to use. Usually I put 8

4. mbq, mrq , mdc are quality parameters ; we recommend the triplet 20-20-20 for Ion AmpliSeq data

5. Out_DIR is the directory to store the read count files ending with .PILEUP.ASEQ 
Since we need a DIR of normal read counts to compute the error levels, we must store all .PILEUP.ASEQ of normal samples in one DIR say NORMAL_ASEQ_DIR 

In a similar fashion, we need a DIR for all tumour read counts to perform variant calling, so we must store all .PILEUP.ASEQ of tumour samples in another DIR say TUMOUR_ASEQ_DIR 

Once we have all .PILEUP.ASEQ files for normal and tumour files stored in two separete DIRs say NORMAL_ASEQ_DIR and TUMOUR_ASEQ_DIR we are ready and we can run the AmpliSolve program. 





