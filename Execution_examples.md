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

AmpliSolve requires a pre-processing step that utilizes ASEQ. To make the execution faster, this is a separate step that can be executed in parallel using a bash script for example in a cluster. ASEQ software is available at https://demichelislab.unitn.it/doku.php?id=public:aseq , but we also provide a simplified version of ASEQ named computeCounts in our Versions_used_during_development DIR.

We recommend to place all available normal samples (germline) in a separete DIR say NORMAL_BAM , and all tumour samples in another DIR named TUMOUR_BAM.



