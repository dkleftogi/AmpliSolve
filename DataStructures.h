#ifndef DATASTRUCTURES_H_INCLUDED
#define DATASTRUCTURES_H_INCLUDED

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "samtools/sam.h"

/////////////////////////////////////////////////////////////////////////////////////
// Structures Definition
/////////////////////////////////////////////////////////////////////////////////////

// Running threads list
struct node
{
   int val;
   int thread;
   struct node *next;
};

typedef struct node threads_alloc_list;


// Contains pileup info of a SNP
struct snps_t
{
    char *chr;
    uint32_t pos;
    char *dbsnp;
    char *maf;
    char *ref;
    char *ref_default;
    char *alt;
    char *alt_default;
    int A;
    int C;
    int G;
    int T;
    int As;
    int Cs;
    int Gs;
    int Ts;
    double ASE_pvalue;
    double FISHER_pvalue;
    int control_ref;
    int control_alt;
};

// Collects pileup results for all SNPs
struct snps_info
{
    int length;
    struct snps_t  **info;
};

// Pileup region coordinates
typedef struct {
        int beg, end;
        samfile_t *in;
} tmpstruct_t;

// Collects info to retrieve SNPs info during alignment pileup
typedef struct
{
        int index;
        struct snps_info *snps;
        tmpstruct_t *tmp;
        bam_pileup1_t *plp[8000];
        int counter;
} snpsstruct_t;

// Used to pass argument to threaded PileUp function
struct args_thread
{
    int start;
    int end;
    struct snps_info *snps;
    char* bamfile;
    bam_index_t *idx;
    struct node *map;
};

// Contains info of a gene
struct genes_t
{
    char *gene;
    char *chr;
    uint32_t from;
    uint32_t to;
    double *samples;
    int *samples_snps;
    uint32_t *init_intervals;
    uint32_t *end_intervals;
    int intervals_number;
};

// Collects info of all the genes
struct genes_info
{
    int length;
    struct genes_t **info;
};

typedef struct
{
	int k, x, y, end;
} cstate_t;

struct input_args
{
    char *mode; // 1-PileUp,2-ASE
    char *bed;
    char *genes;
    char *transcripts;
    char *vcf;
    char *vcflist;
    char *bam;
    char *bamlist;
    char **bamfiles;
    bam_index_t **bamindexes;
    char **vcffiles;
    int bamfiles_number;
    int vcffiles_number;
    int cores;
    int mbq;
    int mrq;
    int mrn;
    int granularity;
    int gt_mode;
    double HT_pvalue;
    double FT_pvalue;
    int HT_zeros;
    char *outdir;
    int gene_stat;
    double snp_stat;
    double HT_perc;
    int HT_binom;
    int pileup_blocks;
    int pileup_mode;
};

#endif // DATASTRUCTURES_H_INCLUDED
