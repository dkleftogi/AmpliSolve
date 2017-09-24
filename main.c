/*
This code has been taken from ASEQ software implemented by Romanel A.

The code has been incorporated to AmpliSolve software with modifications on the way it runs.

*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/stat.h>
#include "samtools/sam.h"
#include "Utilities.h"
#include "DataStructures.h"
#include "Statistics.h"
#include "ParseParameters.h"
#ifdef _WIN32
# include <windows.h>
#define sleep(x) Sleep(1000 * x)
#endif

// basic definitions
int BASE_QUALITY=0;
int READ_QUALITY=0;
int READ_NUMBER=1;
int PRINT_MODE=0;

// version
int VERSION1=1;int VERSION2=1;int VERSION3=11;

// chromosomes indexing structure
int chr_indexing[25];

/////////////////////////////////////////////////////////////////////////////////////
// Auxiliar functions definition
/////////////////////////////////////////////////////////////////////////////////////

// Increments the base of SNP pileup given the codified value of a base
void incBase(struct snps_t *elem, int val)
{
    if (val==1)
        elem->A++;
    else if (val==2)
        elem->C++;
    else if (val==4)
        elem->G++;
    else if (val==8)
        elem->T++;
}

// Increments the base of SNP pileup given the codified value of a base
void incBaseStrand(struct snps_t *elem, int val, int strand)
{
    if (val==1 && strand)
        elem->As++;
    else if (val==2 && strand)
        elem->Cs++;
    else if (val==4 && strand)
        elem->Gs++;
    else if (val==8 && strand)
        elem->Ts++;
}

// Return 1 if position is in interval
int isInInterval(struct snps_t *snp, struct genes_t *transcript)
{
    int i;
    for(i=0;i<transcript->intervals_number;i++)
    {
    	if (snp->pos>=transcript->init_intervals[i] && snp->pos<=transcript->end_intervals[i])
	    	return(1);
    }
    return(0);
}


// Returns read count of reference allele
int getRef(struct snps_t *elem)
{
    if (strncmp(elem->ref,"A",1)==0)
        return(elem->A);
    else if (strncmp(elem->ref,"C",1)==0)
        return(elem->C);
    else if (strncmp(elem->ref,"G",1)==0)
        return(elem->G);
    else if (strncmp(elem->ref,"T",1)==0)
        return(elem->T);
    return(0);
}

int getSum(struct snps_t *elem)
{
    return(elem->A+elem->C+elem->G+elem->T);
}

int isRefSpecified(struct snps_t *elem)
{
    if (strncmp(elem->ref,".",1)==0 || strncmp(elem->ref,"N",1)==0)
        return(0);
    return(1);
}

int isAltSpecified(struct snps_t *elem)
{
    if (strncmp(elem->alt,".",1)==0 || strncmp(elem->alt,"N",1)==0)
        return(0);
    return(1);
}

// Returns read count of alternative allele
int getAlt(struct snps_t *elem)
{
    int max=0;
    char reference[2];
    if (strncmp(elem->alt,"A",1)==0)
        return(elem->A);
    else if (strncmp(elem->alt,"C",1)==0)
        return(elem->C);
    else if (strncmp(elem->alt,"G",1)==0)
        return(elem->G);
    else if (strncmp(elem->alt,"T",1)==0)
        return(elem->T);
        
    // Reference is specified
    if(strncmp(elem->ref,"A",1)==0)
    {
    	sprintf(elem->alt,"C");
    	max = elem->C;
    	if(max<elem->G) { max = elem->G; sprintf(elem->alt,"G"); } 
	if(max<elem->T) { max = elem->T; sprintf(elem->alt,"T"); } 
	if( (max==elem->C)+(max==elem->G)+(max==elem->T)>1 ) { sprintf(elem->alt,"N"); return(0); }
	return(max);
    }
    
    if(strncmp(elem->ref,"C",1)==0)
    {
    	sprintf(elem->alt,"A");
    	max = elem->A;
    	if(max<elem->G) { max = elem->G; sprintf(elem->alt,"G"); }
	if(max<elem->T) { max = elem->T; sprintf(elem->alt,"T"); }
	if( (max==elem->A)+(max==elem->G)+(max==elem->T)>1 ) { sprintf(elem->alt,"N"); return(0); }
	return(max);
    }
    
    if(strncmp(elem->ref,"G",1)==0)
    {
    	sprintf(elem->alt,"C");
    	max = elem->C;
    	if(max<elem->A) { max = elem->A; sprintf(elem->alt,"A"); }
	if(max<elem->T) { max = elem->T; sprintf(elem->alt,"T"); }
	if( (max==elem->C)+(max==elem->A)+(max==elem->T)>1 ) { sprintf(elem->alt,"N"); return(0); }
	return(max);
    }
    
    if(strncmp(elem->ref,"T",1)==0)
    {
    	sprintf(elem->alt,"C");
    	max = elem->C;
    	if(max<elem->G) { max = elem->G; sprintf(elem->alt,"G"); }
	if(max<elem->A) { max = elem->A; sprintf(elem->alt,"A"); }
	if( (max==elem->C)+(max==elem->G)+(max==elem->A)>1 ) { sprintf(elem->alt,"N"); return(0); }
	return(max);
    }
        
    return(0);
}

// Returns the chromosome position in the global chromosome indexing
int getChrPos(char *chr)
{
    if (strncmp(chr,"chr",3)==0)
    {
        char *tmp = chr;
        chr = (char *)malloc(strlen(chr)-2);
        strcpy(chr,tmp+3);
    }

    if (isInteger(chr)==1)
        return(-1);

    int val = atoi(chr);
    if (val>0 && val<=22)
        return(val-1);
    if (strncmp(chr,"X",1)==0)
        return(22);
    if (strncmp(chr,"Y",1)==0)
        return(23);
    if (strncmp(chr,"M",1)==0)
        return(24);
    return(-1);
}

void clearChrIndexing()
{
    int pos;
    for(pos=0;pos<25;pos++)
    {
        chr_indexing[pos] = -1;
    }
}

void resetVCF(struct snps_info *snps)
{
    int i;
    for(i=0;i<snps->length;i++)
    {
        snps->info[i]->A = snps->info[i]->C = snps->info[i]->G = snps->info[i]->T = 0;
        snps->info[i]->As = snps->info[i]->Cs = snps->info[i]->Gs = snps->info[i]->Ts = 0;
        snps->info[i]->ASE_pvalue = 0;
        snps->info[i]->FISHER_pvalue = 1;
        snps->info[i]->control_alt = -1;
        snps->info[i]->control_ref = -1;
	free(snps->info[i]->alt);
        snps->info[i]->alt =(char*)malloc(sizeof(char)*strlen(snps->info[i]->alt_default)+1);
        strcpy(snps->info[i]->alt,snps->info[i]->alt_default);
        free(snps->info[i]->ref);
        snps->info[i]->ref =(char*)malloc(sizeof(char)*strlen(snps->info[i]->ref_default)+1);
        strcpy(snps->info[i]->ref,snps->info[i]->ref_default);
    }
}

void clearVCF(struct snps_info *snps)
{
    int i;
    for(i=0;i<snps->length;i++)
    {
        free(snps->info[i]->alt);
        free(snps->info[i]->chr);
        free(snps->info[i]->dbsnp);
        free(snps->info[i]->maf);
        free(snps->info[i]->ref);
        free(snps->info[i]->alt_default);
        free(snps->info[i]->ref_default);
        free(snps->info[i]);
    }
    free(snps);
}

/////////////////////////////////////////////////////////////////////////////////////
// Load BED and SNP structures
/////////////////////////////////////////////////////////////////////////////////////

// Loads the gene/transcript data structure
struct genes_info* loadGenesTranscripts(char *file_name, int samples_number, int mode)
{
    FILE *file = fopen(file_name,"r");
    int i;

    if  (file==NULL)
    {
        fprintf(stderr,"\nFile %s not present!!!\n",file_name);
        exit(1);
    }

    // get number of lines
    int number_of_lines = 0;

    char line [MAX_READ_BUFF];
    while(fgets(line,sizeof(line),file) != NULL )
    {
       int i=0;
       while (isspace(line[i])) {i++;}
       if (line[i]!='#')
        number_of_lines++;
    }

    rewind(file);

    // create the overall structure
    struct genes_info *genes = (struct genes_info *)malloc(sizeof(struct genes_info));
    genes->info = malloc(sizeof(struct genes_t *)*number_of_lines);
    genes->length = number_of_lines;

    int index = 0;
    int control = 0;
    int columns = 0;
    char sep[] = "\t";
    char sep_intervals[] = ",";

    int line_numb = 1;

    while(fgets(line,sizeof(line),file) != NULL )
    {
        if (control == 0)
        {
            // count the number of columns
            for(i=0;i<strlen(line);i++)
            {
                if (strncmp(&line[i],sep,1)==0)
                    columns++;
            }
            columns++;
            control = 1;
        }

        // tokenize a line
        char *str_tokens[columns];
        char *pch;
        pch = strtok(line,"\t");
        i=0;
        while (pch != NULL)
        {
            str_tokens[i++]=pch;
            pch = strtok(NULL,"\t");
        }
        
        if(mode==1 && i<6)
        {
            fprintf(stderr,"ERROR: at line %d the number of columns is not correct. We require 6 columns (CHR,FROM,TO,ID,INIT_INTERVALS,END_INTERVALS).\n",line_numb);
            exit(1);
        }

        if (mode==0 && i<4)
        {
            fprintf(stderr,"ERROR: at line %d the number of columns is not correct. We require 4 columns (CHR,FROM,TO,ID).\n",line_numb);
            exit(1);
        }

        struct genes_t *current_elem = (struct genes_t *)malloc(sizeof(struct genes_t));
        int size;
        current_elem->chr = (char*)malloc(strlen(str_tokens[0])+1);
        strcpy(current_elem->chr,str_tokens[0]);
        current_elem->from = atoi(str_tokens[1]);
        current_elem->to = atoi(str_tokens[2]);
        if (str_tokens[3][strlen(str_tokens[3])-1]=='\n')
            size = strlen(str_tokens[3])-1;
        else
            size = strlen(str_tokens[3]);
        current_elem->gene = (char*)malloc(size+1);
        strncpy(current_elem->gene,str_tokens[3],size);
        current_elem->gene[size]='\0';
        
        current_elem->samples = (double *)malloc(sizeof(double)*samples_number);
        current_elem->samples_snps = (int *)malloc(sizeof(int)*samples_number);
        for(i=0;i<samples_number;i++)
        {
                current_elem->samples[i] = 0;
                current_elem->samples_snps[i] = 0;
        }
        
        ///////////////////////////////////
        // If transcripts add the intervals
        ///////////////////////////////////
        
        if(mode==1)
        {
        
		int inits=0;
		char *current = str_tokens[4];
		for(i=0;i<strlen(current);i++)
		{
			if (current[i]==',')
				inits++;
		}
		inits++;
	
		int ends=0;
		current = str_tokens[5];
		for(i=0;i<strlen(current);i++)
		{
			if (current[i]==',')
				ends++;
		}
		ends++;
	
		if(ends!=inits)
		{
		    fprintf(stderr,"ERROR: Intervals at line %d inits (%d) and end (%d) intervals does not match.\n",line_numb,inits,ends);
		    exit(1);
		}
		
		current_elem->intervals_number = inits;
		current_elem->init_intervals = malloc(inits*sizeof(int));
		current_elem->end_intervals = malloc(ends*sizeof(int));
		
		char *inits_tokens[inits];
	    	char *ends_tokens[ends];
		int x,y;
		
		pch = strtok(str_tokens[4],",");
		i=0;
		while (pch != NULL)
		{
		    inits_tokens[i++]=pch;
		    pch = strtok(NULL,",");
		}
		
		pch = strtok(str_tokens[5],",");
		i=0;
		while (pch != NULL)
		{
		    ends_tokens[i++]=pch;
		    pch = strtok(NULL,",");
		}
		
		for (i=0;i<inits;i++)
		{
			current_elem->init_intervals[i] = atoi(inits_tokens[i]);
			current_elem->end_intervals[i] = atoi(ends_tokens[i]);
		}
		
		// Check intervals consistency (should be monotonic) and in coordinates
		for (i=0;i<inits;i++)
		{
			if (current_elem->init_intervals[i]<current_elem->from || current_elem->init_intervals[i]>current_elem->to || 
			    current_elem->end_intervals[i]<current_elem->from || current_elem->end_intervals[i]>current_elem->to)
			{
				fprintf(stderr,"ERROR: Interval [%d,%d - %d,%d] at line %d outside transcript region.\n",line_numb,
					current_elem->from,current_elem->to,current_elem->init_intervals[i],current_elem->end_intervals[i]);
		    		exit(1);
			}
		}
		for (i=1;i<inits;i++)
		{
			if (current_elem->init_intervals[i] <= current_elem->init_intervals[i-1] || 
			    current_elem->end_intervals[i] <= current_elem->end_intervals[i-1] )
			{
				fprintf(stderr,"ERROR: Intervals at line %d are not ordered.\n",line_numb);
		    		exit(1);
			}
		}
        
        }
        
        ///////////////////////////////////
        
        genes->info[index] = current_elem;

        index++;
        line_numb++;
    }

    return(genes);
}

// Load the SNPs data structure from BED file
struct snps_info* loadSNPsBED(char *file_name, int granularity, char * mode)
{
    FILE *file = fopen(file_name,"r");
    int chr_pos,i,control,line_numb,number_of_lines,index,columns;

    clearChrIndexing();

    if  (file==NULL)
    {
        fprintf(stderr,"\nFile %s not present!!!\n",file_name);
        exit(1);
    }

    number_of_lines = 0;
    line_numb = 1;
    control = 0;
    columns = 0;
    char sep[] = "\t";
    char *current;

    char line [MAX_READ_BUFF];
    fgets(line,sizeof(line),file);
    if (hasInteger(line)==1)
        rewind(file);

    while(fgets(line,sizeof(line),file) != NULL )
    {
        if (control == 0)
        {
            // count the number of columns
            for(i=0;i<strlen(line);i++)
            {
                if (strncmp(&line[i],sep,1)==0)
                    columns++;
            }
            columns++;
            control = 1;
        }
        
        // tokenize a line
        char *str_tokens[columns];
        char *pch;
        pch = strtok(line,"\t");
        i=0;
        while (pch != NULL)
        {
            str_tokens[i++]=pch;
            pch = strtok(NULL,"\t");
        }

	//fprintf(stderr,"-- %d %d\n",atoi(str_tokens[2]),atoi(str_tokens[1]));

        if (i<3)
        {
            fprintf(stderr,"ERROR: at line %d the number of columns is not correct.\n",line_numb);
            exit(1);
        }

        number_of_lines += (atoi(str_tokens[2])-atoi(str_tokens[1]))/granularity;
        line_numb++;
    }

    fprintf(stderr," %d ",number_of_lines);

    // create the overall SNPs structure
    struct snps_info *snps = (struct snps_info *)malloc(sizeof(struct snps_info));
    snps->info = malloc(sizeof(struct snps_t *)*number_of_lines);
    snps->length = number_of_lines;

    rewind(file);
    fgets(line,sizeof(line),file);
    if (hasInteger(line)==1)
        rewind(file);
    index=0;

    while(fgets(line,sizeof(line),file) != NULL )
    {
        char *str_tokens[columns];
        char *pch;
        pch = strtok(line,"\t");
        
        i=0;
        while (pch != NULL)
        {
            str_tokens[i++]=pch;
            pch = strtok(NULL,"\t");
        }

	i=0;
        while(i<(atoi(str_tokens[2])-atoi(str_tokens[1]))/granularity)
        {
            struct snps_t *current_elem = (struct snps_t *)malloc(sizeof(struct snps_t));

            current_elem->chr =(char*)malloc(sizeof(char)*strlen(str_tokens[0])+1);
            strcpy(current_elem->chr,str_tokens[0]);
            current_elem->pos = atoi(str_tokens[1])+(i*granularity)+1; // +1 becase BED is 0-based
            current_elem->dbsnp = (char*)malloc(2*sizeof(char));

            if (columns==3)
            {
            	sprintf(current_elem->dbsnp,".");
            	
            } else
            	{
            		current_elem->dbsnp = (char*)malloc(sizeof(char)*strlen(str_tokens[2])+1);
            		strcpy(current_elem->dbsnp,str_tokens[3]);
            		if (current_elem->dbsnp[strlen(current_elem->dbsnp)-1]=='\n')
            		{
            			current_elem->dbsnp[strlen(current_elem->dbsnp)-1] = '\0';
            		}
            	}
            current_elem->ref = (char*)malloc(2*sizeof(char));
            sprintf(current_elem->ref,".");
            current_elem->ref_default = (char*)malloc(2*sizeof(char));
            sprintf(current_elem->ref_default,".");
            current_elem->alt = (char*)malloc(2*sizeof(char));
            sprintf(current_elem->alt,".");
            current_elem->alt_default = (char*)malloc(2*sizeof(char));
            sprintf(current_elem->alt_default,".");
            current_elem->maf = (char*)malloc(2*sizeof(char));
            sprintf(current_elem->maf,".");
            current_elem->A = current_elem->C = current_elem->G = current_elem->T = 0;
            current_elem->As = current_elem->Cs = current_elem->Gs = current_elem->Ts = 0;
	        current_elem->ASE_pvalue = 0;
            current_elem->FISHER_pvalue = 1;
            current_elem->control_alt = -1;
            current_elem->control_ref = -1;
            snps->info[index] = current_elem;

            chr_pos = getChrPos(current_elem->chr);
            if(chr_indexing[chr_pos]>0 && current_elem->pos<snps->info[chr_indexing[chr_pos]]->pos && strncmp(mode,"ASE",3) == 0)
            {
            	fprintf(stderr,"\nERROR: In ASE mode VCF input file should be ordered by chromosome and position (%i,%i)!!\n",current_elem->pos,chr_indexing[chr_pos]);
        	exit(1);
            }
            if (chr_pos>=0 && chr_indexing[chr_pos]<0)
            {
                chr_indexing[chr_pos] = index;
                if(chr_pos<24 && strncmp(mode,"ASE",3) == 0)
                {
                	int k=chr_pos+1;
                	while(k<25)
                	{
                		if(chr_indexing[k]>0)
                		{
                			fprintf(stderr,"ERROR: In ASE mode VCF input file should be ordered by chromosome and position!\n");
        				exit(1);
                		}
                		k++;
                	}
                }
	    }
            index++;
            i++;
        }
    }
    fclose (file);

    return snps;
}

// Load the SNPs data structure
struct snps_info* loadSNPsList(char *file_name, char *mode)
{
    FILE *file = fopen(file_name,"r");
    int chr_pos,i;

    if  (file==NULL)
    {
        fprintf(stderr,"\nFile %s not present!!!\n",file_name);
        exit(1);
    }

    // initialize chr_indexing
    clearChrIndexing();

    // get number of lines
    int number_of_lines = 0;

    char line [MAX_READ_BUFF];
    fgets(line,sizeof(line),file);
    if (hasInteger(line)==1)
        rewind(file);
    while(fgets(line,sizeof(line),file) != NULL )
    {
       i=0;
       while (isspace(line[i])) {i++;}
       if (line[i]!='#')
        number_of_lines++;
    }

    //fprintf(stderr," %d ",number_of_lines);

    rewind(file);

    // create the overall SNPs structure
    struct snps_info *snps = (struct snps_info *)malloc(sizeof(struct snps_info));
    snps->info = malloc(sizeof(struct snps_t *)*number_of_lines);
    snps->length = number_of_lines;

    int control = 0;
    int columns = 0;

    // initialize structure
    int index= 0;
    char sep[] = "\t";
    int line_numb = 1;

    fgets(line,sizeof(line),file);
    if (hasInteger(line)==1)
        rewind(file);
    while(fgets(line,sizeof(line),file) != NULL )
    {
        int i=0;
        while (isspace(line[i])) {i++;}
        if (line[i]!='#')
        {
            if (control == 0)
            {
                // count the number of columns
                for(i=0;i<strlen(line);i++)
                {
                    if (strncmp(&line[i],sep,1)==0)
                        columns++;
                }
                columns++;
                control = 1;
            }

            // tokenize a line
            char *str_tokens[columns];
            char *pch;
            pch = strtok(line,"\t");
            i=0;
            while (pch != NULL)
            {
                str_tokens[i++]=pch;
                pch = strtok (NULL, "\t");
            }

            if (i<8)
            {
                fprintf(stderr,"ERROR: at line %d the number of columns is not correct. VCF format requires at least 8 columns.\n",line_numb);
                exit(1);
            }

            struct snps_t *current_elem = (struct snps_t *)malloc(sizeof(struct snps_t));

            // extract MAF
            pch = strtok(str_tokens[7],";");
            current_elem->maf = NULL;
            current_elem->control_alt = -1;
            current_elem->control_ref = -1;

            while (pch != NULL)
            {
                if (strncmp(pch,"GMAF=",5)==0)
                {
                   current_elem->maf = (char*)malloc(sizeof(char)*(strlen(pch)-5)+1);
                   strncpy(current_elem->maf,pch+5,strlen(pch)-5);
                   current_elem->maf[strlen(pch)-5] = '\0';
                }
                if (strncmp(pch,"ASEQ_REF=",9)==0)
                {
                   current_elem->control_ref = atoi(pch+9);
                }
                if (strncmp(pch,"ASEQ_ALT=",9)==0)
                {
                   current_elem->control_alt = atoi(pch+9);
                }
                pch = strtok(NULL, ";");
            }

            if(current_elem->maf==NULL)
            {
                current_elem->maf = (char*)malloc(sizeof(char)+1);
                strcpy(current_elem->maf,".");
            }
            else
            {
            	if (current_elem->maf[strlen(current_elem->maf)-1]=='\n')
            	{
            		current_elem->maf[strlen(current_elem->maf)-1] = '\0';
            	}
            }

            current_elem->chr =(char*)malloc(sizeof(char)*strlen(str_tokens[0])+1);
            strcpy(current_elem->chr,str_tokens[0]);
            current_elem->pos = atoi(str_tokens[1]);
            current_elem->dbsnp = (char*)malloc(sizeof(char)*strlen(str_tokens[2])+1);
            strcpy(current_elem->dbsnp,str_tokens[2]);
            current_elem->ref = (char*)malloc(sizeof(char)*strlen(str_tokens[3])+1);
            strcpy(current_elem->ref,str_tokens[3]);
            current_elem->ref_default = (char*)malloc(sizeof(char)*strlen(str_tokens[3])+1);
            strcpy(current_elem->ref_default,str_tokens[3]);
            current_elem->alt = (char*)malloc(sizeof(char)*strlen(str_tokens[4])+1);
            strcpy(current_elem->alt,str_tokens[4]);
            current_elem->alt_default = (char*)malloc(sizeof(char)*strlen(str_tokens[3])+1);
            strcpy(current_elem->alt_default,str_tokens[4]);
            current_elem->A = current_elem->C = current_elem->G = current_elem->T = 0;
            current_elem->As = current_elem->Cs = current_elem->Gs = current_elem->Ts = 0;
            current_elem->ASE_pvalue = 0;
            current_elem->FISHER_pvalue = 1;
            snps->info[index] = current_elem;

            //TODO: controllare che reference e alternative siano validi

            chr_pos = getChrPos(current_elem->chr);
            if(chr_indexing[chr_pos]>0 && current_elem->pos<snps->info[chr_indexing[chr_pos]]->pos && strncmp(mode,"ASE",3) == 0)
            {
            	fprintf(stderr,"\nERROR: In ASE mode VCF input file should be ordered by chromosome (1,2,3,...) and position (%s,%i,%i)!!\n",current_elem->chr,current_elem->pos,chr_indexing[chr_pos]);
        	exit(1);
            }
            if (chr_pos>=0 && chr_indexing[chr_pos]<0)
            {
                chr_indexing[chr_pos] = index;
                if(chr_pos<24 && strncmp(mode,"ASE",3) == 0)
                {
                	int k=chr_pos+1;
                	while(k<25)
                	{
                		if(chr_indexing[k]>0)
                		{
                			fprintf(stderr,"ERROR: In ASE mode VCF input file should be ordered by chromosome (1,2,3,...) and position!\n");
        				exit(1);
                		}
                		k++;
                	}
                }
	    }

            index++;
            line_numb++;
        }
    }
    fclose (file);

    //printChrIndexing();

    return snps;
}

/////////////////////////////////////////////////////////////////////////////////////
// pileup SAM functions implementation
/////////////////////////////////////////////////////////////////////////////////////

// Callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
        int base,skip=0,op,l;
        snpsstruct_t *info = (snpsstruct_t *)data;
        cstate_t s = g_cstate_null;
        s.end = bam_calend(&b->core, bam1_cigar(b))-1;

        bam_pileup1_t elem;
        elem.b = b;

        if (b->core.flag & BAM_DEF_MASK)
            return 0;

        resolve_cigar_mine(&elem,info->snps->info[info->index]->pos-1,&s,&skip);

        int pos_my = elem.qpos;
        const bam1_core_t *c = &b->core;
        unsigned char *qual  = bam1_qual(b);
        
        int32_t read_init;
        int32_t read_end;

        /*
        // used only for debug purposes
        int n=0;
        char *qseq = (char *) malloc(b->core.l_qseq+1);
        unsigned char *se = bam1_seq(b);

        for(n=0;n<(b->core.l_qseq);n++)
        {
           int v = bam1_seqi(se,n);
           qseq[n] = bam_nt16_rev_table[v];
           if (n==pos_my)
            printf(" ");
           printf("%c",qseq[n]);
            if (n==pos_my)
            printf(" ");
        }
        qseq[n] = '\0';
        printf(" - %d - %d - %d\n",b->core.tid,b->core.pos,info->snps->info[info->index]->pos-1);
        */


        if (!(qual[pos_my] < BASE_QUALITY || c->qual < READ_QUALITY || elem.is_del != 0 || elem.indel != 0 || elem.is_refskip != 0 || skip == 1))
        {
            base = bam1_seqi(bam1_seq(b),pos_my);
            incBase(info->snps->info[info->index],base);
            incBaseStrand(info->snps->info[info->index],base,bam1_strand(b));
        }

        //bam_plbuf_t *buf = (bam_plbuf_t*)data;
        //bam_plbuf_push(b, buf);

        return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
// thread function
/////////////////////////////////////////////////////////////////////////////////////

// Function to be passed to the thread that executes the pileup of a number of SNPs
void *PileUp(void *args)
{
    struct args_thread *foo = (struct args_thread *)args;

    int start = foo->start;
    int end = foo->end;
    struct snps_info *snps = foo->snps;
    char* bamfile = foo->bamfile;
    bam_index_t *idx = foo->idx;

    tmpstruct_t tmp;

    tmp.beg = 0; tmp.end = 0x7fffffff;
    tmp.in = samopen(bamfile, "rb", 0);
    if (tmp.in == 0) {
            fprintf(stderr, "Fail to open BAM file %s\n", bamfile);
            return NULL;
    }

    char s[200];
    char stmp[100];
    int i;
    int ref;
    //bam_plbuf_t *buf;
    snpsstruct_t info;

    for(i=start;i<=end;i++)
    {
         s[0] = '\0';
         sprintf(stmp,"%s",snps->info[i]->chr);strcat(s,stmp);strcat(s,":");
         sprintf(stmp,"%d",snps->info[i]->pos);strcat(s,stmp);strcat(s,"-");
         sprintf(stmp,"%d",snps->info[i]->pos);strcat(s,stmp);

         bam_parse_region(tmp.in->header,s,&ref,&tmp.beg,&tmp.end); // parse the region
         if (ref < 0) {
            fprintf(stderr,"Invalid region!!! %d %s:%d ~ %s ~ %d:%d-%d\n",i,snps->info[i]->chr,snps->info[i]->pos,snps->info[i]->dbsnp,ref,tmp.beg,tmp.end);
            continue;
         }

         info.index = i;
         info.snps = snps;
         info.tmp = &tmp;

         //buf = bam_plbuf_init(pileup_func, &info); // initialize pileup
         bam_fetch(tmp.in->x.bam,idx,ref,tmp.beg,tmp.end,&info,fetch_func);
         //bam_plbuf_push(0, buf); // finalize pileup
         //bam_plbuf_destroy(buf);
    }

    samclose(tmp.in);
    
    if(foo->map!=NULL)
    	foo->map->val=1;
    
    return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////
// printing utility functions
/////////////////////////////////////////////////////////////////////////////////////

void printChrIndexing()
{
    int pos;
    fprintf(stderr,"\nChromosome indexing:");
    for(pos=0;pos<25;pos++)
    {
        fprintf(stderr,"\n- chr %d: %d",pos+1,chr_indexing[pos]);
    }
    fprintf(stderr,"\n");
}

// print analysis
void printSNPsInfo(struct snps_info *snps,FILE *outfile)
{
    int i,a,r,cont,rd,cov;
    float af;
    fprintf(outfile,"chr\tpos\tdbsnp\tMAF\tref\talt\tA\tC\tG\tT\tRD\tArs\tCrs\tGrs\tTrs\n");
    //fprintf(outfile,"chr\tpos\tdbsnp\tMAF\tref\talt\tA\tC\tG\tT\tRD\tArs\tCrs\tGrs\tTrs\taf\tcov\n");
    for(i=0;i<snps->length;i++)
    {
    	af=cov=-1;
    	cont=0;
    	
    	if (isRefSpecified(snps->info[i])==1 && isAltSpecified(snps->info[i])==1)
    	{
    		a = getAlt(snps->info[i]);
        	r = getRef(snps->info[i]);
        	cov=a+r;
        	if(cov>0)
        	{
        		af=((float)a)/((float)cov);
        	} else
        	{
        		af=0;
        	}	
        	if (a+r<READ_NUMBER)
        		cont=1;
        } 
        else
        {
        	a = getSum(snps->info[i]);
        	if (a<READ_NUMBER)
        		cont=1;
        }
        
        if(cont==1)
        	continue;
        
        rd = snps->info[i]->A+snps->info[i]->C+snps->info[i]->G+snps->info[i]->T;
        fprintf(outfile,"%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
               snps->info[i]->chr,
               snps->info[i]->pos,
               snps->info[i]->dbsnp,
               snps->info[i]->maf,
               snps->info[i]->ref,
               snps->info[i]->alt,
               snps->info[i]->A,snps->info[i]->C,snps->info[i]->G,snps->info[i]->T,
               rd,snps->info[i]->As,snps->info[i]->Cs,snps->info[i]->Gs,snps->info[i]->Ts);
    }
}

void printSNPsInfoSTDOUT(struct snps_info *snps)
{
    int i,a,r,cont,rd;
    fprintf(stdout,"chr\tpos\tdbsnp\tMAF\tref\talt\tA\tC\tG\tT\tRD\tArs\tCrs\tGrs\tTrs\n");
    for(i=0;i<snps->length;i++)
    {
    	cont=0;
    	
    	if (isRefSpecified(snps->info[i])==1 && isAltSpecified(snps->info[i])==1)
    	{
    		a = getAlt(snps->info[i]);
        	r = getRef(snps->info[i]);
        	if (a+r<READ_NUMBER)
        		cont=1;
        } 
        else
        {
        	a = getSum(snps->info[i]);
        	if (a<READ_NUMBER)
        		cont=1;
        }
        
        if(cont==1)
        	continue;
        	
        rd = snps->info[i]->A+snps->info[i]->C+snps->info[i]->G+snps->info[i]->T;

        fprintf(stdout,"%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
               snps->info[i]->chr,
               snps->info[i]->pos,
               snps->info[i]->dbsnp,
               snps->info[i]->maf,
               snps->info[i]->ref,
               snps->info[i]->alt,
               snps->info[i]->A,snps->info[i]->C,snps->info[i]->G,snps->info[i]->T,
               rd,snps->info[i]->As,snps->info[i]->Cs,snps->info[i]->Gs,snps->info[i]->Ts);
    }
}


void printHelp()
{
    fprintf(stderr,"Usage: \n mode PILEUP ==> ./computeCounts [vcf=string] [bam=string] [threads=int] [mbq=int] [mrq=int] [mdc=int] [out=string]\n");
}

int main(int argc, char **argv)
{   
   char *tmp_string = NULL;
        threads_alloc_list *curr, *head, *tail, *map, *map_prev, *tmp;
        
        int end_control = 0;
        int offset = 0;

        //printf("This is my AmpliSolve preprocessing code wrapped around ASEQ\n");
        //fprintf(stderr, "ASEQ version modified %d.%d.%d\n\n",VERSION1,VERSION2,VERSION3);

        //modify this code based on my pre-processing plan
        if (argc == 1)
        {
            printHelp();
            return 1;
        }

        //keep this part of code and modify it according to my pre-processing plan
        //fprintf(stderr,"1. Checking input parameters...");
        struct input_args *arguments;
        //arguments = getInputArgs(my_vcf,my_bam,my_threads,my_mbq,my_mrq,my_mdc,my_out);
        arguments = getInputArgs(argc,argv);
        int check = checkInputArgs(arguments);
        if( check==1 )
        {
              return 1;
        }
          
        //fprintf(stderr,"OK\n");

        
        // create outdir
        char *outdir = (char *)malloc(3);
        sprintf(outdir,"./");
        int control = 0;


        if (arguments->outdir != NULL)
        {
            if (isKeyword(arguments->outdir)==0)
            {
                free(outdir);
                outdir = arguments->outdir;
            }
            else
            {
                if (strncmp(arguments->mode,"PILEUP",6)==0 && arguments->vcffiles_number == 1 && arguments->bamfiles_number==1 && isKeyword(arguments->outdir)==1)
                {
                    PRINT_MODE = 1;
                }
                else
                {
                    //fprintf(stderr,"\nWARNING: STDOUT keyword cannot be used in this mode and will be ignored.\n");
                }
                    
            }
        }

        
        #ifdef _WIN32
        int result_code = mkdir(outdir);
        #else
        mode_t process_mask = umask(0);
        int result_code = mkdir(outdir, S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH);
        umask(process_mask);
        #endif

        // get current directory
        char cwd[1024];
        getcwd(cwd,sizeof(cwd));

        if (chdir(outdir)<0)
        {
            //fprintf(stderr,"ERROR from main: The out option is not a valid directory and is different from keyword STDOUT.\n");
            return(1);
        }
        chdir(&cwd);

        
        READ_QUALITY = arguments->mrq;
        BASE_QUALITY = arguments->mbq;
        READ_NUMBER = arguments->mrn;



        struct snps_info *snps = NULL;
        //struct genes_info *genes = NULL;

        // load SNPs list if only one
        if (arguments->vcf != NULL)
        {
            //fprintf(stderr,"2. Loading Panel VCF file %s...",arguments->vcf);
            snps  = loadSNPsList(arguments->vcf,arguments->mode);
            //fprintf(stderr,"OK\n");
        }
        
        int i,j,snps_per_core,snps_per_block;

        //fprintf(stderr,"2. Mode %s:\n\tthreads: %d\n\tmin base quality: %d\n\tmin read quality: %d\n\tmin depth of coverage: %d\n",arguments->mode,arguments->cores,arguments->mbq,arguments->mrq,arguments->mrn);
        
        
        if(arguments->pileup_mode == 0)
        {
            //fprintf(stderr,"\tpileup mode: static\n");
        }
            
        else
        {
            //fprintf(stderr,"\tpileup mode: dynamic (block size %d)\n",arguments->pileup_blocks);
        }
            

        FILE *outfile;
        char *outfile_name = NULL;
        int slash;
        char **names = malloc(sizeof(char *)*arguments->bamfiles_number);


        for(j=0;j<arguments->bamfiles_number;j++)
        {
                // cut name on .bam or if not .bam take as it is
                if (tmp_string != NULL)
                    free(tmp_string);
                slash = lastSlash(arguments->bamfiles[j])+1;
                if (strncmp(arguments->bamfiles[j]+(strlen(arguments->bamfiles[j])-4),".bam",4) == 0)
                {
                    tmp_string = (char *)malloc(strlen(arguments->bamfiles[j])-slash-3);
                    strncpy(tmp_string,arguments->bamfiles[j]+slash,strlen(arguments->bamfiles[j])-slash-4);
                    tmp_string[strlen(arguments->bamfiles[j])-slash-4] = '\0';
                }
                else
                {
                    tmp_string = (char *)malloc(strlen(arguments->bamfiles[j])-slash+1);
                    strcpy(tmp_string,arguments->bamfiles[j]+slash);
                }
                    names[j] = (char *)malloc(strlen(tmp_string)+1);
                    strcpy(names[j],tmp_string);
                if (arguments->vcf == NULL && arguments->bed == NULL)
                {
                    //fprintf(stderr,"- Loading VCF file %s...",arguments->vcffiles[j]);
                    snps  = loadSNPsList(arguments->vcffiles[j],arguments->mode);
                    //fprintf(stderr,"OK\n");
                }


                pthread_t threads[arguments->cores];
                struct args_thread args[arguments->cores];

                snps_per_core = ceil(snps->length/arguments->cores)+1;

                //fprintf(stderr,"3. Pileup computation of BAM file %s...",arguments->bamfiles[j]);
                i=0;
            
                    if(arguments->pileup_mode == 0)
                    {
                       // Fixed mode
                      while(i<arguments->cores)
                      {
                         args[i].start = i*snps_per_core;
                         args[i].end = (i+1)*snps_per_core-1;
                         args[i].idx = arguments->bamindexes[j];
                         args[i].snps = snps;
                         args[i].bamfile = arguments->bamfiles[j];
                         args[i].map = NULL;

                         if(args[i].end>=(snps->length-1))
                         args[i].end = snps->length-1;
                         //fprintf(stderr,"\n%d %d\n",args[i].start,args[i].end);
                         pthread_create(&threads[i],NULL,PileUp,(void*)(&args[i]));
                         i++;
                      }
                      for(i=0;i<arguments->cores;i++)
                      {
                         pthread_join(threads[i],NULL);
                         //fprintf(stderr,"*");
                      }
                    }
                    else
                    {
                       snps_per_block = arguments->pileup_blocks;
                       if(snps_per_block>snps_per_core)
                       snps_per_block = snps_per_core;
                    
                      // Dynamic mode
                      end_control = 0;
                      offset = 0;
                      head = NULL;
                      while(i<arguments->cores)
                      {
                         curr = (threads_alloc_list *)malloc(sizeof(threads_alloc_list));
                         curr->val = 0;
                         curr->thread = i;
                         curr->next = NULL;
                         if(head==NULL)
                         {
                             head = curr;
                         } 
                         else
                         {
                             curr->next = head;
                             head = curr;
                         }
                             args[i].start = offset*snps_per_block;
                             args[i].end = (offset+1)*snps_per_block-1;
                             args[i].idx = arguments->bamindexes[j];
                             args[i].snps = snps;
                             args[i].bamfile = arguments->bamfiles[j];
                             args[i].map = curr;
                             //fprintf(stderr,"\n%d %d\n",args[i].start,args[i].end);

                             if(args[i].end>=(snps->length-1))
                             {
                                 args[i].end = snps->length-1;
                                 end_control = 1;
                             }

                         pthread_create(&threads[i],NULL,PileUp,(void*)(&args[i]));
                         i++;   
                         offset++;  
                      }

                  int index = 0;
                  while(head!=NULL)
                  {
                         map = head;
                         map_prev = head;

                         while(map!=NULL)
                         {
                            if(map->val == 1)
                            {
                               // Delete current
                               i = map->thread;
                               if(map==head)
                               {
                                  tmp = map;
                                  head = map->next;
                                  map_prev = map = head;
                                  free(tmp);
                               }
                               else
                               {
                                  map_prev->next = map->next;
                                  tmp = map;
                                  map = map->next;
                                  free(tmp);
                               }
    
                               if(end_control==0)
                               {
                                 // Start new
                                curr = (threads_alloc_list *)malloc(sizeof(threads_alloc_list));
                                curr->val = 0;
                                curr->thread = i;
                                curr->next = NULL;

                                if(head==NULL)
                                {
                                   head = curr;
                                } 
                                else
                                {
                                   curr->next = head;
                                   head = curr;
                                }

                                  args[i].start = offset*snps_per_block;
                                  args[i].end = (offset+1)*snps_per_block-1;
                                  args[i].idx = arguments->bamindexes[j];
                                  args[i].snps = snps;
                                  args[i].bamfile = arguments->bamfiles[j];
                                  args[i].map = curr;

                                  if(args[i].end>=(snps->length-1))
                                  {
                                      args[i].end = snps->length-1;
                                      end_control = 1;
                                  }

                                  //fprintf(stderr,"\n%d %d\n",args[i].start,args[i].end);

                                  pthread_create(&threads[i],NULL,PileUp,(void*)(&args[i]));

                                  offset++;
                               }
                            } 
                            else
                            {
                             map_prev = map;
                             map = map->next;
                        }
                 }
            sleep(1);
        }
            
            }

            bam_index_destroy(arguments->bamindexes[j]);

            //fprintf(stderr,"...OK\n");

            chdir(outdir);

            if (PRINT_MODE == 0)
            {
                if(outfile_name!=NULL)
                free(outfile_name);
                outfile_name = (char*)malloc(strlen(tmp_string)+strlen(arguments->mode)+7);
                sprintf(outfile_name,"%s.%s.ASEQ",tmp_string,arguments->mode);

                fprintf(stderr,"\t\tcomputeCounts Output file %s in folder %s...",outfile_name,outdir);
                outfile = fopen(outfile_name,"w");
                printSNPsInfo(snps,outfile);
                fclose(outfile);
                fprintf(stderr,"OK\n");
            }
            else
            {
                printSNPsInfoSTDOUT(snps);
            }

            chdir(&cwd);

            if (arguments->vcf != NULL || arguments->bed != NULL)
                resetVCF(snps);
            else
                clearVCF(snps);
        }
        return 0;
    return 0;
}

