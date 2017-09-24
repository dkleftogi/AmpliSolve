#ifndef PARSEPARAMETERS_H_INCLUDED
#define PARSEPARAMETERS_H_INCLUDED

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "DataStructures.h"
#include "Utilities.h"

int MAX_READ_BUFF=1000000;

/////////////////////////////////////////////////////////////////////////////////////
// Parse arguments
/////////////////////////////////////////////////////////////////////////////////////

//char *my_vcf,char *my_bam,int my_threads,int my_mbq,int my_mrq,int my_mdc,char *my_out
struct input_args *getInputArgs(int argc, char **argv)
{
    int i;
    struct input_args *arguments = (struct input_args *)malloc(sizeof(struct input_args));
    arguments->cores = 1;
    arguments->mbq = 1;
    arguments->mrq = 1;
    arguments->mode = (char*)malloc(sizeof(char)*7);
    sprintf(arguments->mode,"PILEUP");
    arguments->vcf = NULL;
    arguments->bam = NULL;
    arguments->outdir=NULL;
    arguments->pileup_blocks = 1000;
    arguments->pileup_mode = 0;
    char *tmp=NULL;

    //arguments->vcf= (char*)malloc(sizeof(char)*strlen(my_vcf)+1);
    //sprintf(arguments->vcf,"%s",my_vcf);
    //subSlash(arguments->vcf);

    //arguments->bam = (char*)malloc(sizeof(char)*strlen(my_bam)+1);
    //sprintf(arguments->bam,"%s",my_bam);
    //subSlash(arguments->bam);

    //arguments->outdir = (char*)malloc(sizeof(char)*strlen(my_out)+1);
    //sprintf(arguments->outdir,"%s",my_out);
    //subSlash(arguments->outdir);

    //arguments->cores = my_threads;
    //arguments->mbq = my_mbq;
    //arguments->mrq = my_mrq;
    //arguments->mrn = my_mdc;


    for(i=1;i<argc;i++)
    {
        if( strncmp(argv[i],"mode=",5) == 0 )
        {
           free(arguments->mode);
           arguments->mode = (char*)malloc(sizeof(char)*strlen(argv[i])+1);
           strcpy(arguments->mode,argv[i]+5);
        }
        else if( strncmp(argv[i],"vcf=",4) == 0 )
        {
            arguments->vcf= (char*)malloc(sizeof(char)*strlen(argv[i])+1);
            strcpy(arguments->vcf,argv[i]+4);
			subSlash(arguments->vcf);
        }
        else if( strncmp(argv[i],"bam=",4) == 0 )
        {
            arguments->bam = (char*)malloc(sizeof(char)*strlen(argv[i])+1);
            strcpy(arguments->bam,argv[i]+4);
			subSlash(arguments->bam);
        }
        else if( strncmp(argv[i],"threads=",8) == 0 )
        {
            tmp = (char*)malloc(strlen(argv[i])-7);
            strcpy(tmp,argv[i]+8);
            arguments->cores = atoi(tmp);
            free(tmp);
        }
        else if( strncmp(argv[i],"mbq=",4) == 0 )
        {
            tmp = (char*)malloc(strlen(argv[i])-3);
            strcpy(tmp,argv[i]+4);
            arguments->mbq = atoi(tmp);
            free(tmp);
        }
        else if( strncmp(argv[i],"mrq=",4) == 0 )
        {
            tmp = (char*)malloc(strlen(argv[i])-3);
            strcpy(tmp,argv[i]+4);
            arguments->mrq = atoi(tmp);
            free(tmp);
        }
        else if( strncmp(argv[i],"mdc=",4) == 0 )
        {
            tmp = (char*)malloc(strlen(argv[i])-3);
            strcpy(tmp,argv[i]+4);
            arguments->mrn = atoi(tmp);
            free(tmp);
        }
        else if(strncmp(argv[i],"pilblock=",9) == 0 )
        {
            tmp = (char*)malloc(strlen(argv[i])-8);
            strcpy(tmp,argv[i]+9);
            arguments->pileup_blocks = atoi(tmp);
            free(tmp);
        }
        else if( strncmp(argv[i],"out=",4) == 0 )
        {
            arguments->outdir = (char*)malloc(strlen(argv[i])-3);
            strcpy(arguments->outdir,argv[i]+4);
			subSlash(arguments->outdir);
        }
        else
        {
            fprintf(stderr,"ERROR: input parameters not valid.\n");
            exit(1);
        }
    }

    return arguments;
}


int checkInputArgs(struct input_args *arguments)
{
    if( strncmp(arguments->mode,"PILEUP",6) != 0 )
    {
        fprintf(stderr,"ERROR checkInputArgs: mode not correct.\n");
        return 1;
    }

    if( strncmp(arguments->mode,"PILEUP",6) == 0 )
    {
        if ( arguments->vcf == NULL | arguments->bam == NULL )
        {
             fprintf(stderr,"ERROR from checkInputArgs function: vcf or bam files not present.\n");
             return 1;
        }
    } 
    if (arguments->cores <=0)
    {
        fprintf(stderr,"ERROR from checkInputArgs function: the number of threads is not valid.\n");
        return 1 ;
    }

    if (arguments->mbq <0)
    {
        fprintf(stderr,"ERROR from checkInputArgs function: the min base squality is not valid.\n");
        return 1 ;
    }

    if (arguments->mrq <0)
    {
        fprintf(stderr,"ERROR from checkInputArgs function: the mkin read quality is not valid.\n");
        return 1 ;
    }

    if (arguments->bam != NULL)
    {
        arguments->bamfiles = malloc(sizeof(char*));
        arguments->bamindexes = malloc(sizeof(bam_index_t *));
        arguments->bamfiles[0] = arguments->bam;

        if (FileExists(arguments->bamfiles[0])==1)
        {
            fprintf(stderr,"WARNING from checkInputArgs function: BAM file %s not present. \n",arguments->bamfiles[0]);
            return 1;
        }

        arguments->bamindexes[0] = bam_index_load(arguments->bamfiles[0]);
        if (arguments->bamindexes[0] == 0)
        {
                arguments->bamindexes[0] = bam_index_build(arguments->bamfiles[0]);
                fprintf(stderr,"WARNING from checkInputArgs function: bam indexing file not present. Will be created now...\n");
                arguments->bamindexes[0] = bam_index_load(arguments->bamfiles[0]);
        	       if (arguments->bamindexes[0] == 0)
        	       {
        		      fprintf(stderr,"ERROR from checkInputArgs function: unable to create index for BAM file.\n",arguments->bamfiles[0]);
        		      return 1;
        	       }
        }
        arguments->bamfiles_number = 1;

    }

    if (arguments->bam == NULL && arguments->bamlist != NULL)
    {
        FILE *file = fopen(arguments->bamlist,"r");

        if  (file==NULL)
        {
            fprintf(stderr,"ERROR from checkInputArgs function: BAM files list not present.\n");
            return 1;
        }

        arguments->bamfiles_number = 0;

        char line [MAX_READ_BUFF];
        while(fgets(line,sizeof(line),file) != NULL )
        {
            arguments->bamfiles_number++;
        }

        rewind(file);
        arguments->bamfiles = malloc(sizeof(char*)*arguments->bamfiles_number);
        arguments->bamindexes = malloc(sizeof(bam_index_t *)*arguments->bamfiles_number);
        int counter=0;
        int size=0;

        while(fgets(line,sizeof(line),file) != NULL )
        {
            if (line[strlen(line)-1]=='\n')
                size = strlen(line);
            else
                size = strlen(line)+1;
            arguments->bamfiles[counter] = (char*)malloc(size);
            strncpy(arguments->bamfiles[counter],line,size-1);
            (arguments->bamfiles[counter])[size-1] = '\0';
			subSlash(arguments->bamfiles[counter]);

            if (FileExists(arguments->bamfiles[counter])==1)
            {
                fprintf(stderr,"ERROR from checkInputArgs function: BAM file %s not present. \n",arguments->bamfiles[counter]);
                return 1;
            }

            // load BAM index
            arguments->bamindexes[counter] = bam_index_load(arguments->bamfiles[counter]);
            if (arguments->bamindexes[counter] == 0)
            {
                    arguments->bamindexes[counter] = bam_index_build(arguments->bamfiles[counter]);
                    fprintf(stderr,"WARNING from checkInputArgs function: BAM indexing of file %s not present. Will be created now.\n",arguments->bamfiles[counter]);
                    arguments->bamindexes[counter] = bam_index_load(arguments->bamfiles[counter]);
            	    if (arguments->bamindexes[counter] != 0)
            	    {
            	    	counter++;
            	    }
            	    else
            	    {
            	    	 fprintf(stderr,"WARNING from checkInputArgs function: unable to create BAM index for file %s.\n",arguments->bamfiles[counter]);
            	    }
            }
            else
            {
            	counter++;
            }
        }
        
        if (counter==0)
        {
        	 fprintf(stderr,"ERROR from checkInputArgs function: index creation failed for all BAM files.\n",arguments->bamfiles[counter]);
        	 return 1;
        }
        
    }

    if (arguments->vcf != NULL)
    {
        //arguments->vcffiles = malloc(sizeof(char*));
        //arguments->vcffiles[0] = arguments->vcf;

        //if (FileExists(arguments->vcffiles[0])==1)
        //{
        //   fprintf(stderr,"WARNING: vcf file %s not present. \n",arguments->vcffiles[0]);
        //    return 1;
        //}

        arguments->vcffiles_number = 1;

        //if (arguments->vcflist != NULL)
        //{
            //fprintf(stderr,"WARNING from checkInputArgs function: vcf files list will be ignored.\n");
        //}
        //if (arguments->bed != NULL)
        //{
            //fprintf(stderr,"WARNING here: bed file list will be ignored.\n");
        //}
    }
    
    if (arguments->vcf == NULL && arguments->bed != NULL)
    {
    	arguments->vcffiles_number = 1;

        if (arguments->vcflist != NULL)
        {
            fprintf(stderr,"WARNING from checkInputArgs function: vcf files list will be ignored.\n");
        }
    }

    if (arguments->vcf == NULL && arguments->vcflist != NULL && arguments->bed == NULL)
    {
        FILE *file = fopen(arguments->vcflist,"r");

        if  (file==NULL)
        {
            fprintf(stderr,"ERROR from checkInputArgs function: vcf files list %s not present.\n",arguments->vcflist);
            return 1;
        }

        arguments->vcffiles_number = 0;

        char line [MAX_READ_BUFF];
        while(fgets(line,sizeof(line),file) != NULL )
        {
            arguments->vcffiles_number++;
        }

        rewind(file);
        arguments->vcffiles = malloc(sizeof(char*)*arguments->vcffiles_number);
        int counter=0;
        int size = 0;

        while(fgets(line,sizeof(line),file) != NULL )
        {
            if (line[strlen(line)-1]=='\n')
                size = strlen(line);
            else
                size = strlen(line)+1;
            arguments->vcffiles[counter] = (char*)malloc(size);
            strncpy(arguments->vcffiles[counter],line,size-1);
            (arguments->vcffiles[counter])[size-1] = '\0';
			subSlash(arguments->vcffiles[counter]);

			if (FileExists(arguments->vcffiles[counter])==1)
            {
                fprintf(stderr,"WARNING from checkInputArgs function: vcf file %s not present. \n",arguments->vcffiles[counter]);
                return 1;
            }

            counter++;
        }
    }

    if (arguments->bamlist != NULL && arguments->vcflist != NULL && arguments->vcffiles_number != arguments->bamfiles_number)
    {
        fprintf(stderr,"ERROR from checkInputArgs function: lengths of vcf and bam files lists are not matching.\n");
        return 1;
    }

    return 0;
}


#endif // PARSEPARAMETERS_H_INCLUDED
