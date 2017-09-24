/******************************************************************************************************
									AmpliSolve: ctDNA analysis tool
BEGIN COPYRIGHT NOTICE

     AmpliSolve code -- (c) 2017 Dimitrios Kleftogiannis -- ICR -- www.icr.ac.uk

     Copyright 2017 Dimitrios Kleftogiannis Licensed under the
     Educational Community License, Version 2.0 (the "License"); you may
     not use this file except in compliance with the License. You may
     obtain a copy of the License at

     https://opensource.org/licenses/ECL-2.0

     Unless required by applicable law or agreed to in writing,
     software distributed under the License is distributed on an "AS IS"
     BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
     or implied. See the License for the specific language governing
     permissions and limitations under the License.

     Published reports of research using this code (or a modified version) should cite the 
     article that describes AmpliSolve tool and all of the programs that depends, see below. 
     
     Comments and bug reports are welcome.  
     Email to dimitrios.kleftogiannis@kicr.ac.uk 
     I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
     You are free to modify, extend or distribute this code, as long 
     as this copyright notice is included whole and unchanged. 

END COPYRIGHT NOTICE
 
UTILITY
  This program is a variant caller program for amplicon-based NGS experiments. It has been tested on Ion Torrent platform
  The program takes as input the arguments (see below for details) and produces a typical VCF and a flat file
  with more info about the variants.


INPUT ARGUMENTS
    1. pre-processing file          : to obtain this file run the pro-processing program AmpliSolvePreProVC.cpp
                              
    2. tumour bam file              : a bam file

    3. normal bam file            : a bam file that matches the tumour bam file

    4. threads                      : number of threads used for computing the read counts (same as in ASEQ)

    5. minimum base quality         : parameter for computing the read counts (default 20); better to use same setting as in AmpliSolvePreProVC.cpp

    6. minimum read quality         : parameter for computing the read counts (default 20)better to use same setting as in AmpliSolvePreProVC.cpp

    7. minimum coverage             : parameter for computing the read counts (default 20)better to use same setting as in AmpliSolvePreProVC.cpp

    8. output dir                   : directory for storing all results of execution
    

DEPENDENCIES
    The program does depend on the following programs:
    
    1. SAM tools http://samtools.sourceforge.net/
    
    2. Bed tools http://bedtools.readthedocs.io/en/latest/

    Make sure that both are intalled and also included in the PATH of your system.
    If not you may need a command like: PATH=$PATH:/your/path/to/Samtools
   
    3. ASEQ https://demichelislab.unitn.it/doku.php?id=public:aseq

    In fact the read count procedure has been wrapper around aseq software that utilizes
    internally samtools libraries. We acknowledge ASEQ implementation. 

    The makefile we provide compiles the codes and produces
    binary computeCounts that computes the read counts per position in the panel for a bam file.
    Program computCounts can be executed alone and incorportated into other pipelines. Please see our
    repository for an implementation that runs fast can can be executed in a cluster to speed up the process.
    
    4. Boost libraries.
    You have simply to download them from http://www.boost.org/ and give the path uppon compilation.
    In the boost web-site there are also more instructions and details.

    More details about AmpliSolve, descriptions and toy execution examples are provided in the Manual file. 

    Obviously this is a C++ program, thus you need the C++ compiler to be installed in your computer.
    The AmpliSolve program has been developed in a Mac OS computer with El Capitan version 10.11.5
    The C++ compiler used for development is the clang C, C++ and Objective-C compiler
    The program works for Unix-like systems.  

COMPILATION

Type the following if you want to compile only this code:

cc -I /Users/dkleftog/boost_1_61_0 -o AmpliSolveVC AmpliSolveVC.cpp -std=c++0x -lstdc++

If you want to install AmpliSolve tool type:

make -f XXXXX BOOST_FLAG=/Users/dkleftog/boost_1_61_0

and XXX is the makefile for your operating system

and BOOST_FLAG is the directory you have stored the boost libraries.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <time.h>
#include <ctime>
#include <sys/time.h>
#include <inttypes.h>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <map>
#include <iomanip>
#include <boost/math/distributions/hypergeometric.hpp>

//color the output
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[36m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"

//definitions of read buffers and other variables
#define MAXIMUM_READ_LENGTH 10000
#define MINIMUM_READ_LENGTH 1000
#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290

//functions
void printUsage();
void storeInputFile(char *file_name,char *dummyVCF);
void produceReadCounts_both(char *input_vcf, char *input_tumour_bam, char *input_germline_bam,char *input_threads,char *input_mbq,char *input_mrq, char *input_mdc, char *output_dir);
void produceReadCounts_tumour(char *input_vcf, char *input_tumour_bam,char *input_threads,char *input_mbq,char *input_mrq, char *input_mdc, char *output_dir);
void storeCountList(char *list_name, char *COUNT_DIR,std::unordered_map<std::string,std::string> &Hash);
void storeGermlineStatistics(std::unordered_map<std::string,std::string> &FILE_Hash,std::unordered_map<std::string,std::string> &Value_Hash_Info);
//functions for calling variants
double fisherTest(int a,int b,int c, int d);
double kf_lgamma(double z);
double kf_erfc(double x);
static double _kf_gammap(double s, double z);
static double _kf_gammaq(double s, double z);
double kf_gammap(double s, double z);
double kf_gammaq(double s, double z);
long double mutationRulesPoissonQualityScore(int supporting_reads,int RD,float AF_error);
void callVariants(std::unordered_map<std::string,std::string> &ReferenceBase, std::unordered_map<std::string,std::string> &DuplicatePosition,std::unordered_map<std::string,std::string> &TumourFileList,std::unordered_map<std::string,std::string> &Thresholds,char *output_dir);
void find_kmer_down(std::unordered_map<std::string,std::string> &Hash, char *chromosome, int position, char *output);
void find_kmer_up(std::unordered_map<std::string,std::string> &Hash, char *chromosome, int position,char *output);
int homopolymerTest(char *down, char *up,char sub);

//data structures
//store reference bases
std::unordered_map<std::string,std::string> ReferenceBase_Hash;
std::unordered_map<std::string,std::string>::iterator got_ReferenceBase_Hash;

//store duplicate positions 
std::unordered_map<std::string,std::string> DuplicatePosition_Hash;
std::unordered_map<std::string,std::string>::iterator got_DuplicatePosition_Hash;

std::unordered_map<std::string,std::string> Thresholds_Hash_Analytic;
std::unordered_map<std::string,std::string>::iterator got_Thresholds_Hash_Analytic;

//store the MAX AF for all germline positions
std::unordered_map<std::string,std::string> Germline_Max_Hash;
std::unordered_map<std::string,std::string>::iterator got_Germline_Max_Hash;

std::unordered_map<std::string,std::string> GermlineFileList_Hash;
std::unordered_map<std::string,std::string>::iterator got_GermlineFileList_Hash;

//storing the tumour files with complete paths and filenames
std::unordered_map<std::string,std::string> TumourFileList_Hash;
std::unordered_map<std::string,std::string>::iterator got_TumourFileList_Hash;

std::unordered_map<std::string,std::string> Germline_Hash_forPatients;
std::unordered_map<std::string,std::string>::iterator got_Germline_Hash_forPatients;

int main(int argc, char **argv)
{
        //initialize the input arguments, read the inputs from the command line and check for correctness
        char *user_prepro_file;
        user_prepro_file=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        char *user_tumour_bam;
        user_tumour_bam=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        char *user_germline_bam;
        user_germline_bam=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        char *user_threads;
        user_threads=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        char *user_mbq;
        user_mbq=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        char *user_mrq;
        user_mrq=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        char *user_mdc;
        user_mdc=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        char *user_output_dir;
        user_output_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        if(argc!=9)
        {
        std::cout<<"************************************************************************************************************************************"<<std::endl;
        std::cout<<ANSI_COLOR_RED<<"                                        Your input arguments are not correct"<<ANSI_COLOR_RESET<<std::endl;
        std::cout<<"                         Copyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n"<<std::endl;
        printUsage();
        return 0;
        }
        else
        {
          //store the input variable only if the number of input variables is correct
            //also this part makes some corrections for threads, mbq, mrq and mdc if not correct (i.e., when user gives <0)
            char *prepro_file;
            prepro_file=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            char *tumour_bam;
            tumour_bam=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            char *germline_bam;
            germline_bam=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            char *threads;
            threads=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            char *mbq;
            mbq=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            char *mrq;
            mrq=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            char *mdc;
            mdc=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            char *output_dir;
            output_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            strcpy(user_prepro_file,argv[1]);
            memset(prepro_file,0,MINIMUM_READ_LENGTH);
            sscanf(user_prepro_file,"preprocessing_file=%s",prepro_file);

            strcpy(user_tumour_bam,argv[2]);
            memset(tumour_bam,0,MINIMUM_READ_LENGTH);
            sscanf(user_tumour_bam,"tumour_bam=%s",tumour_bam);

            strcpy(user_germline_bam,argv[3]);
            memset(germline_bam,0,MINIMUM_READ_LENGTH);
            sscanf(user_germline_bam,"germline_bam=%s",germline_bam);

            strcpy(user_threads,argv[4]);
            memset(threads,0,MINIMUM_READ_LENGTH);
            sscanf(user_threads,"threads=%s",threads);

            strcpy(user_mbq,argv[5]);
            memset(mbq,0,MINIMUM_READ_LENGTH);
            sscanf(user_mbq,"mbq=%s",mbq);

            strcpy(user_mrq,argv[6]);
            memset(mrq,0,MINIMUM_READ_LENGTH);
            sscanf(user_mrq,"mrq=%s",mrq);

            strcpy(user_mdc,argv[7]);
            memset(mdc,0,MINIMUM_READ_LENGTH);
            sscanf(user_mdc,"mdc=%s",mdc);

            strcpy(user_output_dir,argv[8]);
            memset(output_dir,0,MINIMUM_READ_LENGTH);
            sscanf(user_output_dir,"output_dir=%s",output_dir);

            //convert values to int numbers
            int threads_int=std::atoi(threads);
            int mbq_int=std::atoi(mbq);
            int mrq_int=std::atoi(mrq);
            int mdc_int=std::atoi(mdc);

            //flag the germline bam values
            char flag_germline_bam[]="not_available";
            int no_germline=-1;
            std::cout<<"************************************************************************************************************************************\n"<<std::endl;
            std::cout<<"                          AmpliSolve: Variant Calling program for amplicon-based NGS data\n"<<std::endl;
            std::cout<<"                       Copyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n"<<std::endl;
            std::cout<<"Execution started under the following parameters:"<<std::endl;
            std::cout<<"\t1. Pre-processing file                : "<<ANSI_COLOR_GREEN<<prepro_file<<ANSI_COLOR_RESET<<std::endl;
            std::cout<<"\t2. Tumour bam                         : "<<ANSI_COLOR_GREEN<<tumour_bam<<ANSI_COLOR_RESET<<std::endl;
            
            if(strcmp(germline_bam,flag_germline_bam)==0)
            {
                no_germline=-1;
                std::cout<<"\t3. Germline bam                       : "<<ANSI_COLOR_RED<<"NO matched germline BAM available"<<ANSI_COLOR_RESET<<std::endl;
            }
            else
            {
                std::cout<<"\t3. Germline bam                       : "<<ANSI_COLOR_GREEN<<germline_bam<<ANSI_COLOR_RESET<<std::endl;
                no_germline=1;
            }
     
             if(threads_int<=0)
            {
                threads_int=2;
                std::cout<<"\t4. threads                            : "<<ANSI_COLOR_RED<<"User gave: "<<threads<<ANSI_COLOR_RESET<<" The value is converted to 2"<<std::endl;
            }
            else
            {
                std::cout<<"\t4. threads                            : "<<ANSI_COLOR_GREEN<<threads<<ANSI_COLOR_RESET<<std::endl;
            }
            if(mbq_int<=0)
            {
                mbq_int=20;
                std::cout<<"\t5. mbq                                : "<<ANSI_COLOR_RED<<"User gave: "<<mbq<<ANSI_COLOR_RESET<<". The value is converted to 20"<<std::endl;
            }
            else
            {
                std::cout<<"\t5. mbq                                : "<<ANSI_COLOR_GREEN<<mbq<<ANSI_COLOR_RESET<<std::endl;
            }
            
            if(mrq_int<=0)
            {
                mrq_int=20;
                std::cout<<"\t6. mrq                                    : "<<ANSI_COLOR_RED<<"User gave: "<<mrq<<ANSI_COLOR_RESET<<". The value is converted to 20"<<std::endl;
            }
            else
            {
                 std::cout<<"\t6. mrq                                : "<<ANSI_COLOR_GREEN<<mrq<<ANSI_COLOR_RESET<<std::endl;
            }

            if(mdc_int<=0)
            {
                mdc_int=20;
                std::cout<<"\t7. mdc                                : "<<ANSI_COLOR_RED<<"User gave: "<<mdc<<ANSI_COLOR_RESET<<". The value is converted to 20"<<std::endl;
            }
            else
            {
                 std::cout<<"\t7. mdc                                : "<<ANSI_COLOR_GREEN<<mdc<<ANSI_COLOR_RESET<<std::endl;
            }

            std::cout<<"\t8. output_dir                         : "<<ANSI_COLOR_GREEN<<output_dir<<ANSI_COLOR_RESET<<std::endl;

            if(no_germline==1)
            {
              //we have germlines so continue....

              char dummyVCF[50];
              memset(dummyVCF,0,50);
              srand (time(NULL));
              int t0=(rand()%100)+(rand()%100);
              sprintf(dummyVCF,"dummyVCF_%d.vcf",t0);


              std::cout<<"\nRunning function storeInputFile: ";
              storeInputFile(prepro_file,dummyVCF);

              std::cout<<"Running function produceReadCounts for both tumour and matched germline"<<std::endl;
              //produceReadCounts_both(dummyVCF,tumour_bam,germline_bam,threads,mbq,mrq,mdc,output_dir);

              //make a trick and get the actual file names of tumour and germline
              char tmp[100];
              
              //germline
              std::string str = germline_bam;
              std::string res = str.substr( str.find_last_of("/") + 1 );
              char filename_germline[100];
              memset(filename_germline,0,100);
              memset(tmp,0,100);
              sprintf(tmp,"%s",res.c_str());
              strncpy(filename_germline,tmp,strlen(tmp)-4);
              //tumour
              str = tumour_bam;
              res = str.substr( str.find_last_of("/") + 1 );
              char filename_tumour[100];
              memset(filename_tumour,0,100);
              memset(tmp,0,100);
              sprintf(tmp,"%s",res.c_str());
              strncpy(filename_tumour,tmp,strlen(tmp)-4);

              std::cout<<"\t\tTUMOUR file    : "<<output_dir<<"/"<<filename_tumour<<".PILEUP.ASEQ"<<std::endl;
              std::cout<<"\t\tGERMLINE file  : "<<output_dir<<"/"<<filename_germline<<".PILEUP.ASEQ"<<std::endl;

              std::ofstream output;

              char germline_count_list_name[50];
              memset(germline_count_list_name,0,50);
              srand (time(NULL));
              int t1=(rand()%100)+(rand()%100);
              sprintf(germline_count_list_name,"germline_count_list_%d.txt",t1);
              output.open(germline_count_list_name);
              output<<output_dir<<"/"<<filename_germline<<".PILEUP.ASEQ"<<std::endl;
              output.close();
              
              char tumour_count_list_name[50];
              memset(tumour_count_list_name,0,50);
              srand (time(NULL));
              int t2=(rand()%100)+(rand()%100);
              sprintf(tumour_count_list_name,"tumour_count_list_%d.txt",t2);
              output.open(tumour_count_list_name);
              output<<output_dir<<"/"<<filename_tumour<<".PILEUP.ASEQ"<<std::endl;
              output.close();

              //read and store the list of files for germline
              storeCountList(germline_count_list_name, output_dir,GermlineFileList_Hash);
              //do the same for tumour
              //this is just because we need to extend the processing for more samples in the future....
              storeCountList(tumour_count_list_name, output_dir,TumourFileList_Hash);
              //store the germlines 
              storeGermlineStatistics(GermlineFileList_Hash,Germline_Hash_forPatients);
              callVariants(ReferenceBase_Hash,DuplicatePosition_Hash,TumourFileList_Hash,Thresholds_Hash_Analytic,output_dir);   
            }
            else
            {
                //do the same as before but omit the parts that store the germline files...because there is no germline
                
                char dummyVCF[50];
                memset(dummyVCF,0,50);
                srand (time(NULL));
                int t0=(rand()%100)+(rand()%100);
                sprintf(dummyVCF,"dummyVCF_%d.vcf",t0);

                std::cout<<"\nRunning function storeInputFile: ";
                storeInputFile(prepro_file,dummyVCF);

                std::cout<<"Running function produceReadCounts for tumour sample "<<std::endl;
                produceReadCounts_tumour(dummyVCF,tumour_bam,threads,mbq,mrq,mdc,output_dir);

                //make a trick and get the actual file names of tumour and germline
                char tmp[100];
              
                //tumour
                std::string str = tumour_bam;
                std::string res = str.substr( str.find_last_of("/") + 1 );
                char filename_tumour[100];
                memset(filename_tumour,0,100);
                memset(tmp,0,100);
                sprintf(tmp,"%s",res.c_str());
                strncpy(filename_tumour,tmp,strlen(tmp)-4);

                std::cout<<"\t\tTUMOUR file    : "<<output_dir<<"/"<<filename_tumour<<".PILEUP.ASEQ"<<std::endl;
              
                std::ofstream output;
                
                char tumour_count_list_name[50];
                memset(tumour_count_list_name,0,50);
                srand (time(NULL));
                int t2=(rand()%100)+(rand()%100);
                sprintf(tumour_count_list_name,"tumour_count_list_%d.txt",t2);
                output.open(tumour_count_list_name);
                output<<output_dir<<"/"<<filename_tumour<<".PILEUP.ASEQ"<<std::endl;
                output.close();

                //just do this
                std::string tmp_1="test";
                Germline_Hash_forPatients.insert(std::make_pair(tmp_1,tmp_1));

                //do the same for tumour
                //this is just because we need to extend the processing for more samples in the future....
                storeCountList(tumour_count_list_name, output_dir,TumourFileList_Hash);
                
                callVariants(ReferenceBase_Hash,DuplicatePosition_Hash,TumourFileList_Hash,Thresholds_Hash_Analytic,output_dir); 
            }
        }
        TumourFileList_Hash.clear();
        GermlineFileList_Hash.clear();
        Germline_Hash_forPatients.clear();
        Thresholds_Hash_Analytic.clear();
        Germline_Max_Hash.clear();
        std::cout<<"\n************************************************************************************************************************************"<<std::endl;

}

//print program's usage
void printUsage()
{
  std::cout<<"Please type the following: "<<std::endl;
    std::cout<<"\n./AmpliSolveVC"<<ANSI_COLOR_GREEN<<" preprocessing_file="<<ANSI_COLOR_RESET<<"/your/file/returned/by/AmpliSolvePreProVC"<<ANSI_COLOR_GREEN<<" tumour_bam="<<ANSI_COLOR_RESET<<"/your/tumour/bam/file"<<ANSI_COLOR_GREEN<<" germline_bam="<<ANSI_COLOR_RESET<<"/your/germline/bam/file"<<ANSI_COLOR_GREEN<<" threads="<<ANSI_COLOR_RESET<<"number_of_threads"<<ANSI_COLOR_GREEN<<" mbq="<<ANSI_COLOR_RESET<<"minimum_base_quality"<<ANSI_COLOR_GREEN<<" mrq="<<ANSI_COLOR_RESET<<"minimum_read_quality"<<ANSI_COLOR_GREEN<<" mdc="<<ANSI_COLOR_RESET<<"minimum_coverage"<<ANSI_COLOR_GREEN<<" output_dir="<<ANSI_COLOR_RESET<<"/dir/to/store/all/outputs"<<std::endl;
    std::cout<<"\nRemember that: "<<std::endl;
    std::cout<<"This code internally uses a binary made by ASEQ codes that are based on SAMtools"<<std::endl;
    std::cout<<"If no matched germline bam file is available just type not_available in the germline_bam argument."<<std::endl;
    std::cout<<"\nExecution example:"<<std::endl;
    std::cout<<"The program takes 8 input arguments, so please fill them similar to the following command:"<<std::endl;
    std::cout<<"time ./AmpliSolveVC preprocessing_file=/Users/dkleftog/Downloads/aseq-v1.1.11-source/testOutput.txt tumour_bam=/Users/dkleftog/Desktop/tumour_data_samples_testing/CA36_plasma.bam germline_bam=/Users/dkleftog/Desktop/germline_data_samples_prepro_testing/CA36_germline_344.bam threads=4 mbq=20 mrq=20 mdc=20 output_dir=/Users/dkleftog/Desktop/AmpliSolveResults"<<std::endl;
    std::cout<<"************************************************************************************************************************************"<<std::endl;
}

//parse the input file with all info we need for variant calling 
void storeInputFile(char *file_name,char *dummyVCF)
{

        //contruct again the dummy VCF
        std::ofstream output;
        char dummy_BED_complete_name[MINIMUM_READ_LENGTH];
        sprintf(dummy_BED_complete_name,"%s",dummyVCF);
        output.open(dummy_BED_complete_name);

    char lineCharArray[MAXIMUM_READ_LENGTH];
    FILE *input;
    int ret=-1;
    int count=0;

    char flag[]="YES";

    if((input=fopen(file_name,"r"))==NULL)
    {
        printf("Error from storeInputFile function: Cannot open %s\n",file_name);
        exit(0);
    }
    else
    {
        //read the file line by line
        input=fopen(file_name,"r");
        //read the header
        memset(lineCharArray,0,MINIMUM_READ_LENGTH);
        ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
        while(!feof(input))
        {
            //read line by line
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            //then extract the positions
            char chrom[50];
            char position[50];
            char reference[50];
            char duplicate[50];
            char thres_A[50];
            char thres_C[50];
            char thres_G[50];
            char thres_T[50];
            char germ_max_A[50];
            char germ_max_C[50];
            char germ_max_G[50];
            char germ_max_T[50];

            if(ret!=EOF)
            {
                memset(chrom,0,50);
                memset(position,0,50);
                memset(reference,0,50);
                memset(duplicate,0,50);
                memset(thres_A,0,50);
                memset(thres_C,0,50);
                memset(thres_G,0,50);
                memset(thres_T,0,50);
                memset(germ_max_A,0,50);
                memset(germ_max_C,0,50);
                memset(germ_max_G,0,50);
                memset(germ_max_T,0,50);


                sscanf(lineCharArray,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",chrom,position,reference,duplicate,thres_A,thres_C,thres_G,thres_T,germ_max_A,germ_max_C,germ_max_G,germ_max_T);
                
                //if you want to process only one chromosome specified by the user, change 1 with strcmp(chromosome,chrom)==0 and add chromosome as argument
                if(1)
                {
                    //at this point we have the position we need and we are storing in the hash
                    //the key is the concatenation of chrom_position and the value is just the position
                    
                    //store the reference
                    char key[50];
                    memset(key,0,50);
                    sprintf(key,"%s_%s",chrom,position);
                    ReferenceBase_Hash.insert(std::make_pair<std::string,std::string>(key,reference));

                    //store the duplicates
                    if(strcmp(duplicate,flag)==0)
                    {

                      DuplicatePosition_Hash.insert(std::make_pair<std::string,std::string>(key,position));
                    }
                    else
                    {
                      //dont store it because it is not a duplicate
                    }

                    
                    char thres_key[50];
                    //store thresholds for A
                    memset(thres_key,0,50);
                    sprintf(thres_key,"%s_%s_A",chrom,position);
                    Thresholds_Hash_Analytic.insert(std::make_pair(thres_key,thres_A)); 

                    //store thresholds for A
                    memset(thres_key,0,50);
                    sprintf(thres_key,"%s_%s_C",chrom,position);
                    Thresholds_Hash_Analytic.insert(std::make_pair(thres_key,thres_C)); 

                    //store thresholds for A
                    memset(thres_key,0,50);
                    sprintf(thres_key,"%s_%s_G",chrom,position);
                    Thresholds_Hash_Analytic.insert(std::make_pair(thres_key,thres_G)); 

                    //store thresholds for A
                    memset(thres_key,0,50);
                    sprintf(thres_key,"%s_%s_T",chrom,position);
                    Thresholds_Hash_Analytic.insert(std::make_pair(thres_key,thres_T)); 


                    char germ_key[50];
                    //store max germline AF
                    memset(germ_key,0,50);
                    sprintf(germ_key,"%s_%s_A",chrom,position);
                    Germline_Max_Hash.insert(std::make_pair(germ_key,germ_max_A)); 

                    //store thresholds for A
                    memset(germ_key,0,50);
                    sprintf(germ_key,"%s_%s_C",chrom,position);
                    Germline_Max_Hash.insert(std::make_pair(germ_key,germ_max_C)); 

                    //store thresholds for A
                    memset(germ_key,0,50);
                    sprintf(germ_key,"%s_%s_G",chrom,position);
                    Germline_Max_Hash.insert(std::make_pair(germ_key,germ_max_G)); 

                    //store thresholds for A
                    memset(germ_key,0,50);
                    sprintf(germ_key,"%s_%s_T",chrom,position);
                    Germline_Max_Hash.insert(std::make_pair(germ_key,germ_max_T)); 

                    count++;

                    output<<chrom<<"\t"<<position<<"\t.\t.\t.\t.\t.\t."<<std::endl;
                }
            }
        }
    }
    fclose(input);
    std::cout<<" Running function storeInputFile: all info for panel reference bases stored with success "<<count<<ANSI_COLOR_RESET<<std::endl;
    //std::cout<<" \n\tReference Hash     : "<<ReferenceBase_Hash.size()<<std::endl;
    //std::cout<<" \n\tDuplicate Hash     : "<<DuplicatePosition_Hash.size()<<std::endl;
    //std::cout<<" \n\tThreshold Hash     : "<<Thresholds_Hash_Analytic.size()<<std::endl;
    //std::cout<<" \n\tGermline max Hash  : "<<Germline_Max_Hash.size()<<std::endl;
    output.close();
} 

void produceReadCounts_both(char *input_vcf, char *input_tumour_bam, char *input_germline_bam,char *input_threads,char *input_mbq,char *input_mrq, char *input_mdc, char *output_dir)
{

    char command[500];
    memset(command,0,500);
    int count=0;
    memset(command,0,500);
    sprintf(command,"./computeCounts vcf=%s bam=%s threads=%s mbq=%s mrq=%s mdc=%s out=%s",input_vcf,input_tumour_bam,input_threads,input_mbq,input_mrq,input_mdc,output_dir);
    system(command);
    memset(command,0,500);
    sprintf(command,"./computeCounts vcf=%s bam=%s threads=%s mbq=%s mrq=%s mdc=%s out=%s",input_vcf,input_germline_bam,input_threads,input_mbq,input_mrq,input_mdc,output_dir);      
    system(command);
}

void produceReadCounts_tumour(char *input_vcf, char *input_tumour_bam,char *input_threads,char *input_mbq,char *input_mrq, char *input_mdc, char *output_dir)
{
    char command[500];
    memset(command,0,500);
    int count=0;
    memset(command,0,500);
    sprintf(command,"./computeCounts vcf=%s bam=%s threads=%s mbq=%s mrq=%s mdc=%s out=%s",input_vcf,input_tumour_bam,input_threads,input_mbq,input_mrq,input_mdc,output_dir);
    system(command);
}

//store the list of files for processing
void storeCountList(char *list_name, char *COUNT_DIR,std::unordered_map<std::string,std::string> &Hash)
{
  //open the BAM_list file
  char lineCharArray[MINIMUM_READ_LENGTH];
  FILE *input;
  int ret=-1;
  int count=0;
  int size_DIR_name=strlen(COUNT_DIR);
  //check if exists
  if((input=fopen(list_name,"r"))==NULL)
  {
    printf("Error: Cannot open %s\n",list_name);
    exit(0);
  }
  else
  {
    //open the file
    input=fopen(list_name,"r");
    //read it all
    while(!feof(input))
    {
      memset(lineCharArray,0,MINIMUM_READ_LENGTH);
      ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
      
      if(ret!=EOF)
      {
        char file_name[MINIMUM_READ_LENGTH];
        memset(file_name,0,MINIMUM_READ_LENGTH);
        int size_current=-1;
        size_current=strlen(lineCharArray);
        //at this point be very careful...number 13 is based on the assumption that
        //the file names are PILEUP.ASEQ
        strncpy(file_name,lineCharArray+size_DIR_name+1,size_current-size_DIR_name-13);
        Hash.insert(std::make_pair<std::string,std::string>(lineCharArray,file_name));
        count++;
      }
    }
  }
  fclose(input);

  std::cout<<"\nRunning function storeList: "<<ANSI_COLOR_GREEN<< list_name <<ANSI_COLOR_RESET<<" stored with success. It contains "<<ANSI_COLOR_GREEN<<Hash.size()<<ANSI_COLOR_RESET<<" samples"<<std::endl;

  //for(std::unordered_map<std::string,std::string>::iterator it=Hash.begin(); it!=Hash.end();++it)
  //{
    //std::cout<<it->first<<" - "<<it->second<<std::endl;
  //}
}

//parse all germline files and store AFs from forward and reverse strand.
//this information is used for the threshold computation
void storeGermlineStatistics(std::unordered_map<std::string,std::string> &FILE_Hash,std::unordered_map<std::string,std::string> &Value_Hash_Info)
{    
    

    char chrom_chrX[]="chrX";
    char chrom_X[]="X";

    char chrom_chrY[]="chrY";
    char chrom_Y[]="Y";

    char chrom_chrM[]="chrM";
    char chrom_M[]="M";
    
    int ret=-1;
    char lineCharArray[MINIMUM_READ_LENGTH];
    //we parse one by one the original count files
     int count=0;
    for(std::unordered_map<std::string,std::string>::iterator it=FILE_Hash.begin(); it!=FILE_Hash.end();++it)
    {

        char filename[50];
        memset(filename,0,50);
        sprintf(filename,"%s",it->second.c_str());
        char patient_ID[50];
        memset(patient_ID,0,50);
        char file_type[50];
        memset(file_type,0,50);
        sscanf(filename,"%[^_]_%[^_]",patient_ID,file_type);
       
        count++;
        //current file that we read
        char current_file[MINIMUM_READ_LENGTH];
        memset(current_file,0,MINIMUM_READ_LENGTH);
        strcpy(current_file,it->first.c_str());
        //debug
        //cout<<current_file<<endl;
        FILE *input;

        if((input=fopen(current_file,"r"))==NULL)
        {
            printf("\tError from storeGermlineStatistics function: Cannot open %s\n",current_file);
            exit(0);
        }
        else
        { 
            //read the file line by line
            input=fopen(current_file,"r");
            //read the first line which is the header and then read it all till the end
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            while(!feof(input))
            {
                //read line by line
                memset(lineCharArray,0,MINIMUM_READ_LENGTH);
                ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
                
                //then get the ranges and extract the positions
                if(ret!=EOF)
                {
                    //define variables to read the input
                    char chrom[50];
                    char position[50];
                    char dbsnp[50];
                    char maf[50];
                    char ref[50];
                    char alt[50];
                    int base_A=-1;
                    int base_C=-1;
                    int base_G=-1;
                    int base_T=-1;
                    int RD=-1;
                    int base_Ar;
                    int base_Cr;
                    int base_Gr;
                    int base_Tr;
                    memset(chrom,0,50);
                    memset(position,0,50);
                    memset(dbsnp,0,50);
                    memset(maf,0,50);
                    memset(ref,0,50);
                    memset(alt,0,50);
                    
                    sscanf(lineCharArray,"%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",chrom,position,dbsnp,maf,ref,alt,&base_A,&base_C,&base_G,&base_T,&RD,&base_Ar,&base_Cr,&base_Gr,&base_Tr);
                    
                    if(1)
                    {
                            //reads in forward strands
                            int base_Afw=base_A-base_Ar;
                            int base_Cfw=base_C-base_Cr;
                            int base_Gfw=base_G-base_Gr;
                            int base_Tfw=base_T-base_Tr;

                            base_A=base_A;
                            base_C=base_C;
                            base_G=base_G;
                            base_T=base_T;
                            RD=RD;
                            base_Ar=base_Ar;
                            base_Cr=base_Cr;
                            base_Gr=base_Gr;
                            base_Tr=base_Tr;
                            base_Afw=base_Afw;
                            base_Cfw=base_Cfw;
                            base_Gfw=base_Gfw;
                            base_Tfw=base_Tfw;

                            //total reads in forward and backward starnds
                            int FW=(base_A-base_Ar)+(base_C-base_Cr)+(base_G-base_Gr)+(base_T-base_Tr);
                            int BW=base_Ar+base_Cr+base_Gr+base_Tr;

                            if((FW+BW)!=RD)
                            {
                                std::cout<<"malakia paizei edo"<<std::endl;
                            }

                            //AF in fw strands
                            float AF_Afw=0;
                            float AF_Cfw=0;
                            float AF_Gfw=0;
                            float AF_Tfw=0;

                            if(FW==0)
                            {
                                AF_Afw=0;
                                AF_Cfw=0;
                                AF_Gfw=0;
                                AF_Tfw=0;
                            }
                            else
                            {
                                AF_Afw=float(base_Afw)/float(FW);
                                AF_Cfw=float(base_Cfw)/float(FW);
                                AF_Gfw=float(base_Gfw)/float(FW);
                                AF_Tfw=float(base_Tfw)/float(FW);
                            }

                            //AF in bw strands
                            float AF_Abw=0;
                            float AF_Cbw=0;
                            float AF_Gbw=0;
                            float AF_Tbw=0;

                            if(BW==0)
                            {
                                AF_Abw=0;
                                AF_Cbw=0;
                                AF_Gbw=0;
                                AF_Tbw=0;
                            }
                            else
                            {
                                AF_Abw=float(base_Ar)/float(BW);
                                AF_Cbw=float(base_Cr)/float(BW);
                                AF_Gbw=float(base_Gr)/float(BW);
                                AF_Tbw=float(base_Tr)/float(BW);
                            }

                            float AF_A=0;
                            float AF_C=0;
                            float AF_G=0;
                            float AF_T=0;
                            AF_A=float(base_A)/float(RD);
                            AF_C=float(base_C)/float(RD);
                            AF_G=float(base_G)/float(RD);
                            AF_T=float(base_T)/float(RD);

                            //store each key to a separate Hash
                            //the key is always the concatenation of chrom and position : e.g., chrX_12413523
                            char value[MINIMUM_READ_LENGTH];
                            char key[MINIMUM_READ_LENGTH];
                    
                            //at this point I need to store the actual germline of the file using the patient ID
                            //and compute the max germline AF that is less than the cut-offs
                            memset(key,0,MINIMUM_READ_LENGTH);
                            sprintf(key,"%s_%s_A",chrom,position);
                            memset(value,0,MINIMUM_READ_LENGTH);
                            sprintf(value,"%f_%d_%d",AF_A,base_A,RD);
                            Value_Hash_Info.insert(std::make_pair(key,value));

                            //at this point I need to store the actual germline of the file using the patient ID
                            //and compute the max germline AF that is less than the cut-offs
                            memset(key,0,MINIMUM_READ_LENGTH);
                            sprintf(key,"%s_%s_C",chrom,position);
                            memset(value,0,MINIMUM_READ_LENGTH);
                            sprintf(value,"%f_%d_%d",AF_C,base_C,RD);
                            Value_Hash_Info.insert(std::make_pair(key,value));

                             //at this point I need to store the actual germline of the file using the patient ID
                            //and compute the max germline AF that is less than the cut-offs
                            memset(key,0,MINIMUM_READ_LENGTH);
                            sprintf(key,"%s_%s_G",chrom,position);
                            memset(value,0,MINIMUM_READ_LENGTH);
                            sprintf(value,"%f_%d_%d",AF_G,base_G,RD);
                            Value_Hash_Info.insert(std::make_pair(key,value));

                             //at this point I need to store the actual germline of the file using the patient ID
                            //and compute the max germline AF that is less than the cut-offs
                            memset(key,0,MINIMUM_READ_LENGTH);
                            sprintf(key,"%s_%s_T",chrom,position);
                            memset(value,0,MINIMUM_READ_LENGTH);
                            sprintf(value,"%f_%d_%d",AF_T,base_T,RD);
                            Value_Hash_Info.insert(std::make_pair(key,value)); 
                    }
                }
            }
        }
        fclose(input);
        //std::cout<<"\t\tDebug Parsing: "<<current_file<<" Storing:"<<Value_Hash_Info.size()<<std::endl;
    }
}

//this function computes the quality score given by the formula Q=-10*log10(pvalue), where p value is P(>=k) approximated by the Gamma function
long double mutationRulesPoissonQualityScore(int supporting_reads,int RD,float AF_error)
{
    long double Q=0;
    long double pvalue=0;
    long double p_limit=0.0000000001;
    long double m=0;
     
    //std::cout<<"Supporting Reads: "<<supporting_reads<<"\t"<<"Read Depth: "<<RD<<"\t"<<"AF_error rate: "<<AF_error<<"M value: "<<m<<std::endl;

    //need to do something with the NaNs
    if(AF_error == -1)
    {   
        Q=-888;
        //std::cout<<"Supporting Reads: "<<supporting_reads<<"\t"<<"Read Depth: "<<RD<<"\t"<<"AF_error rate: "<<AF_error<<"\t"<<"M value: "<<m<<"\t"<<"Q value"<<Q<<std::endl;
        return Q;
    }
    else
    {
        if(AF_error == 0)
        {
            //this is my magic error rate
            AF_error=0.0010008;
        }

        if(supporting_reads == 0)
        {
            pvalue=1;
        }
        else
        {
            m=double(RD)*AF_error;
            pvalue=1-kf_gammaq(supporting_reads,m);
        }

        if(pvalue < p_limit)
        {
            pvalue=p_limit;
            Q=-10*log10(pvalue);
        }
        else if( pvalue == 1)
        {   
            Q=0; 
        }
        else
        {
            Q=-10*log10(pvalue);
        }
        //std::cout<<"Supporting Reads: "<<supporting_reads<<"\t"<<"Read Depth: "<<RD<<"\t"<<"AF_error rate: "<<AF_error<<"M value: "<<m<<"\t"<<"Q value"<<Q<<std::endl;
        return Q;
    }
}

//implementation of Fisher exact test (two-tailed p value) using the hypergeometric distibutuion
double fisherTest(int a,int b,int c, int d)
{
    using namespace boost::math;
    unsigned N = a + b + c + d;
    unsigned r = a + c;
    unsigned n = c + d;
    unsigned max_for_k = std::min(r, n); 
    unsigned min_for_k = (unsigned)std::max(0, int(r + n - N));
    hypergeometric_distribution<> hgd(r, n, N); 
    double cutoff = pdf(hgd, c);
    double tmp_p = 0.0;
    for(int k = min_for_k;k < max_for_k + 1;k++) 
    {
        double p = pdf(hgd, k);
        if(p <= cutoff) tmp_p += p;
    }
    return tmp_p;
}

//function addopted by https://github.com/lh3/samtools/blob/master/bcftools/kfunc.c
double kf_lgamma(double z)
{
    double x = 0;
    x += 0.1659470187408462e-06 / (z+7);
    x += 0.9934937113930748e-05 / (z+6);
    x -= 0.1385710331296526     / (z+5);
    x += 12.50734324009056      / (z+4);
    x -= 176.6150291498386      / (z+3);
    x += 771.3234287757674      / (z+2);
    x -= 1259.139216722289      / (z+1);
    x += 676.5203681218835      / z;
    x += 0.9999999999995183;
    return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

//function addopted by https://github.com/lh3/samtools/blob/master/bcftools/kfunc.c
double kf_erfc(double x)
{
    const double p0 = 220.2068679123761;
    const double p1 = 221.2135961699311;
    const double p2 = 112.0792914978709;
    const double p3 = 33.912866078383;
    const double p4 = 6.37396220353165;
    const double p5 = .7003830644436881;
    const double p6 = .03526249659989109;
    const double q0 = 440.4137358247522;
    const double q1 = 793.8265125199484;
    const double q2 = 637.3336333788311;
    const double q3 = 296.5642487796737;
    const double q4 = 86.78073220294608;
    const double q5 = 16.06417757920695;
    const double q6 = 1.755667163182642;
    const double q7 = .08838834764831844;
    double expntl, z, p;
    z = fabs(x) * M_SQRT2;
    if (z > 37.) return x > 0.? 0. : 2.;
    expntl = exp(z * z * - .5);
    if (z < 10. / M_SQRT2) // for small z
        p = expntl * ((((((p6 * z + p5) * z + p4) * z + p3) * z + p2) * z + p1) * z + p0)
            / (((((((q7 * z + q6) * z + q5) * z + q4) * z + q3) * z + q2) * z + q1) * z + q0);
    else p = expntl / 2.506628274631001 / (z + 1. / (z + 2. / (z + 3. / (z + 4. / (z + .65)))));
    return x > 0.? 2. * p : 2. * (1. - p);
}

//function addopted by https://github.com/lh3/samtools/blob/master/bcftools/kfunc.c
// regularized lower incomplete gamma function, by series expansion
static double _kf_gammap(double s, double z)
{
    double sum, x;
    int k;
    for (k = 1, sum = x = 1.; k < 100; ++k) {
        sum += (x *= z / (s + k));
        if (x / sum < KF_GAMMA_EPS) break;
    }
    return exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
}

//function addopted by https://github.com/lh3/samtools/blob/master/bcftools/kfunc.c
// regularized upper incomplete gamma function, by continued fraction
static double _kf_gammaq(double s, double z)
{
    int j;
    double C, D, f;
    f = 1. + z - s; C = f; D = 0.;
    // Modified Lentz's algorithm for computing continued fraction
    // See Numerical Recipes in C, 2nd edition, section 5.2
    for (j = 1; j < 100; ++j) {
        double a = j * (s - j), b = (j<<1) + 1 + z - s, d;
        D = b + a * D;
        if (D < KF_TINY) D = KF_TINY;
        C = b + a / C;
        if (C < KF_TINY) C = KF_TINY;
        D = 1. / D;
        d = C * D;
        f *= d;
        if (fabs(d - 1.) < KF_GAMMA_EPS) break;
    }
    return exp(s * log(z) - z - kf_lgamma(s) - log(f));
}
//function addopted by https://github.com/lh3/samtools/blob/master/bcftools/kfunc.c
double kf_gammap(double s, double z)
{
    return z <= 1. || z < s? _kf_gammap(s, z) : 1. - _kf_gammaq(s, z);
}
//function addopted by https://github.com/lh3/samtools/blob/master/bcftools/kfunc.c
double kf_gammaq(double s, double z)
{
    return z <= 1. || z < s? 1. - _kf_gammap(s, z) : _kf_gammaq(s, z);
}

//this is the actual function that calls the variants
void callVariants(std::unordered_map<std::string,std::string> &ReferenceBase, std::unordered_map<std::string,std::string> &DuplicatePosition,std::unordered_map<std::string,std::string> &TumourFileList,std::unordered_map<std::string,std::string> &Thresholds,char *output_dir)
{
    std::cout<<"\nRunning function callVariants:"<<std::endl;

    char lineCharArray[MINIMUM_READ_LENGTH];
    int count=0;
    FILE *input;
    int ret=-1;
    int result=1;

    char case_A[]="A";
    char case_C[]="C";
    char case_G[]="G";
    char case_T[]="T";

    char Flag_Dup[50];
    double p=-1;
    char Flag_Fisher[50];

    char Flag_Tier[50];

    long double Q_fw=0;
    long double Q_bw=0;

    //iterators
    std::unordered_map<std::string,std::string>::iterator got_ReferenceBase;
    std::unordered_map<std::string,std::string>::iterator got_Thresholds;
    std::unordered_map<std::string,std::string>::iterator got_DuplicatePosition;
    
   

    for(std::unordered_map<std::string,std::string>::iterator it=TumourFileList.begin(); it!=TumourFileList.end();++it)
    {

      std::ofstream output;
      char output_filename[500];
      memset(output_filename,0,500);
      sprintf(output_filename,"%s/AmpliSolve_VC_Summary_%s.txt",output_dir,it->second.c_str());
      output.open(output_filename);
      //output<<"Filename\tChrom\tPosition\tSubtitution\tRD\tRD_fw\tRD_bw\tAF\tReads_fw\tReads_bw\tAF_fw\tAF_bw\tPoissionPvalue\tAmpliconEdge_StrandBias\tFisherPvalue\tTier"<<std::endl;
      //output<<"Filename\tChrom\tPosition\tSubtitution\tRD\tRD_fw\tRD_bw\tAF\tReads_fw\tReads_bw\tQ_fw\tQ_bw"<<std::endl;
      output<<"Filename\tChrom\tPosition\tSubtitution\tRD\tRD_fw\tRD_bw\tAF\tReads_fw\tReads_bw\tAF_fw\tAF_bw\tAmpliconEdge_StrandBias\tFisherPvalue\tQscore_fw\tQscore_bw\tReadTier\tGermlineInfo\tMaxGermlineAF\t10merDownstream\t10merUpstream\tHomopolymerFlag"<<std::endl;
       
        std::ofstream vcf_output;
        char vcf_file[100];
        memset(vcf_file,0,100);
        sprintf(vcf_file,"%s/%s.vcf",output_dir,it->second.c_str());
        char cat_dup_fisher[50];
        double max_germ=0;
        //print the date and time
        time_t now = time(0);
        // convert now to string form
        char* dt = ctime(&now);

        vcf_output.open(vcf_file);
        vcf_output<<"##fileformat=VCFv4.0\n##fileDate="<<dt<<"##source=AmpliSolve_VariantCalling\n##reference=??\n##phasing=??\n##INFO=<ID=RD,Number=1,Type=Integer,Description='Total Read Depth'>\n##INFO=<ID=AF,Number=.,Type=Float,Description='Allele Frequency'>\n##INFO=<ID=SR,Number=1,Type=String,Description='Supporting Reads'>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"<<std::endl;


        char filename_ID[150];
        char patient_ID[50];
        char sample_ID[50];
        memset(filename_ID,0,150);
        memset(sample_ID,0,50);
        memset(patient_ID,0,50);
        sprintf(filename_ID,"%s",it->second.c_str());
        sscanf(filename_ID,"%[^_]_%[^_]",patient_ID,sample_ID);

        char read_file[MINIMUM_READ_LENGTH];
        memset(read_file,0,MINIMUM_READ_LENGTH);
        sprintf(read_file,"%s",it->first.c_str());
        if((input=fopen(read_file,"r"))==NULL)
        {
          printf("\tError from callVariantsAnalytic:  Cannot open %s\n",read_file);
          exit(0);
        }
        else
        {
          //open the file
          input=fopen(read_file,"r");
          count++;
          //std::cout<<"\t\tParsing: "<<read_file<<std::endl;

          if(count % 1==0)
          {
          std::cout<<"\tParsed successfully "<<ANSI_COLOR_GREEN<<count<<"/"<<TumourFileList.size()<<ANSI_COLOR_RESET<<" plasma samples"<<std::endl;
          }


          //read the first line which is the header
          memset(lineCharArray,0,MINIMUM_READ_LENGTH);
          ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
          while(!feof(input))
          {
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            
            if(ret!=EOF)
            {
                    //define variables to read the input
                    char chrom[50];
                    char position[50];
                    char dbsnp[50];
                    char maf[50];
                    char ref[50];
                    char alt[50];
                    int base_A=-1;
                    int base_C=-1;
                    int base_G=-1;
                    int base_T=-1;
                    int RD=-1;
                    int base_Ar=-1;
                    int base_Cr=-1;
                    int base_Gr=-1;
                    int base_Tr=-1;
                    memset(chrom,0,50);
                    memset(position,0,50);
                    memset(dbsnp,0,50);
                    memset(maf,0,50);
                    memset(ref,0,50);
                    memset(alt,0,50);
                    sscanf(lineCharArray,"%s %s %s %s %s %s %d %d %d %d %d %d %d %d %d",chrom,position,dbsnp,maf,ref,alt,&base_A,&base_C,&base_G,&base_T,&RD,&base_Ar,&base_Cr,&base_Gr,&base_Tr);

                    if(1)
                    {

                        //Compute info as IS

                    //total reads in forward and backward starnds
                    int FW=(base_A-base_Ar)+(base_C-base_Cr)+(base_G-base_Gr)+(base_T-base_Tr);
                    int BW=base_Ar+base_Cr+base_Gr+base_Tr;
                    if((FW+BW)!=RD)
                    {
                        std::cout<<"malakia paizei edo"<<std::endl;
                    }
                    //reads in forward strands
                    int base_Afw=base_A-base_Ar;
                    int base_Cfw=base_C-base_Cr;
                    int base_Gfw=base_G-base_Gr;
                    int base_Tfw=base_T-base_Tr;
                    //AF in fw strands
                    float AF_Afw=0;
                    float AF_Cfw=0;
                    float AF_Gfw=0;
                    float AF_Tfw=0;
                    if(FW==0)
                    {
                        AF_Afw=0;
                        AF_Cfw=0;
                        AF_Gfw=0;
                        AF_Tfw=0;
                    }
                    else
                    {
                        AF_Afw=float(base_Afw)/float(FW);
                        AF_Cfw=float(base_Cfw)/float(FW);
                        AF_Gfw=float(base_Gfw)/float(FW);
                        AF_Tfw=float(base_Tfw)/float(FW);
                    }
                    //AF in bw strands
                    float AF_Abw=0;
                    float AF_Cbw=0;
                    float AF_Gbw=0;
                    float AF_Tbw=0;
                    if(BW==0)
                    {
                        AF_Abw=0;
                        AF_Cbw=0;
                        AF_Gbw=0;
                        AF_Tbw=0;
                    }
                    else
                    {
                        AF_Abw=float(base_Ar)/float(BW);
                        AF_Cbw=float(base_Cr)/float(BW);
                        AF_Gbw=float(base_Gr)/float(BW);
                        AF_Tbw=float(base_Tr)/float(BW);
                    }
                    //total freq
                    float AF_A=0;
                    float AF_C=0;
                    float AF_G=0;
                    float AF_T=0;
                    AF_A=float(base_A)/float(RD);
                    AF_C=float(base_C)/float(RD);
                    AF_G=float(base_G)/float(RD);
                    AF_T=float(base_T)/float(RD);
                    //extra variable needed for Fisher Exact test
                    int RD_reverse=base_Tr+base_Gr+base_Cr+base_Ar;

                    //use the pseduo counts for comparison with the germline cutoffs and then
                    //when reporting results use the normal values...those without extra count

                    //first find the reference base
                    char ref_key[50];
                    memset(ref_key,0,50);

                    sprintf(ref_key,"%s_%s",chrom,position);
                    
                    got_ReferenceBase=ReferenceBase.find(ref_key);
                    got_DuplicatePosition=DuplicatePosition.find(ref_key);

                    memset(Flag_Dup,0,50);
                    if(got_DuplicatePosition==DuplicatePosition.end())
                    {
                        //no duplicate position
                        sprintf(Flag_Dup,"NO");
                    }
                    else
                    {
                        sprintf(Flag_Dup,"YES");
                    }

                    memset(Flag_Fisher,0,50);
                    memset(Flag_Tier,0,50);
                    
                    if(got_ReferenceBase==ReferenceBase.end())
                    {
                        //just check that the code works...and all the reference base positions are inside the Hash
                        //debug
                        std::cout<<"malakia"<<std::endl;    
                    }
                    else
                    {
                        char reference[MINIMUM_READ_LENGTH];
                        memset(reference,0,MINIMUM_READ_LENGTH);
                        strcpy(reference,got_ReferenceBase->second.c_str());

                        char prompt_key[50];
                        char threshold_values[50];

                        char AF_fw_char[50];
                        char AF_bw_char[50];

                        float AF_fw=-1;;
                        float AF_bw=-1;;
                       

                        if(strcmp(reference,case_A)==0)
                        {
                          //reference is A so try C,G,T

                          //C
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_C",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);
                                
                                //results --> need to compute the poisson values and report flag and pvalue
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_C-base_Cr,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Cr,RD_reverse,AF_bw);
                                //Q score we set...
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_C-base_Cr,base_Cr);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Cfw<5 || base_Cr<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_C",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_C",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                        sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);

                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    
                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_C>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'C')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;

                                                  }
                                                  else
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                    }

                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"A->C"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_C<<"\t"<<base_Cfw<<"\t"<<base_Cr<<"\t"<<AF_Cfw<<"\t"<<AF_Cbw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'C')<<std::endl;
                                    //write the variant
                                    //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"A->C"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_C<<"\t"<<base_Cfw<<"\t"<<base_Cr<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;
                                }
                            }

                            //G
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_G",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);

                                //result
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_G-base_Gr,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Gr,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_G-base_Gr,base_Gr);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Gfw<5 || base_Gr<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_G",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_G",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                       sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;
                                    
                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_G>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'G')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;

                                                  }
                                                  else
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                    }

                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"A->G"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_G<<"\t"<<base_Gfw<<"\t"<<base_Gr<<"\t"<<AF_Gfw<<"\t"<<AF_Gbw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'G')<<std::endl;
                                    //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"A->G"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_G<<"\t"<<base_Gfw<<"\t"<<base_Gr<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }

                            //T
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_T",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);

                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_T-base_Tr,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Tr,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_T-base_Tr,base_Tr);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Tfw<5 || base_Tr<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_T",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_T",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                        sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;
                                    
                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_T>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'T')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;

                                                  }
                                                  else
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"A"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                    }

                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"A->T"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_T<<"\t"<<base_Tfw<<"\t"<<base_Tr<<"\t"<<AF_Tfw<<"\t"<<AF_Tbw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'T')<<std::endl;
                                    //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"A->T"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_T<<"\t"<<base_Tfw<<"\t"<<base_Tr<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }   
                    
                        }//end of ref A
                        else if(strcmp(reference,case_C)==0)
                        {
                            //reference is C

                            //A
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_A",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);


                                //results --> 
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_A-base_Ar,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Ar,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_A-base_Ar,base_Ar);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Afw<5 || base_Ar<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_A",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_A",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                       sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_A>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'A')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                       vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;

                                                  }
                                                  else
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                    }

                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"C->A"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_A<<"\t"<<base_Afw<<"\t"<<base_Ar<<"\t"<<AF_Afw<<"\t"<<AF_Abw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'A')<<std::endl;

                                }
                            }

                              //G
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_G",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);
                                
                                //results ->
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_G-base_Gr,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Gr,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_G-base_Gr,base_Gr);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Gfw<5 || base_Gr<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_G",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_G",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                        sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_G>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'G')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;

                                                  }
                                                  else
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                    }

                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"C->G"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_G<<"\t"<<base_Gfw<<"\t"<<base_Gr<<"\t"<<AF_Gfw<<"\t"<<AF_Gbw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'G')<<std::endl;
                                    //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"C->G"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_G<<"\t"<<base_Gfw<<"\t"<<base_Gr<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }

                              //T
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_T",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);
                                
                                //results --> 
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_T-base_Tr,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Tr,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_T-base_Tr,base_Tr);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Tfw<5 || base_Tr<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_T",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_T",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                        sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_T>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'T')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                    vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;

                                                  }
                                                  else
                                                  {
                                                    vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"C"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                    }

                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"C->T"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_T<<"\t"<<base_Tfw<<"\t"<<base_Tr<<"\t"<<AF_Tfw<<"\t"<<AF_Tbw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'T')<<std::endl;
                                    //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"C->T"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_T<<"\t"<<base_Tfw<<"\t"<<base_Tr<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }

                        }//end of ref C
                        else if(strcmp(reference,case_G)==0)
                        {
                            //reference is G

                             //A
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_A",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);
                                
                                //results -->
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_A-base_Ar,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Ar,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_A-base_Ar,base_Ar);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Afw<5 || base_Ar<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_A",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_A",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                        sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_A>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'A')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                    vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;

                                                  }
                                                  else
                                                  {
                                                    vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                    }

                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"G->A"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_A<<"\t"<<base_Afw<<"\t"<<base_Ar<<"\t"<<AF_Afw<<"\t"<<AF_Abw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'A')<<std::endl;
                                   //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"G->A"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_A<<"\t"<<base_Afw<<"\t"<<base_Ar<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }

                            //C
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_C",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);
                                
                                //results -->
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_C-base_Cr,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Cr,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_C-base_Cr,base_Cr);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Cfw<5 || base_Cr<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_C",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_C",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                        sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_C>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'C')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                    vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                                  }
                                                  else
                                                  {
                                                    vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                    }

                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"G->C"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_C<<"\t"<<base_Cfw<<"\t"<<base_Cr<<"\t"<<AF_Cfw<<"\t"<<AF_Cbw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'C')<<std::endl;
                                    //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"G->C"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_C<<"\t"<<base_Cfw<<"\t"<<base_Cr<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }

                                //T
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_T",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);
                                
                                //results -->
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_T-base_Tr,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Tr,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_T-base_Tr,base_Tr);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Tfw<5 || base_Tr<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_T",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_T",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                       sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_T>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'T')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                    vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;

                                                  }
                                                  else
                                                  {
                                                    vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"G"<<"\t"<<"T"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_T<<";"<<RD<<";"<<base_Tfw+base_Tr<<std::endl;
                                    }


                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"G->T"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_T<<"\t"<<base_Tfw<<"\t"<<base_Tr<<"\t"<<AF_Tfw<<"\t"<<AF_Tbw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'T')<<std::endl;
                                   //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"G->T"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_T<<"\t"<<base_Tfw<<"\t"<<base_Tr<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }

                            
                        }//end of ref G
                        else if(strcmp(reference,case_T)==0)
                        {
                          //reference is T
                            
                               //A
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_A",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);
                                
                                //results ->
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_A-base_Ar,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Ar,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_A-base_Ar,base_Ar);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Afw<5 || base_Ar<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_A",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_A",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                        sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_A>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'A')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;

                                                  }
                                                  else
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"A"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_A<<";"<<RD<<";"<<base_Afw+base_Ar<<std::endl;
                                    }



                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"T->A"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_A<<"\t"<<base_Afw<<"\t"<<base_Ar<<"\t"<<AF_Afw<<"\t"<<AF_Abw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'A')<<std::endl;
                                    //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"T->A"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_A<<"\t"<<base_Afw<<"\t"<<base_Ar<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }

                               //C
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_C",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);

                                //results -->
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_C-base_Cr,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Cr,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_C-base_Cr,base_Cr);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Cfw<5 || base_Cr<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }
                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_C",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_C",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                        sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_C>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'C')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;

                                                  }
                                                  else
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"C"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_C<<";"<<RD<<";"<<base_Cfw+base_Cr<<std::endl;
                                    }

                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"T->C"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_C<<"\t"<<base_Cfw<<"\t"<<base_Cr<<"\t"<<AF_Cfw<<"\t"<<AF_Cbw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'C')<<std::endl;
                                    //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"T->C"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_C<<"\t"<<base_Cfw<<"\t"<<base_Cr<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }

                               //G
                          memset(prompt_key,0,50);
                          sprintf(prompt_key,"%s_G",ref_key);
                          got_Thresholds=Thresholds.find(prompt_key);
                          if(got_Thresholds==Thresholds.end())
                          {
                            std::cout<<"malakia"<<std::endl;
                          }
                          else
                          {
                                memset(threshold_values,0,50);
                                memset(AF_fw_char,0,50);
                                memset(AF_bw_char,0,50);
                              
                                sprintf(threshold_values,"%s",got_Thresholds->second.c_str());
                                sscanf(threshold_values,"%[^_]_%[^_]",AF_fw_char,AF_bw_char);
                                AF_fw=std::stof(AF_fw_char);
                                AF_bw=std::stof(AF_bw_char);
                                
                                //results -->
                                Q_fw=0;
                                Q_bw=0;
                                Q_fw=mutationRulesPoissonQualityScore(base_G-base_Gr,RD-RD_reverse,AF_fw);
                                Q_bw=mutationRulesPoissonQualityScore(base_Gr,RD_reverse,AF_bw);
                                if(FW>=100 && BW>=100 && Q_fw>=5 && Q_bw>=5)
                                {
                                    //you have a variant so do the fisher test and write it
                                    p=-1;
                                    p=fisherTest(RD-RD_reverse,RD_reverse,base_G-base_Gr,base_Gr);
                                    if(p<=0.05)
                                    {
                                        sprintf(Flag_Fisher,"YES");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Fisher,"NO");
                                    }

                                    if(base_Gfw<5 || base_Gr<5)
                                    {
                                        sprintf(Flag_Tier,"LowQual");
                                    }
                                    else
                                    {
                                        sprintf(Flag_Tier,"HighQual");
                                    }

                                    char GermlineFlag[50];
                                    memset(GermlineFlag,0,50);
                                    char germline_prompt_key[50];
                                    memset(germline_prompt_key,0,50);
                                    sprintf(germline_prompt_key,"%s_%s_G",chrom,position);
                                    //query the data structures
                                    got_Germline_Hash_forPatients=Germline_Hash_forPatients.find(germline_prompt_key);
                                    if(got_Germline_Hash_forPatients==Germline_Hash_forPatients.end())
                                    {
                                            sprintf(GermlineFlag,"-");
                                    }
                                    else
                                    {   
                                            sprintf(GermlineFlag,"%s",got_Germline_Hash_forPatients->second.c_str());
                                    }

                                    
                                    char MaxGermlineFlag[50];
                                    memset(MaxGermlineFlag,0,50);
                                    char max_germ_prompt_key[50];
                                    memset(max_germ_prompt_key,0,50);
                                    sprintf(max_germ_prompt_key,"%s_%s_G",chrom,position);
                                    got_Germline_Max_Hash=Germline_Max_Hash.find(max_germ_prompt_key);
                                    if(got_Germline_Max_Hash==Germline_Max_Hash.end())
                                    {
                                        sprintf(MaxGermlineFlag,"-");
                                    }
                                    else
                                    {   
                                    
                                        sprintf(MaxGermlineFlag,"%s",got_Germline_Max_Hash->second.c_str());
                            
                                        
                                    }

                                    //here we need to add the upstream and downstream 10-mers
                                    int position_prompt_key=0;
                                    position_prompt_key=std::atoi(position);
                                    //store the 10-mers
                                    char down[50];
                                    memset(down,0,50);
                                    char up[50];
                                    memset(up,0,50);
                                    find_kmer_down(ReferenceBase_Hash,chrom,position_prompt_key,down);
                                    find_kmer_up(ReferenceBase_Hash,chrom,position_prompt_key,up);
                                    double Q=double(Q_fw+Q_bw)/2.000;

                                    memset(cat_dup_fisher,0,50);
                                    sprintf(cat_dup_fisher,"%s_%s",Flag_Dup,Flag_Fisher);
                                    max_germ=std::atof(MaxGermlineFlag);
                                    if(strcmp(cat_dup_fisher,"NO_NO")==0)
                                    {
                                        if(strcmp(Flag_Tier,"HighQual")==0)
                                        {
                                            if(AF_G>=max_germ)
                                            {
                                                if(homopolymerTest(down,up,'G')==1)
                                                {
                                                  vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"HomoPolymerRegion"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                                }
                                                else
                                                {

                                                  if(Q_fw>=20 && Q_bw>=20)
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"PASS"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;

                                                  }
                                                  else
                                                  {
                                                      vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"LowQscore"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;

                                                  }

                                                }
                                            }
                                            else
                                            {
                                                 vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"PositionWithHighNoise"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                            }
                                        }
                                        else
                                        {
                                            vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"LowSupportingReads"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                        }
                                    }
                                    else if(strcmp(cat_dup_fisher,"YES_NO")==0)
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"AmpliconEdge"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                    }
                                    else
                                    {
                                        vcf_output<<chrom<<"\t"<<position<<"\t"<<"-"<<"\t"<<"T"<<"\t"<<"G"<<"\t"<<Q<<"\t"<<"StrandBias"<<"\t"<<AF_G<<";"<<RD<<";"<<base_Gfw+base_Gr<<std::endl;
                                    }

                                    //write the variant
                                    output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"T->G"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_G<<"\t"<<base_Gfw<<"\t"<<base_Gr<<"\t"<<AF_Gfw<<"\t"<<AF_Gbw<<"\t"<<Flag_Dup<<"_"<<Flag_Fisher<<"\t"<<p<<"\t"<<std::setprecision(4)<<Q_fw<<"\t"<<std::setprecision(4)<<Q_bw<<"\t"<<Flag_Tier<<"\t"<<GermlineFlag<<"\t"<<MaxGermlineFlag<<"\t"<<down<<"\t"<<up<<"\t"<<homopolymerTest(down,up,'G')<<std::endl;
                                    //output<<it->second<<"\t"<<chrom<<"\t"<<position<<"\t"<<"T->G"<<"\t"<<RD<<"\t"<<FW<<"\t"<<BW<<"\t"<<AF_G<<"\t"<<base_Gfw<<"\t"<<base_Gr<<"\t"<<std::setprecision(10)<<Q_fw<<"\t"<<std::setprecision(10)<<Q_bw<<std::endl;

                                }
                            }

                        }//end of ref T
                        else
                        {
                          std::cout<<"malakia"<<std::endl;
                        }   
                    }  

                    }                 
            }//end of if after reading a line
          }
        }
        vcf_output.close();
        output.close();
    }
    
}

void find_kmer_down(std::unordered_map<std::string,std::string> &Hash, char *chromosome, int position,char *output)
{

    char prompt_key[50];

    int down_pos_10=position-10;
    int down_pos_9=position-9;
    int down_pos_8=position-8;
    int down_pos_7=position-7;
    int down_pos_6=position-6;
    int down_pos_5=position-5;
    int down_pos_4=position-4;
    int down_pos_3=position-3;
    int down_pos_2=position-2;
    int down_pos_1=position-1;
    
    std::string down;
    
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_10;
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_9;
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_8;
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_7;
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_6;
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_5;
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_4;
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_3;
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_2;
    std::unordered_map<std::string,std::string>::iterator got_Base_down_pos_1;
    
    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_10);
    got_Base_down_pos_10=Hash.find(prompt_key);
    if(got_Base_down_pos_10==Hash.end())
    {
        down="-|";
    }
    else
    {
        down=got_Base_down_pos_10->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_9);
    got_Base_down_pos_9=Hash.find(prompt_key);
    if(got_Base_down_pos_9==Hash.end())
    {
        down=down+"-|";
    }
    else
    {
        down=down+got_Base_down_pos_9->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_8);
    got_Base_down_pos_8=Hash.find(prompt_key);
    if(got_Base_down_pos_8==Hash.end())
    {
        down=down+"-|";
    }
    else
    {
        down=down+got_Base_down_pos_8->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_7);
    got_Base_down_pos_7=Hash.find(prompt_key);
    if(got_Base_down_pos_7==Hash.end())
    {
        down=down+"-|";
    }
    else
    {
        down=down+got_Base_down_pos_7->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_6);
    got_Base_down_pos_6=Hash.find(prompt_key);
    if(got_Base_down_pos_6==Hash.end())
    {
        down=down+"-";
    }
    else
    {
        down=down+got_Base_down_pos_6->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_5);
    got_Base_down_pos_5=Hash.find(prompt_key);
    if(got_Base_down_pos_5==Hash.end())
    {
        down=down+"-|";
    }
    else
    {
        down=down+got_Base_down_pos_5->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_4);
    got_Base_down_pos_4=Hash.find(prompt_key);
    if(got_Base_down_pos_4==Hash.end())
    {
        down=down+"-|";
    }
    else
    {
        down=down+got_Base_down_pos_4->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_3);
    got_Base_down_pos_3=Hash.find(prompt_key);
    if(got_Base_down_pos_3==Hash.end())
    {
        down=down+"-";
    }
    else
    {
        down=down+got_Base_down_pos_3->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_2);
    got_Base_down_pos_2=Hash.find(prompt_key);
    if(got_Base_down_pos_2==Hash.end())
    {
        down=down+"-|";
    }
    else
    {
        down=down+got_Base_down_pos_2->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,down_pos_1);
    got_Base_down_pos_1=Hash.find(prompt_key);
    if(got_Base_down_pos_1==Hash.end())
    {
        down=down+"-";
    }
    else
    {
        down=down+got_Base_down_pos_1->second;
    }

    sprintf(output,"%s",down.c_str());

}


void find_kmer_up(std::unordered_map<std::string,std::string> &Hash, char *chromosome, int position,char *output)
{

    char prompt_key[50];

    int up_pos_1=position+1;
    int up_pos_2=position+2;
    int up_pos_3=position+3;
    int up_pos_4=position+4;
    int up_pos_5=position+5;
    int up_pos_6=position+6;
    int up_pos_7=position+7;
    int up_pos_8=position+8;
    int up_pos_9=position+9;
    int up_pos_10=position+10;

    std::string up; 

    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_1;
    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_2;
    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_3;
    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_4;
    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_5;
    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_6;
    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_7;
    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_8;
    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_9;
    std::unordered_map<std::string,std::string>::iterator got_Base_up_pos_10;

    
    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_1);
    got_Base_up_pos_1=Hash.find(prompt_key);
    if(got_Base_up_pos_1==Hash.end())
    {
        up="-|";
    }
    else
    {
        up=got_Base_up_pos_1->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_2);
    got_Base_up_pos_2=Hash.find(prompt_key);
    if(got_Base_up_pos_2==Hash.end())
    {
        up=up+"-|";
    }
    else
    {
        up=up+got_Base_up_pos_2->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_3);
    got_Base_up_pos_3=Hash.find(prompt_key);
    if(got_Base_up_pos_3==Hash.end())
    {
        up=up+"-|";
    }
    else
    {
        up=up+got_Base_up_pos_3->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_4);
    got_Base_up_pos_4=Hash.find(prompt_key);
    if(got_Base_up_pos_4==Hash.end())
    {
        up=up+"-|";
    }
    else
    {
        up=up+got_Base_up_pos_4->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_5);
    got_Base_up_pos_5=Hash.find(prompt_key);
    if(got_Base_up_pos_5==Hash.end())
    {
        up=up+"-|";
    }
    else
    {
        up=up+got_Base_up_pos_5->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_6);
    got_Base_up_pos_6=Hash.find(prompt_key);
    if(got_Base_up_pos_6==Hash.end())
    {
        up=up+"-|";
    }
    else
    {
        up=up+got_Base_up_pos_6->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_7);
    got_Base_up_pos_7=Hash.find(prompt_key);
    if(got_Base_up_pos_7==Hash.end())
    {
        up=up+"-|";
    }
    else
    {
        up=up+got_Base_up_pos_7->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_8);
    got_Base_up_pos_8=Hash.find(prompt_key);
    if(got_Base_up_pos_8==Hash.end())
    {
        up=up+"-|";
    }
    else
    {
        up=up+got_Base_up_pos_8->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_9);
    got_Base_up_pos_9=Hash.find(prompt_key);
    if(got_Base_up_pos_9==Hash.end())
    {
        up=up+"-|";
    }
    else
    {
        up=up+got_Base_up_pos_9->second;
    }

    memset(prompt_key,0,50);
    sprintf(prompt_key,"%s_%d",chromosome,up_pos_10);
    got_Base_up_pos_10=Hash.find(prompt_key);
    if(got_Base_up_pos_10==Hash.end())
    {
        up=up+"-";
    }
    else
    {
        up=up+got_Base_up_pos_10->second;
    }

    //std::cout<<"Found the 5-mer: "<<down<<std::endl;
    sprintf(output,"%s",up.c_str());
}

int homopolymerTest(char *down, char *up,char sub)
{
    int count_A=0;
    int count_C=0;
    int count_G=0;
    int count_T=0;

    if(sub=='A')
    {
        count_A=1;
    }
    else if(sub=='C')
    {
        count_C=1;
    }
    else if(sub=='G')
    {
        count_G=1;
    }
    else if(sub=='T')
    {
        count_T=1;
    }

    for(int idx=0;idx<=strlen(down);idx++)
    {
        if(down[idx]=='A')
        {
            count_A++;
        }
        else if (down[idx]=='C')
        {
            count_C++;
        }
        else if (down[idx]=='G')
        {
            count_G++;
        }
        else if (down[idx]=='T')
        {
            count_T++;
        }
    }
    for(int idx=0;idx<=strlen(up);idx++)
    {
        if(up[idx]=='A')
        {
            count_A++;
        }
        else if (up[idx]=='C')
        {
            count_C++;
        }
        else if (up[idx]=='G')
        {
            count_G++;
        }
        else if (up[idx]=='T')
        {
            count_T++;
        }
    }

    //std::cout<<"Returned up     :"<<up<<std::endl;
    //std::cout<<"Returned down   :"<<down<<std::endl;
    //std::cout<<"Found  A: "<<count_A<<"  C: "<<count_C<<"  G: "<<count_G<<"  T: "<<count_T<<std::endl;
    
    //the criterion is taken from Platypus and modified 
    //To avoid calling these variants, we compute a sequence complexity statistic that measures the
    //contribution of the two most frequent nucleotides among the 21 around a site; if this measure
   //exceeds 90%, SNP calls are flagged as suspicious.
   //Which is actually >18 

    if( (count_A+count_C) > 18)
    {
        return 1;
    }
    else if ((count_A+count_G) > 18)
    {
        return 1;
    }
    else if ((count_A+count_T) > 18)
    {
        return 1;
    }
    else if ((count_C+count_G) > 18)
    {
        return 1;
    }
    else if ((count_C+count_T) > 18)
    {
        return 1;
    }
    else if ((count_T+count_G) > 18)
    {
        return 1;
    }
    else
    {
        return 0;
    }
    

}
