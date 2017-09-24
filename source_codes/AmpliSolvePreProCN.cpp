/******************************************************************************************************
					AmpliSolve: detection of SNVs and CN aberations from Ion Torrent
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
     Email to dimitrios.kleftogiannis@icr.ac.uk 
     I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
     You are free to modify, extend or distribute this code, as long 
     as this copyright notice is included whole and unchanged. 

END COPYRIGHT NOTICE
 
UTILITY
  This is a pre-processing program required for estimating CN aberrations in a chromosome.
  The program takes as input some arguments (see below for details) and produces a tab-limited file tha contains
  the coverage ratio per amplicon for the population of all normal samples used as baseline. 
  This information can be further visualized as boxplot or violin plot using R, Matlab, Python etc. tools.
  Please note that this current version does not support data visualization.

INPUT ARGUMENTS
    1. panel design          : a tab limited bed-like file with info about the custom gene panel. Typically this file has 6 columns
                              as follows: 
                                    chrom \t amplicon_start \t amplicon_end \t info \t info \t gene_name
                              
                              if gene name is not available or there is an intronic region the gene_name field is N/A

    2. Germline dir          : a directory with germline bam files. If no germline files are available write: not_available 

    3. mapping quality Q     : mapping quality Q. Default is 20 for Ion Torrent

    4. chromosome            : the chromosome you want to test

    5. output dir            : directory for storing the intermediate results of the execution
    

DEPENDENCIES
    The AmpliSolve program depends on the following programs:
    
    1. SAM tools http://samtools.sourceforge.net/
    
    2. Bed tools http://bedtools.readthedocs.io/en/latest/

    Make sure that both are intalled and also included in the PATH of your system.
    If not you may need a command like: PATH=$PATH:/your/path/to/Samtools

    During the program development and testing Samtools v1.3.1 and bedtools v2.17.9 were used.
   
    3. ASEQ https://demichelislab.unitn.it/doku.php?id=public:aseq

    In fact the read count procedure has been wrapper around aseq software that utilizes
    internally samtools libraries. We acknowledge ASEQ implementation. 

    The provided makefile compiles all and actually produces
    binary computeCounts that computes the read counts per position in the panel for a bam file.
    Program computCounts can be executed alone and incorportated into other pipelines.
    
    4. Boost libraries.
    You have simply to download them from http://www.boost.org/ and give the path uppon compilation.
    In the boost web-site there are also more instructions and details.

    More details about AmpliSolve, descriptions and toy execution examples are provided in our web repository. 

    This is a C++ program, thus you need the C++ compiler to be installed in your computer.
    The AmpliSolve program has been developed in a Mac OS computer with El Capitan version 10.11.5
    The C++ compiler used for development is the clang C, C++ and Objective-C compiler
    The program works for Unix-like systems and has been tested not extensively.  

COMPILATION
    Type the following if you want to compile this code only:

    cc -o  AmpliSolveCN  AmpliSolveCN.cpp -std=c++0x -lstdc++

    If you want to install AmpliSolve tool and generate all binaries this is an example:

    make -f XXXXX BOOST_FLAG=YYYYY

    where XXXXX is the makefile for your operating system (linux or macosx)

    and YYYYY is the full directory of your BOOST libraries. For example /Users/kleftogiannis/boost_1_61_0

*/


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <time.h>
#include <sys/time.h>
#include <inttypes.h>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <map>
#include <cmath>

//color the output
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[36m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"

//definitions of buffers
#define MAXIMUM_READ_LENGTH 10000
#define MINIMUM_READ_LENGTH 1000

//functions
void printUsage();
void generateFolder(char *dir_path);
void generateBAMList(char *dir_path,char *list_name);
void myExec(char *command);
std::string myExec_v2(const char* cmd);
void storeBAMList(char *list_name, char *COUNT_DIR,std::unordered_map<std::string,std::string> &Hash);
int storePanel(char *amplicon_coordinates,std::map<int,std::string> &Hash, char *CHROM);
void computeNumberOfMappedReads(char *panel_design,std::unordered_map<std::string,std::string> &BAM_Hash,int quality,std::unordered_map<std::string,int> &RDS_Hash,char *result_dir);
void computeAmpliconCoverage(std::unordered_map<std::string,std::string> &BAM_Hash,int Q, std::unordered_map<std::string,int> &Cov_Hash, std::map<int,std::string> &Amplicon_Hash,char *result_dir);
void outputTable(std::unordered_map<std::string,std::string> &BAM_Hash,std::map<int,std::string> &Amplicon_Hash,std::unordered_map<std::string,int> &Cov_Hash,std::unordered_map<std::string,int> &RDS_Hash,char *result_dir,char *chr);
void simplePrintFunction(std::unordered_map<std::string,int> &Hash,char *result_dir);

//data structures
//store the germline BAM files with complete paths and filenames
std::unordered_map<std::string,std::string> GermlineBAMFileList_Hash;
std::unordered_map<std::string,std::string>::iterator got_GermlineBAMFileList_Hash;
//store the amplicons: the key is an integer number to keep the order of amplicons
std::map<int,std::string> PanelDesign_Hash;
std::map<int,std::string>::iterator got_PanelDesign_Hash;
//store the number of reads mapped on the amplicons
std::unordered_map<std::string,int> NumReads_Hash;
std::unordered_map<std::string,int>::iterator got_NumReads_Hash;
//store the number of reads mapped on specific amplicons of interest
std::unordered_map<std::string,int> AmpliconNumReads_Hash;
std::unordered_map<std::string,int>::iterator got_AmpliconNumReads_Hash;

//here your main program starts
int main(int argc, char **argv)
{
        //initialize the input arguments, read the inputs from the command line and check for correctness
        char *user_panel_design;
        user_panel_design=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        char *user_germline_dir;
        user_germline_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        char *user_Q;
        user_Q=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        char *user_chromosome;
        user_chromosome=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        char *user_output_dir;
        user_output_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        if(argc!=6)
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
            char *panel_design;
            panel_design=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            char *germline_dir;
            germline_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            char *Q;
            Q=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            char *chromosome;
            chromosome=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            char *output_dir;
            output_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            //finally save the arguments from the input list
            strcpy(user_panel_design,argv[1]);
            memset(panel_design,0,MINIMUM_READ_LENGTH);
            sscanf(user_panel_design,"panel_design=%s",panel_design);

            strcpy(user_germline_dir,argv[2]);
            memset(germline_dir,0,MINIMUM_READ_LENGTH);
            sscanf(user_germline_dir,"germline_dir=%s",germline_dir);

            strcpy(user_Q,argv[3]);
            memset(Q,0,MINIMUM_READ_LENGTH);
            sscanf(user_Q,"Q=%s",Q);

            strcpy(user_chromosome,argv[4]);
            memset(chromosome,0,MINIMUM_READ_LENGTH);
            sscanf(user_chromosome,"chromosome=%s",chromosome);

            strcpy(user_output_dir,argv[5]);
            memset(output_dir,0,MINIMUM_READ_LENGTH);
            sscanf(user_output_dir,"output=%s",output_dir);

            //flag the germline dir value
            char flag_germline_dir[]="not_available";
            int no_germline=-1;
            //save the Q in integer format
            int Q_int=std::atoi(Q);

            char command[MINIMUM_READ_LENGTH];
            std::cout<<"************************************************************************************************************************************\n"<<std::endl;
            std::cout<<"                        Pre-processing program that generates the coverage ratios for all normal samples for specified chromosome\n"<<std::endl;
            std::cout<<"                        Copyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n"<<std::endl;
            std::cout<<"Execution started under the following parameters:"<<std::endl;
            std::cout<<"\t1. Panel design                           : "<<ANSI_COLOR_GREEN<<panel_design<<ANSI_COLOR_RESET<<std::endl;
            if(strcmp(germline_dir,flag_germline_dir)==0)
            {
                no_germline=-1;
                std::cout<<"\t2. Germline bam dir                       : "<<ANSI_COLOR_RED<<"NO germline BAMs available"<<std::endl;
            }
            else
            {
                std::cout<<"\t2. Germline bam dir                       : "<<ANSI_COLOR_GREEN<<germline_dir<<ANSI_COLOR_RESET<<std::endl;
                no_germline=1;
            }
            if(Q_int<=0)
            {
                Q_int=20;
                std::cout<<"\t3. Q                                      : "<<ANSI_COLOR_RED<<"User gave: "<<Q<<ANSI_COLOR_RESET<<" The value has changed to 20 as default for Ion Torrent"<<std::endl;
            }
            else
            {
                std::cout<<"\t3. Q                                      : "<<ANSI_COLOR_GREEN<<Q_int<<ANSI_COLOR_RESET<<std::endl;
            }
            std::cout<<"\t4. Chromosome                             : "<<ANSI_COLOR_GREEN<<chromosome<<ANSI_COLOR_RESET<<std::endl;
            std::cout<<"\t5. output                                 : "<<ANSI_COLOR_GREEN<<output_dir<<ANSI_COLOR_RESET<<std::endl;
            
            //check if germline bams are available
            if(no_germline==1)
            {

                //generate the folder for storing the main results
                generateFolder(output_dir);
                char interm_results_dir[200];
                memset(interm_results_dir,0,200);
                //and the intermediate results
                sprintf(interm_results_dir,"%s/%s_PreProCN_interm_files",output_dir,chromosome);
                generateFolder(interm_results_dir);

                //generate a list with all germline bam files
                char germline_bam_list_name[MINIMUM_READ_LENGTH];
                memset(germline_bam_list_name,0,MINIMUM_READ_LENGTH);
                sprintf(germline_bam_list_name,"%s/germline_bam_list.txt",interm_results_dir);
                std::cout<<"Running function generateBAMList    \t\t";
                generateBAMList(germline_dir,germline_bam_list_name);
                std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;
                
                //read and store the list of germline files 
                std::cout<<"Running function storeBAMList       \t\t";
                storeBAMList(germline_bam_list_name, germline_dir,GermlineBAMFileList_Hash);
                std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;

                //parse the panel design and store the amplicons that correspond to the specified chromosome
                int F=storePanel(panel_design,PanelDesign_Hash,chromosome);
                //check if user's chromosome is in the panel. if not report a message and exit
                if(F<=0)
                {
                    std::cout<<ANSI_COLOR_RED<<"                                         THE PROGRAM WILL TERMINATE: chromosome provided by user is not in the panel    "<<ANSI_COLOR_RESET<<std::endl;
                    std::cout<<"************************************************************************************************************************************"<<std::endl;
                    exit(0);

                }
                else
                {
                    //continue from here the rest of execution
                    //store the number of reads in the BAMs
                    //small change in the order of printing commands
                    computeNumberOfMappedReads(panel_design,GermlineBAMFileList_Hash,Q_int,NumReads_Hash,interm_results_dir);
                    std::cout<<"Running function computeNumberOfMappedReads     \t";
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;
                    std::cout<<"Running function computeAmpliconCoverage     \t";
                    computeAmpliconCoverage(GermlineBAMFileList_Hash,Q_int,AmpliconNumReads_Hash,PanelDesign_Hash,interm_results_dir);
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;
                    
                    //all values we need are store, so we need to parse the two Hash tables populated earlier
                    //and generate the tab limited file with all ratios
                    //columns are the amplicons and rows are the filenames 
                    //the value are X/Y, where X is the reads of the amplicon and Y the overall number of reads
                    std::cout<<"Running function outputTable       \t\t";
                    outputTable(GermlineBAMFileList_Hash,PanelDesign_Hash,AmpliconNumReads_Hash,NumReads_Hash,output_dir,chromosome);
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;
                    std::cout<<"Execution was successful. Results can be found at: "<<ANSI_COLOR_YELLOW<<output_dir<<"/"<<chromosome<<"_germline_ratios.txt"<<ANSI_COLOR_RESET<<std::endl;
                }
            }
            else
            {
                //we dont have germline..so we terminate
                std::cout<<ANSI_COLOR_RED<<"                                            THE PROGRAM WILL TERMINATE: no normal file is available    "<<ANSI_COLOR_RESET<<std::endl;

            }
            std::cout<<"************************************************************************************************************************************\n"<<std::endl;
            return 0;
        }
        //clear the memory
        GermlineBAMFileList_Hash.clear();
        PanelDesign_Hash.clear();   
        NumReads_Hash.clear();
        AmpliconNumReads_Hash.clear();
}

//print program's usage
void printUsage()
{
  std::cout<<"Please type the following: "<<std::endl;
    std::cout<<"\n./AmpliSolvePreProCN"<<ANSI_COLOR_GREEN<<" panel_design="<<ANSI_COLOR_RESET<<"/your/panel/design/bed/file"<<ANSI_COLOR_GREEN<<" germline_dir="<<ANSI_COLOR_RESET<<"/dir/with/germline/bam/files"<<ANSI_COLOR_GREEN<<" Q="<<ANSI_COLOR_RESET<<"mapping_quality_score"<<ANSI_COLOR_GREEN<<" chromosome="<<ANSI_COLOR_RESET<<"chr_to_estimate_CN"<<ANSI_COLOR_GREEN<<" output="<<ANSI_COLOR_RESET<<"/dir/name/to/store/the/results/and/intermediate/results"<<std::endl;
    std::cout<<"\nRemember that: "<<std::endl;
    std::cout<<"This code internally uses SAMtools and Bedtools; make sure that are installed and configured (in your PATH)"<<std::endl;
    std::cout<<"\nIf no germline bam files are available just type germline_dir=not_available ; If germline files are available should be named bam and not BAM (case-sensitive) "<<std::endl;
    std::cout<<"\nExecution example:"<<std::endl;
    std::cout<<"The program takes 5 input arguments, so please fill them similar to the following command:"<<std::endl;
    std::cout<<"./AmpliSolvePreProCN panel_design=/Users/dkleftog/Desktop/AmpliSolve_Execution_Example/myBED/example_panel_design.bed germline_dir=/Users/dkleftog/Desktop/AmpliSolve_Execution_Example/Germline_samples Q=20 chromosome=chrX output=/Users/dkleftog/Desktop/AmpliSolveTesting"<<std::endl;
    std::cout<<"\n\tIt is important to give the arguments in this order. Otherwise the program will crach!"<<std::endl;
    std::cout<<"\nMore info about the input data and detailed execution examples can be found at our web-repository"<<std::endl;
    std::cout<<"************************************************************************************************************************************"<<std::endl;
}

//generate folders to save the main results and the intermediate results
void generateFolder(char *dir_path)
{
    char command[MINIMUM_READ_LENGTH];
    memset(command,0,MINIMUM_READ_LENGTH);
    sprintf(command, "mkdir -p %s",dir_path);
    myExec(command);
        
}

//the function takes as input a directory name (the actual path) and creates a list of BAM files
void generateBAMList(char *dir_path,char *list_name)
{
  char command[MINIMUM_READ_LENGTH];
  memset(command,0,MINIMUM_READ_LENGTH);
  sprintf(command, "ls %s/*.bam > %s",dir_path,list_name);
  myExec(command);
  
}

//function that executes a command internally as a pipe
//It allows the program to continue if and only if the command is successful
//otherwise it terminates the program
void myExec(char *command)
{
  FILE *in;
  if(!(in = popen(command,"r")))
  {
    std::cout<<"\t\tSomething went wrong with pipe opening during: "<<command<<ANSI_COLOR_RED<<" The program will terminate."<<ANSI_COLOR_RESET<<std::endl;
    exit(0);
  }
  int stat=pclose(in); 
  if(WEXITSTATUS(stat)==127)
  {
    std::cout<<ANSI_COLOR_RED<<"                                  Sorry but Amplisolve cannot continue...Please resolve the DEPENDENCIES"<<ANSI_COLOR_RESET<<std::endl;
    std::cout<<"************************************************************************************************************************************"<<std::endl;
    exit(0);
  }
  else if(WEXITSTATUS(stat)==0)
  {
    //let it continue....this is the status we want
  }
  else
  {
    //for all the rest
    std::cout<<"\t\t\nSomething went wrong with: "<<command<<std::endl;
    std::cout<<ANSI_COLOR_RED<<"                                        Sorry but Amplisolve cannot continue..."<<ANSI_COLOR_RESET<<std::endl;
    std::cout<<"************************************************************************************************************************************"<<std::endl;
    exit(0);
  }

}


//this is the safe way of executing external commands, if you want to catch command's return status
std::string myExec_v2(const char* cmd) 
{
    char buffer[MINIMUM_READ_LENGTH];
    std::string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    try 
    {
        while (!feof(pipe)) 
        {
            if (fgets(buffer, MINIMUM_READ_LENGTH, pipe) != NULL)
                result += buffer;
        }
    } 
    catch (...) 
    {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    std::cout<<result<<std::endl;
    return result;
}

//store the list of BAM files for processing
void storeBAMList(char *list_name, char *COUNT_DIR,std::unordered_map<std::string,std::string> &Hash)
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
    printf("\tError from storeBAMList: Cannot open %s\n",list_name);
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
        strncpy(file_name,lineCharArray+size_DIR_name+1,size_current-size_DIR_name-5);
        //the key is tha actual filename
        //the value filed is the name
        Hash.insert(std::make_pair<std::string,std::string>(lineCharArray,file_name));
        count++;
      }
    }
  }
  fclose(input);
  //std::cout<<" Running function storeBAMList: "<<ANSI_COLOR_GREEN<< list_name <<ANSI_COLOR_RESET<<" stored with success. It contains "<<ANSI_COLOR_GREEN<<Hash.size()<<ANSI_COLOR_RESET<<" germline samples"<<std::endl;
}


//store tha amplicons that correspond to the chromosome given by the user
int storePanel(char *amplicon_coordinates,std::map<int,std::string> &Hash, char *CHROM)
{
    char lineCharArray[MAXIMUM_READ_LENGTH];
    FILE *input;
    int ret=-1;
    int count=1;
    char complete_name[MAXIMUM_READ_LENGTH];
    sprintf(complete_name,"%s",amplicon_coordinates);
    if((input=fopen(complete_name,"r"))==NULL)
    {
        printf("Error from function storePanel: Cannot open %s\n",complete_name);
        exit(0);
    }
    else
    {
        //make sure that there is no header
        input=fopen(complete_name,"r");
        while(!feof(input))
        {
            //read line by line
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            //then get the positions
            char chrom[50];
            char pos_start[50];
            char pos_end[50];
            char ampli_name[50];
            char rs_number[50];
            char gene_name[50];
            if(ret!=EOF)
            {
                memset(chrom,0,50);
                memset(pos_start,0,50);
                memset(pos_end,0,50);
                memset(ampli_name,0,50);
                memset(rs_number,0,50);
                memset(gene_name,0,50);
                sscanf(lineCharArray,"%s\t%s\t%s\t%s\t%s\t%s",chrom,pos_start,pos_end,ampli_name,rs_number,gene_name);
                //compare user's chromosome with the actual chromosome name in the file
                if(strcmp(CHROM,chrom)==0)
                {
                    //store only the ones you are interested
                    char value[50];
                    memset(value,0,50);
                    sprintf(value,"%s:%s-%s",chrom,pos_start,pos_end);
                    Hash[count]=value;
                    count++;
                } 
            }
        }
    }
    fclose(input);
    //std::cout<<"\tPanel design has parsed with success. It contains: "<<ANSI_COLOR_GREEN<< Hash.size() <<" amplicons"<<ANSI_COLOR_RESET<<std::endl;
    int a=Hash.size();
    return a;
}

//this function runs internally samtools and bedtools and computes the total number of mapped reads on amplicons
void computeNumberOfMappedReads(char *panel_design,std::unordered_map<std::string,std::string> &BAM_Hash,int quality,std::unordered_map<std::string,int> &RDS_Hash,char *result_dir)
{
    //gerenate a random number unique for each execution
    time_t t;
    int seed;
    srand((unsigned) time(&t));
    seed=floor(rand()%100);
    int count=1;
    //save the command
    char command[MINIMUM_READ_LENGTH];
    //parse the files and generate a file with two columns: file name and number of reads mapped on all amplicons
    for(std::unordered_map<std::string,std::string>::iterator it=BAM_Hash.begin(); it!=BAM_Hash.end();++it)
    {
        
        if(count%2==0)
        {
            std::cout<<"\t\tProcessed so far: "<<count<<" files "<<BAM_Hash.size()<<std::endl;
        }
        count++;
        //computes the number of reads from the BAM
        memset(command,0,MINIMUM_READ_LENGTH);
        sprintf(command,"samtools view -b -q %d -F4 %s | intersectBed -abam stdin -b  %s -u -bed | wc -l > %s/%d_rds_mapped_on_ampl.txt",quality,it->first.c_str(),panel_design,result_dir,seed);
        myExec(command);

        //read this file and 
        char lineCharArray[MAXIMUM_READ_LENGTH];
        FILE *input;
        int ret=-1;
        int count=1;
        char complete_name[MAXIMUM_READ_LENGTH];
        sprintf(complete_name,"%s/%d_rds_mapped_on_ampl.txt",result_dir,seed);
        if((input=fopen(complete_name,"r"))==NULL)
        {
            printf("Error from function computeNumberOfMappedReads: Cannot open %s\n",complete_name);
            exit(0);
        }
        else
        {
            //make sure that there is no header
            input=fopen(complete_name,"r");
            while(!feof(input))
            {
                //read line by line
                memset(lineCharArray,0,MINIMUM_READ_LENGTH);
                ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
                int rds=0;
                if(ret!=EOF)
                {
                    sscanf(lineCharArray,"%d",&rds);
                    RDS_Hash.insert(std::pair<std::string,int>(it->second,rds));
                    count++;
                }
            }
        }
        fclose(input);   
    }
}

//this function computes the coverage per amplicon and stores the results in a Hash
void computeAmpliconCoverage(std::unordered_map<std::string,std::string> &BAM_Hash,int quality, std::unordered_map<std::string,int> &Cov_Hash, std::map<int,std::string> &Amplicon_Hash,char *result_dir)
{
    //gerenate a random number unique for each execution
    time_t t;
    int seed;
    srand((unsigned) time(&t));
    seed=floor(rand()%100);

    //save the command
    char command[MINIMUM_READ_LENGTH];
    //parse the files and compute for every amplicon the coverage per amplicon
    for(std::unordered_map<std::string,std::string>::iterator it=BAM_Hash.begin(); it!=BAM_Hash.end();++it)
    {

        for(std::map<int,std::string>::iterator iit=Amplicon_Hash.begin(); iit!=Amplicon_Hash.end();++iit)
        {
            memset(command,0,MINIMUM_READ_LENGTH);
            sprintf(command,"samtools view -q %d %s %s | wc -l > %s/%d_cov_on_selected_amplicons.txt", quality,it->first.c_str(),iit->second.c_str(),result_dir,seed);
            myExec(command);
            //read this file just produed and store the value to the Hash
            char lineCharArray[MAXIMUM_READ_LENGTH];
            FILE *input;
            int ret=-1;
            int count=1;
            char complete_name[MAXIMUM_READ_LENGTH];
            sprintf(complete_name,"%s/%d_cov_on_selected_amplicons.txt",result_dir,seed);
            if((input=fopen(complete_name,"r"))==NULL)
            {
                printf("Error from function computeNumberOfMappedReads: Cannot open %s\n",complete_name);
                exit(0);
            }
            else
            {
                //make sure that there is no header
                input=fopen(complete_name,"r");
                while(!feof(input))
                {
                    //read line by line
                    memset(lineCharArray,0,MINIMUM_READ_LENGTH);
                    ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
                    int rds=0;
                    if(ret!=EOF)
                    {
                        sscanf(lineCharArray,"%d",&rds);
                        char key[100];
                        memset(key,0,100);
                        sprintf(key,"%s_%d",it->second.c_str(),iit->first);
                        Cov_Hash.insert(std::pair<std::string,int>(key,rds));
                        count++;
                    }
                }
            }
            fclose(input);   
        }
    }

}

//the function generate the output of the program, the results are stored to the output_dir given by the user
void outputTable(std::unordered_map<std::string,std::string> &BAM_Hash,std::map<int,std::string> &Amplicon_Hash,std::unordered_map<std::string,int> &Cov_Hash,std::unordered_map<std::string,int> &RDS_Hash,char *result_dir,char *chr)
{
    //write the file
    std::ofstream output;
    char output_name[MINIMUM_READ_LENGTH];
    memset(output_name,0,MINIMUM_READ_LENGTH);
    sprintf(output_name,"%s/%s_germline_ratios.txt",result_dir,chr);
    output.open(output_name);
    //write the header of the file
    output<<"Filename";
    for(std::map<int,std::string>::iterator it=Amplicon_Hash.begin(); it!=Amplicon_Hash.end();++it)
    {
        output<<"\t"<<it->second;
    }
    output<<std::endl;
    //parse again the hash tables and generate the actual matrix of ratio
    for(std::unordered_map<std::string,std::string>::iterator a=BAM_Hash.begin(); a!=BAM_Hash.end();++a)
    {

        char filename[100];
        memset(filename,0,100);
        sprintf(filename,"%s",a->second.c_str());
        output<<filename;
        //parse the amplicons
        for(std::map<int,std::string>::iterator b=Amplicon_Hash.begin(); b!=Amplicon_Hash.end();++b)
        {
            char prompt_key[100];
            memset(prompt_key,0,100);
            sprintf(prompt_key,"%s_%d",filename,b->first);
            //retrieve the number of reads on the amplicon
            got_AmpliconNumReads_Hash=Cov_Hash.find(prompt_key);
            if(got_AmpliconNumReads_Hash==Cov_Hash.end())
            {
                //something is wrong so write NaN
                output<<"\tNaN";
            }
            else
            {
                //no obtain the number of reads mapped on all amplicons to generate the ratio
                got_NumReads_Hash=RDS_Hash.find(filename);
                if(got_NumReads_Hash==RDS_Hash.end())
                {
                    //something is wrong here so write NaN
                    output<<"\tNaN";
                }
                else
                {
                    int rds=got_NumReads_Hash->second;
                    if(rds==0)
                    {
                        //we cannot devide with 0 so write NaN
                        output<<"\tNaN";
                    }
                    else
                    {
                        //we hae everything so we simply write the ratio
                        output<<"\t"<<float(got_AmpliconNumReads_Hash->second)/float(rds);
                    }
                }
            }

        }
        output<<std::endl;
    }
    output.close();
}


//debuging function
void simplePrintFunction(std::unordered_map<std::string,int> &Hash,char *result_dir)
{
    std::ofstream output;
    char output_name[MINIMUM_READ_LENGTH];
    memset(output_name,0,MINIMUM_READ_LENGTH);
    sprintf(output_name,"%s/test.txt",result_dir);
    output.open(output_name);
    for(std::unordered_map<std::string,int>::iterator it=Hash.begin(); it!=Hash.end();++it)
    {
        output<<it->first<<"\t"<<it->second<<std::endl;
    }
    output.close();
}


