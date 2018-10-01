/************************************************************************************************************************************************************************************************************
				                                AmpliSolveErrorEstimation

                Position-specific, stand-specific and nucleotide-specific error estimation
            for targeted deep-sequencing data (tailored and tested for Ion Torrent AmpliSeq )

BEGIN COPYRIGHT NOTICE

      AmpliSolveVariantCalling code -- (c) 2018 Dimitrios Kleftogiannis -- ICR -- www.icr.ac.uk

      The MIT License; 
      You may not use this file except in compliance with the License. 

      Copyright 2018 Dimitrios Kleftogiannis

      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
      to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
      and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
      FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
      WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

      Published reports of research using this code (or a modified version) should cite the 
      article that describes AmpliSolve, see below for more information. 
     
      Comments and bug reports are welcome. 

      Email to dimitrios.kleftogiannis@icr.ac.uk 
      I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
      You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 

END COPYRIGHT NOTICE
 
UTILITY

  This is a pre-processing program required for running the actual variant calling mode program.

  The program takes as input some arguments (see below for details) and produces a tab-limited file that contains the 
  error estimations per nucleotide (6 thresholds per nucleotide)


INPUT ARGUMENTS
    1. panel design          : a tab limited file with info about the gene panel. Typically this file has 6 columns
                               as follows: 
                                    chrom TAB amplicon_start TAB amplicon_end TAB info TAB info TAB gene_name
                              
                              if gene name is not available or there is an intro region the gene_name field is N/A

    2. reference genome      : the reference genome in fasta format. Has to be indexed (fai file) 

    3. germline count dir    : a DIR (just the full path not ending with '/') that contains the read counts for all available 
                               germline (normal) samples. The read count files follow the ASEQ format, and they must produced using
                               the supplementary program computeCounts (provided). 

                               For default error estimation give "not_available" as an argument with no quotes =)

    4. C value               : provide the C value, which is the normalization factor described in the manuscript. For ctDNA samples
                               we reccomend C=0.002 or C=0.001, but for solid tumour samples any value above 0.01 gives sufficient FDR


    5. Output dir            : directory for storing the intermediate results of this execution


DEPENDENCIES

    The AmpliSolveErrorEstimation program depends on the following programs:
    
    1. SAM tools http://samtools.sourceforge.net/
    

    Make sure that is intalled and also included in the PATH of your system.
    
    If not you may need a command like: PATH=$PATH:/your/path/to/Samtools to add the programs to the working PATH


    2. ASEQ https://demichelislab.unitn.it/doku.php?id=public:aseq
    
    (Supplementary script is provided)

    The computeCounts program is based on the ASEQ software that utilizes
    internally samtools libraries to compute the number of read counts per allele. 

    We acknowledge ASEQ implementation and the modified version is distributed using the original licence. 

    Program computCounts can be executed alone, or incorportated into other pipelines for parallel faster execution (e.g., bash script in a cluster).
    
    3. Boost libraries.
    
    Download them from http://www.boost.org/ and give the path uppon compilation in the makefile (instructions provided).
    
    In the boost web-site there are also more instructions and implementation details.

    More details about AmpliSolve, descriptions and toy execution examples are provided in our web repository. 

    This is a C++ program, thus you need the C++ compiler to be installed in your computer.
    
    The AmpliSolveErrorEstimation program has been developed in a Mac OS computer with El Capitan version 10.11.5
    
    The C++ compiler used for development is the clang C, C++ and Objective-C compiler
    
    The program works for Unix-like systems (e.g., Fedora) and has been tested but not extensively.  

COMPILATION

    Type the following if you want to compile only this code:

    cc -o  AmpliSolveErrorEstimation  AmpliSolveErrorEstimation.cpp -std=c++0x -lstdc++

    If you want to install AmpliSolve tool and generate all programs included type the following:

    make -f XXXXX BOOST_FLAG=YYYYY

    where XXXXX is the appropriate makefile for your operating system (linux or macosx)

    and YYYYY is the full directory of your BOOST libraries. For example /Users/kleftogiannis/boost_1_61_0

    At the end you might want to add all binaries generated with the previous command to your PATH: the command looks like PATH=$PATH:/write/your/dir/to/AmpliSolve/source_codes

************************************************************************************************************************************************************************************************************
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
#include <iomanip>
#include <math.h>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <random>

//color the output
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[36m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"

//definitions of read buffers
#define MAXIMUM_READ_LENGTH 10000
#define MINIMUM_READ_LENGTH 1000

//functions: definition of more functions than needed is included
//this is because the code is under developemnt and might be helpful to add more material

//print help
void printUsage();
//generate a simple VCF file that is required to perform computeCount function
int generateDummyVCFFile(char *panelDesign,char *results_dir);
//generates the panel reference bases in a txt
int generateReferenceBases(char *panelDesign, char *referenceGenome,char *results_dir);
//generate list of BAMs: the code works with BAMS as well...but it is slow
void generateBAMList(char *dir_path,char *list_name);
//store list of BAMs: the code works with BAMS as well...but it is slow
void storeBAMList(char *list_name, char *COUNT_DIR,std::unordered_map<std::string,std::string> &Hash);
//store reference bases
void storeReference(char *reference_bases_name,std::unordered_map<std::string,std::string> &Hash);
//store duplicate positions
void storeDuplicates(char *amplicon_positions,std::unordered_map<std::string,std::string> &Hash);
//store germline statistics for threshold estimation
void storeGermlineStatistics(std::unordered_map<std::string,std::string> &FILE_Hash, std::unordered_multimap<std::string,std::string> &Value_Hash_Thres,std::unordered_map<std::string,double> &Germline_M_Hash);
//run the computeCounts program and generate read counts for all gernline files
//slow execution this is by-passed in the current version
void produceReadCounts(char *input_vcf, char *input_threads,char *input_mbq,char *input_mrq, char *input_mdc, char *output_dir,std::unordered_map<std::string,std::string> &Hash);
//generate the read count list: we note that the files are named witht the extension .ASEQ.PILEUP as the original ASEQ software
void generateCountList(char *dir_path,char *list_name);
//store the list of read count files
void storeCountList(char *list_name, char *COUNT_DIR,std::unordered_map<std::string,std::string> &Hash);
//estimate the noise levels and output the results in flat files
void estimateThresholds(float norm_factor,std::unordered_map<std::string,std::string> &Position_Hash,std::unordered_multimap<std::string,std::string> &Value_Hash,std::unordered_map<std::string,std::string> &Thresholds,std::unordered_map<std::string,std::string> &CountHash,std::unordered_map<std::string,std::string> &RatioHash);
void generateFinalOutput(float C_value_float,char *panelDesign,std::unordered_map<std::string,std::string> &Reference_Hash,std::unordered_map<std::string,std::string> &Duplicate_Hash,std::unordered_map<std::string,std::string> &Thresholds,std::unordered_map<std::string,double> &Germline_M_Hash,char *output_dir,std::unordered_map<std::string,std::string> &CountHash,std::unordered_map<std::string,std::string> &RatioHash);
void generateFinalOutput_default(float C_value_float,char *panelDesign,std::unordered_map<std::string,std::string> &Reference_Hash,std::unordered_map<std::string,std::string> &Duplicate_Hash,char *output_dir);
//safe way of running commands
void myExec(char *command);
//generate folder
void generateFolder(char *dir_path);
//function for experiments and only
void shuffleCountList(char *list_name,char *COUNT_DIR,int index,int sample_size,float Cvalue);
//debug function
void printGermline_Debug(std::unordered_multimap<std::string,std::string> my_map);

//data structures

//storing the germline BAM files with complete paths and filenames
std::unordered_map<std::string,std::string> GermlineBAMFileList_Hash;
std::unordered_map<std::string,std::string>::iterator got_GermlineBAMFileList_Hash;

//store the count files
std::unordered_map<std::string,std::string> GermlineCountFileList_Hash;
std::unordered_map<std::string,std::string>::iterator got_GermlineCountFileList_Hash;

//store reference bases
std::unordered_map<std::string,std::string> ReferenceBase_Hash;
std::unordered_map<std::string,std::string>::iterator got_ReferenceBase_Hash;

//store duplicate positions 
std::unordered_map<std::string,std::string> DuplicatePosition_Hash;
std::unordered_map<std::string,std::string>::iterator got_DuplicatePosition_Hash;

//storing the germline AFs and supporting reads for all positions, all nucleotides per strand for all germlines
//this info is used for the threshold estimation next
std::unordered_multimap<std::string,std::string> GermlineValues_Hash_forThresholds;
std::unordered_multimap<std::string,std::string>::iterator got_GermlineValues_Hash_forThresholds;

//store the MAX AF for all germline positions
std::unordered_map<std::string,double> Germline_Max_Hash;
std::unordered_map<std::string,double>::iterator got_Germline_Max_Hash;

//store the count
std::unordered_map<std::string,std::string> Count_Hash;
std::unordered_map<std::string,std::string>::iterator got_Count_Hash;

//store the ratio Ri/RD
std::unordered_map<std::string,std::string> Ratio_Hash;
std::unordered_map<std::string,std::string>::iterator got_Ratio_Hash;

//store the threshold
std::unordered_map<std::string,std::string> Thresholds_Hash_Analytic;
std::unordered_map<std::string,std::string>::iterator got_Thresholds_Hash_Analytic;

int main(int argc, char **argv)
{
        //initialize the input arguments, read the inputs from the command line and check for correctness
        char *user_panel_design;
        user_panel_design=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        
        char *user_reference_genome;
        user_reference_genome=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        
        char *user_germline_dir;
        user_germline_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        
        char *user_C_value;
        user_C_value=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH)); 

        char *user_output_dir;
        user_output_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        //this version contains two extra parameters for experimentation purposes
        if(argc!=6)
        {
        std::cout<<"************************************************************************************************************************************"<<std::endl;
        std::cout<<ANSI_COLOR_RED<<"                                        Your input arguments are not correct !"<<ANSI_COLOR_RESET<<std::endl;
        std::cout<<"                         Copyright 2018 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n"<<std::endl;
        printUsage();
        return 0;
        }
        else
        {   
            //store the input variable only if the number of input variables is correct
            //also this part makes some corrections for threads, mbq, mrq and mdc if not correct (i.e., when user gives <0)
            char *panel_design;
            panel_design=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            
            char *reference_genome;
            reference_genome=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            
            char *germline_dir;
            germline_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            
            char *C_value;
            C_value=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH)); 

            char *output_dir;
            output_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
    
            //assign the arguments to different variables
            strcpy(user_panel_design,argv[1]);
            memset(panel_design,0,MINIMUM_READ_LENGTH);
            sscanf(user_panel_design,"panel_design=%s",panel_design);

            strcpy(user_reference_genome,argv[2]);
            memset(reference_genome,0,MINIMUM_READ_LENGTH);
            sscanf(user_reference_genome,"reference_genome=%s",reference_genome);

            strcpy(user_germline_dir,argv[3]);
            memset(germline_dir,0,MINIMUM_READ_LENGTH);
            sscanf(user_germline_dir,"germline_dir=%s",germline_dir);

            strcpy(user_C_value,argv[4]);
            memset(C_value,0,MINIMUM_READ_LENGTH);
            sscanf(user_C_value,"C_value=%s",C_value);

            strcpy(user_output_dir,argv[5]);
            memset(output_dir,0,MINIMUM_READ_LENGTH);
            sscanf(user_output_dir,"output_dir=%s",output_dir);

            //convert the original values from string to numbers
            float C_value_float=std::atof(C_value);

            //flag the germline dir value
            char flag_germline_dir[]="not_available";
            //return 1 if germlines are available otherwise execute with the default errors
            int no_germlines=-1;


            std::cout<<"************************************************************************************************************************************\n"<<std::endl;
            std::cout<<"                                Error estimation required for AmpliSolveVariantCalling program \n"<<std::endl;
            std::cout<<"                        Copyright 2018 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n"<<std::endl;
            std::cout<<"Execution started under the following parameters:"<<std::endl;
            std::cout<<"\t1. Panel design                                   : "<<ANSI_COLOR_GREEN<<panel_design<<ANSI_COLOR_RESET<<std::endl;
            std::cout<<"\t2. Reference genome                               : "<<ANSI_COLOR_GREEN<<reference_genome<<ANSI_COLOR_RESET<<std::endl;
            
            if(strcmp(germline_dir,flag_germline_dir)==0)
            {
                no_germlines=-1;
                std::cout<<"\t3. Germline count dir                             : "<<ANSI_COLOR_RED<<"NO germline count files available"<<ANSI_COLOR_RESET<<". Estimation of error is based on default setting (Error=C_value)"<<std::endl;
            }
            else
            {
                std::cout<<"\t3. Germline count dir                             : "<<ANSI_COLOR_GREEN<<germline_dir<<ANSI_COLOR_RESET<<std::endl;
                no_germlines=1;
            }

            if(C_value_float<=0)
            {
              //this is the default C value in case the users give something non acceptable
                C_value_float=0.002;
                std::cout<<"\t4. C value                                         : "<<ANSI_COLOR_RED<<"User gave: "<<C_value_float<<ANSI_COLOR_RESET<<". The value is converted to 0.002"<<std::endl;
            }
            else 
            {
                 std::cout<<"\t4. C value                                        : "<<ANSI_COLOR_GREEN<<C_value_float<<ANSI_COLOR_RESET<<std::endl;
            }

            std::cout<<"\t5. Output dir                                     : "<<ANSI_COLOR_GREEN<<output_dir<<ANSI_COLOR_RESET<<std::endl;

            //at this point we have all input arguments and we can proceed further

            if(no_germlines==1)
            {

                //here I add an internal flag for my debuging procedure and because the work is in progress...
                int internal_flag=1;
                if(internal_flag==0)
                {
                    //this is a dummy if-statement
                    //the code is never gonna execute...used this if  during developemnt

                }
                else
                {
                    //store the intermediate results
                    char interm_results_dir[2000];
                    memset(interm_results_dir,0,2000);
                    sprintf(interm_results_dir,"%s/AmpliSolveErrorEstimation_interm_files",output_dir);
                    generateFolder(interm_results_dir);
                    std::cout<<"\nRunning function generateReferenceBases: ";
                    int seed_RefDup=0;
                    seed_RefDup=generateReferenceBases(panel_design,reference_genome,interm_results_dir);
                    
                    //this file name has generated before...
                    char reference_bases[MINIMUM_READ_LENGTH];
                    memset(reference_bases,0,MINIMUM_READ_LENGTH);
                    sprintf(reference_bases,"%s/%d_panelReferenceBases.txt",interm_results_dir,seed_RefDup);
                    //read and store the reference bases
                    //std::cout<<"Running function storeReference: ";
                    storeReference(reference_bases,ReferenceBase_Hash);

                    char duplicate_bases[MINIMUM_READ_LENGTH];
                    memset(duplicate_bases,0,MINIMUM_READ_LENGTH);
                    sprintf(duplicate_bases,"%s/%d_ampliconDuplicatedPositions.txt",interm_results_dir,seed_RefDup);
                    //read and store the duplicate positions
                    //std::cout<<"Running function storeDuplicates: ";
                    storeDuplicates(duplicate_bases,DuplicatePosition_Hash);
                    
                    char GERMLINE_DIR[MINIMUM_READ_LENGTH];
                    memset(GERMLINE_DIR,0,MINIMUM_READ_LENGTH);
                    sprintf(GERMLINE_DIR,"%s",germline_dir);
                    
                    //generate the list of ASEQ files and store the file names 
                    char germline_count_list_name[MINIMUM_READ_LENGTH];
                    memset(germline_count_list_name,0,MINIMUM_READ_LENGTH);
                    sprintf(germline_count_list_name,"%s/%d_germline_count_list_original.txt",interm_results_dir,seed_RefDup);
                    generateCountList(GERMLINE_DIR,germline_count_list_name);
                    
                    storeCountList(germline_count_list_name,GERMLINE_DIR,GermlineCountFileList_Hash);
                    
                    //store the germline statistics you need for estimating thresholds
                    std::cout<<"Running function storeGermlineStatistics:"<<std::endl;
                    storeGermlineStatistics(GermlineCountFileList_Hash,GermlineValues_Hash_forThresholds,Germline_Max_Hash);
                    //compute thresholds
                    std::cout<<"Running function estimateThresholds: ";
                    estimateThresholds(C_value_float,ReferenceBase_Hash,GermlineValues_Hash_forThresholds,Thresholds_Hash_Analytic,Count_Hash,Ratio_Hash);
                    //generate and write the output file
                    generateFinalOutput(C_value_float,panel_design,ReferenceBase_Hash,DuplicatePosition_Hash,Thresholds_Hash_Analytic,Germline_Max_Hash,output_dir,Count_Hash,Ratio_Hash);
                    
                    char output_name1[MINIMUM_READ_LENGTH];
                    memset(output_name1,0,MINIMUM_READ_LENGTH);
                    sprintf(output_name1,"%s/positionSpecific_%.4f.txt",output_dir,C_value_float);
                    std::cout<<"\nAmpliSolveErrorEstimation execution was successful. Results can be found at: "<<ANSI_COLOR_YELLOW<<output_name1<<ANSI_COLOR_RESET<<std::endl;

                    //clean the memory...
                    GermlineCountFileList_Hash.clear();
                    GermlineValues_Hash_forThresholds.clear();
                    Germline_Max_Hash.clear();
                    Thresholds_Hash_Analytic.clear();
                    Count_Hash.clear();
                    Ratio_Hash.clear();

                }

            }
            else
            {

                //store the results of the program and the intermediate files
                char interm_results_dir[2000];
                memset(interm_results_dir,0,2000);
                //and the intermediate results
                sprintf(interm_results_dir,"%s/AmpliSolveErrorEstimation_interm_files",output_dir);
                generateFolder(interm_results_dir);
                
                std::cout<<"Running function generateReferenceBases: "<<std::endl;
                int seed_RefDup=0;
                seed_RefDup=generateReferenceBases(panel_design,reference_genome,interm_results_dir);

                 //this file name has generated before...
                char reference_bases[MINIMUM_READ_LENGTH];
                memset(reference_bases,0,MINIMUM_READ_LENGTH);
                sprintf(reference_bases,"%s/%d_panelReferenceBases.txt",interm_results_dir,seed_RefDup);
                //read and store the reference bases
                std::cout<<"Running function storeReference: ";
                storeReference(reference_bases,ReferenceBase_Hash);
                
                char duplicate_bases[MINIMUM_READ_LENGTH];
                memset(duplicate_bases,0,MINIMUM_READ_LENGTH);
                sprintf(duplicate_bases,"%s/%d_ampliconDuplicatedPositions.txt",interm_results_dir,seed_RefDup);
                //read and store the duplicate positions
                std::cout<<"Running function storeDuplicates: ";
                storeDuplicates(duplicate_bases,DuplicatePosition_Hash);

                generateFinalOutput_default(C_value_float,panel_design,ReferenceBase_Hash,DuplicatePosition_Hash,output_dir);
                char output_name1[MINIMUM_READ_LENGTH];
                memset(output_name1,0,MINIMUM_READ_LENGTH);
                sprintf(output_name1,"%s/positionSpecificNoise_default.txt",output_dir);
                std::cout<<"\nAmpliSolveErrorEstimation execution was successful. Results can be found at: "<<ANSI_COLOR_YELLOW<<output_name1<<ANSI_COLOR_RESET<<std::endl;
            }

        }
        //clean the memory before you go home...
        GermlineBAMFileList_Hash.clear();
        GermlineCountFileList_Hash.clear();
        ReferenceBase_Hash.clear();
        DuplicatePosition_Hash.clear();
        GermlineValues_Hash_forThresholds.clear();
        Germline_Max_Hash.clear();
        Thresholds_Hash_Analytic.clear();
        Count_Hash.clear();
        Ratio_Hash.clear();
        std::cout<<"\n************************************************************************************************************************************"<<std::endl;
}

//print program's usage
void printUsage()
{
  std::cout<<"Please type the following: "<<std::endl;
    std::cout<<"\n./AmpliSolveErrorEstimation"<<ANSI_COLOR_GREEN<<" panel_design="<<ANSI_COLOR_RESET<<"/your/panel/design/file"<<ANSI_COLOR_GREEN<<" reference_genome="<<ANSI_COLOR_RESET<<"/your/reference/genome"<<ANSI_COLOR_GREEN<<" germline_dir="<<ANSI_COLOR_RESET<<"/your/dir/with/germline/read/count/files"<<ANSI_COLOR_GREEN<<" C_value="<<ANSI_COLOR_RESET<<"/the/error/rescaling/factor"<<ANSI_COLOR_GREEN<<" output_dir="<<ANSI_COLOR_RESET<<"/dir/to/store/the/intermediate/results"<<std::endl;
    std::cout<<"\nRemember that: "<<std::endl;
    std::cout<<"If NO germline count files are available, type germline_dir=not_available ; If germline files are available must follow the ASEQ data format and all files ending with .ASEQ"<<std::endl;
    std::cout<<"\nExecution example:"<<std::endl;
    std::cout<<"The program takes 5 input arguments, so please fill them similar to the following command:"<<std::endl;
    std::cout<<"./AmpliSolveErrorEstimation panel_design=/Example/panel_design.bed reference_genome=/Example/hg19_chr.fa germline_dir=/Example/Germline_ASEQ_DIR C_value=0.002 output_dir=/Example_Testing"<<std::endl;
    std::cout<<"\n\tIt is important to give the arguments in this order. Otherwise the program will crach !"<<std::endl;
    std::cout<<"\nMore info about the input data format and more execution examples can be found at our web-repository"<<std::endl;
    std::cout<<"************************************************************************************************************************************"<<std::endl;

}


//the function takes as input a directory name (the actual path)
//and creates a list of BAM files belong to the folder
void generateBAMList(char *dir_path,char *list_name)
{
  char command[MINIMUM_READ_LENGTH];
  memset(command,0,MINIMUM_READ_LENGTH);
  sprintf(command, "ls %s/*.bam > %s",dir_path,list_name);
  myExec(command);
    
}

//be careful here because the files must be ending as .ASEQ
void generateCountList(char *dir_path,char *list_name)
{
  char command[MINIMUM_READ_LENGTH];
  memset(command,0,MINIMUM_READ_LENGTH);
  sprintf(command, "ls %s/*.ASEQ > %s",dir_path,list_name);
  myExec(command);
    
}

//utility function that prints all thresholds; used mainly for debbuging...
void printThresholds(std::unordered_map<std::string,std::string> my_map)
{    
     char filename[50];
     memset(filename,0,50);
     sprintf(filename,"thresholds.txt");
     std::ofstream output;
     output.open(filename);
    for(std::unordered_map<std::string,std::string>::iterator it=my_map.begin(); it!=my_map.end();++it)
    {
      output<<it->first<<"\t" <<it->second<<std::endl;
    }
    
    output.close();   
}

//the implementation is very naive and takes ~5 min for a panel of ~40,000 positions
int generateReferenceBases(char *panelDesign, char *referenceGenome,char *results_dir)
{

    time_t t;
    int seed;
    srand((unsigned) time(&t));
    seed=floor(rand()%1000);

    //read buffer
    char lineCharArray[MINIMUM_READ_LENGTH];
    //file pointer to open the panel design and read it
    FILE *input;
    int ret=-1;
    //use this variable to generate and execute commands
    char command[500];

    int ampliconNumber=0;
    int positionCount=0;

    //write the position info
    std::ofstream output;
    char tmp_file[MINIMUM_READ_LENGTH];
    memset(tmp_file,MINIMUM_READ_LENGTH,0);
    sprintf(tmp_file,"%s/%d_tmp_info.txt",results_dir,seed);
    output.open(tmp_file);

    int count=0;
    //check if file exists before you open it
    if((input=fopen(panelDesign,"r"))==NULL)
    {
        printf("Error from generateReferenceBases: Cannot open file: %s\n",panelDesign);
        exit(0);
    }
    else
    {
        //open the file it is safe now
        input=fopen(panelDesign,"r");
        while(!feof(input))
        {
            //read line by line the amplicons from the panel design file
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            //get all the fields of the file
            char chrom[50];
            int position_start=-1;
            int position_end=-1;
            char ampliconID[50];
            char transcriptID[50];
            char geneName[50];
            //make sure that are clean
            memset(chrom,0,50);
            memset(ampliconID,0,50);
            memset(transcriptID,0,50);
            memset(geneName,0,50);
            //this is the actual function that reads and stores the fields into the variables
            sscanf(lineCharArray,"%s\t%d\t%d\t%s\t%s\t%s",chrom,&position_start,&position_end,ampliconID,transcriptID,geneName);
            ampliconNumber++;
                
                //we generate all positions in the amplicon and feed them to samtools faidx command
                for(int idx=position_start;idx<=position_end;idx++)
                {   

                    positionCount++;
                    count++;
                    //generate the command
                    memset(command,0,100);
                    sprintf(command,"samtools faidx %s %s:%d-%d >> %s/%d_tmp_ref.txt",referenceGenome,chrom,idx,idx,results_dir,seed);
                    //execute the command
                    myExec(command);
                    //write the position info
                    output<<chrom<<"\t"<<idx<<std::endl;
                }
        }
    }
    output.close();

    //file tmp_ref.txt contains all annotated positions but are in the fasta format we need to grep the actual bases and
    //then merge with the chrom and position info stored separately at tmp_info.txt
    memset(command,0,500);
    sprintf(command,"grep -v '>' %s/%d_tmp_ref.txt > %s/%d_tmp_ref_filtered.txt",results_dir,seed,results_dir,seed);
    myExec(command);
    memset(command,0,500);
    sprintf(command,"paste %s/%d_tmp_info.txt %s/%d_tmp_ref_filtered.txt > %s/%d_panelReferenceBases.txt",results_dir,seed,results_dir,seed,results_dir,seed);
    myExec(command);

    memset(command,0,500);
    sprintf(command,"sort %s/%d_panelReferenceBases.txt | uniq -d | cut -d$'\t' -f1,2 > %s/%d_ampliconDuplicatedPositions.txt",results_dir,seed,results_dir,seed);
    myExec(command);

    //write an output message at the end
    std::cout<<"Reference bases and amplicon duplicated positions have generated"<<"\n\t\t --> Parsed in total "<<ampliconNumber<<" amplicons and annotated "<<positionCount<<" positions."<<std::endl;
    return seed;
}

//the function reads the panel and generates all positions in a VCF-kind of format
int generateDummyVCFFile(char *panelDesign,char *results_dir)
{

    time_t t;
    int seed;
    srand((unsigned) time(&t));
    seed=floor(rand()%1000);

    char lineCharArray[MINIMUM_READ_LENGTH];
    FILE *input;
    int ret=-1;
    //check if file exists
    if((input=fopen(panelDesign,"r"))==NULL)
    {
        printf("Error from function generateDummyVCFFile: Cannot open %s\n",panelDesign);
        exit(0);
    }
    else
    {
        //open the file
        input=fopen(panelDesign,"r");
        //we generate the file to output the dymmy vcf
        //this file will be placed in the same folder as the BED_file
        std::ofstream output;
        char dummy_BED_complete_name[MINIMUM_READ_LENGTH];
        sprintf(dummy_BED_complete_name,"%s/%d_dummyVCF.vcf",results_dir,seed);
        output.open(dummy_BED_complete_name);
        
        while(!feof(input))
        {
            //read line by line
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            //then get the ranges and extract the positions
            char chrom[50];
            char start[50];
            char end[50];
            char ampli_name[50];
            char rs_number[50];
            char gene_name[50];
            if(ret!=EOF)
            {
                memset(chrom,0,50);
                memset(start,0,50);
                memset(end,0,50);
                memset(ampli_name,0,50);
                memset(rs_number,0,50);
                memset(gene_name,0,50);
                sscanf(lineCharArray,"%s\t%s\t%s\t%s\t%s\t%s",chrom,start,end,ampli_name,rs_number,gene_name);
                //debug
                //cout<<chrom<<":"<<start<<":"<<end<<endl;
                
                //at this point we have the ranges and we need to extract all single positions in the ranges
                //this will be done as follows
                int start_num=-1;
                int end_num=-1;
                int idx=-1;
                start_num=std::stoi(start);
                end_num=std::stoi(end);
                for(idx=start_num;idx<=end_num;idx++)
                {
                    //we fill the first positions we need and then we add . . . . Thats why we call it DUMMY =)
                    output<<chrom<<"\t"<<idx<<"\t.\t.\t.\t.\t.\t."<<std::endl;
                }
            }
        }
        //here you close the dummy file
        output.close();
    }

    std::cout<<" Dummy VCF file for the panel has generated -> "<<ANSI_COLOR_GREEN<<results_dir<<"/dummyVCF.vcf"<<ANSI_COLOR_RESET<<std::endl;

    return seed;
}

//store the list of files for processing
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
        //at this point be very careful...number 13 is based on the assumption that
        //the file names are PILEUP.ASEQ
        strncpy(file_name,lineCharArray+size_DIR_name+1,size_current-size_DIR_name-5);
        Hash.insert(std::make_pair<std::string,std::string>(lineCharArray,file_name));
        count++;
      }
    }
  }
  fclose(input);

  std::cout<<" Running function storeBAMList: "<<ANSI_COLOR_GREEN<< list_name <<ANSI_COLOR_RESET<<" stored with success. It contains "<<ANSI_COLOR_GREEN<<Hash.size()<<ANSI_COLOR_RESET<<" germline samples"<<std::endl;

}

//store the list of files for processing
void storeCountList(char *list_name, char *COUNT_DIR,std::unordered_map<std::string,std::string> &Hash)
{


    //store them first in a vector
    std::vector<std::string> V; 

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

}


//store the list of files for processing
void shuffleCountList(char *list_name,char *COUNT_DIR,int index,int sample_size,float Cvalue)
{


  //store them first in a vector
  std::vector<std::string> V; 

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
        V.push_back(lineCharArray);
        count++;
      }
    }
  }
  fclose(input);
  //perform randomization of the list
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  shuffle (V.begin(), V.end(), std::default_random_engine(seed));
  
  //TRAINING: write the results on a new list
  char outputFile[MINIMUM_READ_LENGTH];
  memset(outputFile,0,50);
  sprintf(outputFile,"training_%d_%d_%.7f.txt",index,sample_size,Cvalue);
  std::ofstream output;
  output.open(outputFile);

  //TESTING: write the results on a different list
  char testingFile[MINIMUM_READ_LENGTH];
  memset(testingFile,0,50);
  sprintf(testingFile,"testing_%d_%d_%.7f.txt",index,sample_size,Cvalue);
  std::ofstream output_2;
  output_2.open(testingFile);

  int myCounter=1; 
  for (std::vector<std::string>::const_iterator i = V.begin(); i != V.end(); ++i)
  {
    if(myCounter<=sample_size)
    {
        output << *i <<std::endl;
    }
    else
    {
        output_2 << *i <<std::endl;
    }
    myCounter++;
    
  }


  output.close();
  output_2.close();
    

  std::cout<<"\nRunning function shuffleCountList: "<<ANSI_COLOR_GREEN<< list_name <<ANSI_COLOR_RESET<<" stored shuffled and filtered with success. The new list contains "<<ANSI_COLOR_GREEN<<sample_size<<ANSI_COLOR_RESET<<" samples"<<std::endl;

}


//this function parses the panel reference bases (has to be pre-computed) and stores them in a Hash 
void storeReference(char *reference_bases_name,std::unordered_map<std::string,std::string> &Hash)
{
    char lineCharArray[MAXIMUM_READ_LENGTH];
    FILE *input;
    int ret=-1;
    int count=1;

    if((input=fopen(reference_bases_name,"r"))==NULL)
    {
        printf("Error from storeReference function: Cannot open %s\n",reference_bases_name);
        exit(0);
    }
    else
    {
        //read the file line by line
        input=fopen(reference_bases_name,"r");
        while(!feof(input))
        {
            //read line by line
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            //then extract the positions
            char chrom[50];
            char position[50];
            char base[10];
            if(ret!=EOF)
            {
                memset(chrom,0,50);
                memset(position,0,50);
                memset(base,0,10);
                sscanf(lineCharArray,"%s\t%s\t%s",chrom,position,base);
                
                //if you want to process only one chromosome specified by the user, change 1 with strcmp(chromosome,chrom)==0 and add chromosome as argument
                if(1)
                {
                    //at this point we have the position we need and we are storing in the hash
                    //the key is the concatenation of chrom_position and the value is just the position
                    char key[50];
                    memset(key,0,50);
                    sprintf(key,"%s_%s",chrom,position);
                    Hash.insert(std::make_pair<std::string,std::string>(key,base));
                    count++;
                }
            }
        }
    }
    fclose(input);
    std::cout<<"Running function storeReference: panel reference bases stored with success "<<ANSI_COLOR_GREEN<<Hash.size()<<ANSI_COLOR_RESET<<std::endl;
}

//this function parses the panel duplicate positions (have to be pre-computed) and stores them in a Hash 
void storeDuplicates(char *amplicon_positions,std::unordered_map<std::string,std::string> &Hash)
{

    char lineCharArray[MAXIMUM_READ_LENGTH];
    FILE *input;
    int ret=-1;
    int count=0;
    
    if((input=fopen(amplicon_positions,"r"))==NULL)
    {
        printf("Error from storeDuplicates function: Cannot open %s\n",amplicon_positions);
        exit(0);
    }
    else
    {
        //read the file line by line
        input=fopen(amplicon_positions,"r");
        while(!feof(input))
        {
            //read line by line
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            //then get the positions
            char chrom[50];
            char position[50];
            if(ret!=EOF)
            {
                memset(chrom,0,50);
                memset(position,0,50);
                sscanf(lineCharArray,"%s %s",chrom,position);
                
                //if you want to process only one chromosome specified by the user, change 1 with strcmp(chromosome,chrom)==0 and add chromosome as argument
                if(1)
                {
                    char key[50];
                    memset(key,0,50);
                    sprintf(key,"%s_%s",chrom,position);
                    Hash.insert(std::make_pair<std::string,std::string>(key,position));
                    count++;
                } 
            }
        }
    }
    fclose(input);
    std::cout<<"Running function storeDuplicates: panel duplicate positions stored with success "<<ANSI_COLOR_GREEN<<Hash.size()<<ANSI_COLOR_RESET<<std::endl;
}

//TODO this can be parallelized using OpenMP
void produceReadCounts(char *input_vcf, char *input_threads,char *input_mbq,char *input_mrq, char *input_mdc, char *output_dir,std::unordered_map<std::string,std::string> &Hash)
{

    char command[500];
    memset(command,0,500);
    int count=0;
    for(std::unordered_map<std::string,std::string>::iterator it=Hash.begin(); it!=Hash.end();++it)
    {
            memset(command,0,500);
            //binary computeCounts is the same as aseq; use exactly the same command
            sprintf(command,"computeCounts vcf=%s bam=%s threads=%s mbq=%s mrq=%s mdc=%s out=%s",input_vcf,it->first.c_str(),input_threads,input_mbq,input_mrq,input_mdc,output_dir);
            //std::cout<<command<<std::endl;
            count++;
            myExec(command);
            if(count % 5==0)
            {
                std::cout<<"\t\tWe have generated read counts for "<<ANSI_COLOR_GREEN<<count<<"/"<<Hash.size()<<ANSI_COLOR_RESET<<" germline samples"<<std::endl;
            }
    }
}



//parse all germline files and store VAFs from forward and reverse strand.
//this information is used for the threshold computation
void storeGermlineStatistics(std::unordered_map<std::string,std::string> &FILE_Hash, std::unordered_multimap<std::string,std::string> &Value_Hash_Thres,std::unordered_map<std::string,double> &Germline_M_Hash)
{    

    //this part of the code is used to avoid un-expected problems with chromosome naming
    //the code does not work if you name chrX and chrY as 23 and 24...
    //it accepts either X or chrX or 2, chr2 as long as this is consistent to all of your input files

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

     //simple iterator
     std::unordered_map<std::string,double>::iterator got_Germline_M_Hash;

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

                  //NOTE THAT THIS IS THE ASEQ PILEUP FORMAT
                  //THE SAMPLES MUST FOLLOW THIS FORMAT, OTHERWISE THE CODE DOES NOT WORK....
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
                    
                    //store all of the chrom, but if you want to focus on one just change 1 with strcmp(chromosome,chrom)==
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
                    
                            // key: chr_position_nucleotide --> e.g., chrX_12143125_A
                            memset(key,0,MINIMUM_READ_LENGTH);
                            sprintf(key,"%s_%s_A",chrom,position);
                            //value: Afw_FW_Abw_BW for the analytical
                            memset(value,0,MINIMUM_READ_LENGTH);
                            sprintf(value,"%d_%d_%d_%d",base_Afw,FW,base_Ar,BW);
                            Value_Hash_Thres.insert(std::make_pair(key,value));
                            
                            //now store the extra germline statistics we need....
                            if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0 )
                            {
                              //this setting is used for all the experiments in Ion Torrent; might change in another platform
                                    if(AF_A<=0.05 && FW>=100 && BW>=100)
                                    {
                                        char key_for_max[50];
                                        memset(key_for_max,0,50);
                                        sprintf(key_for_max,"%s_%s_A",chrom,position);

                                        got_Germline_M_Hash=Germline_M_Hash.find(key_for_max);
                                        if(got_Germline_M_Hash==Germline_M_Hash.end())
                                        {
                                            double value=-888;
                                            Germline_M_Hash[key_for_max]=value;
                                        }
                                        else
                                        {
                                            double value=got_Germline_M_Hash->second;
                                            if(value<=AF_A)
                                            {
                                                Germline_M_Hash[key_for_max]=AF_A;
                                            }
                                        }
                                    }
                            }
                            else
                            {
                                    if(AF_A<=0.05 && FW>=100 && BW>=100)
                                    {
                                        char key_for_max[50];
                                        memset(key_for_max,0,50);
                                        sprintf(key_for_max,"%s_%s_A",chrom,position);

                                        got_Germline_M_Hash=Germline_M_Hash.find(key_for_max);
                                        if(got_Germline_M_Hash==Germline_M_Hash.end())
                                        {
                                            double value=-888;
                                            Germline_M_Hash[key_for_max]=value;
                                        }
                                        else
                                        {
                                            double value=got_Germline_M_Hash->second;
                                            if(value<=AF_A)
                                            {
                                                Germline_M_Hash[key_for_max]=AF_A;
                                            }
                                        }
                                    }
                            }



                            memset(key,0,MINIMUM_READ_LENGTH);
                            sprintf(key,"%s_%s_C",chrom,position);
                            memset(value,0,MINIMUM_READ_LENGTH);
                            sprintf(value,"%d_%d_%d_%d",base_Cfw,FW,base_Cr,BW);
                            Value_Hash_Thres.insert(std::make_pair(key,value));

                            //now store the extra germline statistics we need....
                            if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
                            {
                                    if(AF_C<=0.05 && FW>=100 && BW>=100)
                                    {
                                        char key_for_max[50];
                                        memset(key_for_max,0,50);
                                        sprintf(key_for_max,"%s_%s_C",chrom,position);

                                        got_Germline_M_Hash=Germline_M_Hash.find(key_for_max);
                                        if(got_Germline_M_Hash==Germline_M_Hash.end())
                                        {
                                            double value=0;
                                            Germline_M_Hash[key_for_max]=value;
                                        }
                                        else
                                        {
                                            double value=got_Germline_M_Hash->second;
                                            if(value<=AF_C)
                                            {
                                                Germline_M_Hash[key_for_max]=AF_C;
                                            }
                                        }
                                    }
                            }
                            else
                            {
                                    if(AF_C<=0.05 && FW>=100 && BW>=100)
                                    {
                                        char key_for_max[50];
                                        memset(key_for_max,0,50);
                                        sprintf(key_for_max,"%s_%s_C",chrom,position);

                                        got_Germline_M_Hash=Germline_M_Hash.find(key_for_max);
                                        if(got_Germline_M_Hash==Germline_M_Hash.end())
                                        {
                                            double value=0;
                                            Germline_M_Hash[key_for_max]=value;
                                        }
                                        else
                                        {
                                            double value=got_Germline_M_Hash->second;
                                            if(value<=AF_C)
                                            {
                                                Germline_M_Hash[key_for_max]=AF_C;
                                            }
                                        }
                                    }
                            }

                            memset(key,0,MINIMUM_READ_LENGTH);
                            sprintf(key,"%s_%s_G",chrom,position);
                            memset(value,0,MINIMUM_READ_LENGTH);
                            sprintf(value,"%d_%d_%d_%d",base_Gfw,FW,base_Gr,BW);
                            Value_Hash_Thres.insert(std::make_pair(key,value));

                            //now store the extra germline statistics we need....
                            if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
                            {
                                    if(AF_G<=0.05 && FW>=100 && BW>=100)
                                    {
                                        char key_for_max[50];
                                        memset(key_for_max,0,50);
                                        sprintf(key_for_max,"%s_%s_G",chrom,position);

                                        got_Germline_M_Hash=Germline_M_Hash.find(key_for_max);
                                        if(got_Germline_M_Hash==Germline_M_Hash.end())
                                        {
                                            double value=0;
                                            Germline_M_Hash[key_for_max]=value;
                                        }
                                        else
                                        {
                                            double value=got_Germline_M_Hash->second;
                                            if(value<=AF_G)
                                            {
                                                Germline_M_Hash[key_for_max]=AF_G;
                                            }
                                        }
                                    }
                            }
                            else
                            {
                                    if(AF_G<=0.05 && FW>=100 && BW>=100)
                                    {
                                        char key_for_max[50];
                                        memset(key_for_max,0,50);
                                        sprintf(key_for_max,"%s_%s_G",chrom,position);

                                        got_Germline_M_Hash=Germline_M_Hash.find(key_for_max);
                                        if(got_Germline_M_Hash==Germline_M_Hash.end())
                                        {
                                            double value=0;
                                            Germline_M_Hash[key_for_max]=value;
                                        }
                                        else
                                        {
                                            double value=got_Germline_M_Hash->second;
                                            if(value<=AF_G)
                                            {
                                                Germline_M_Hash[key_for_max]=AF_G;
                                            }
                                        }
                                    }
                            }


                            memset(key,0,MINIMUM_READ_LENGTH);
                            sprintf(key,"%s_%s_T",chrom,position);
                            memset(value,0,MINIMUM_READ_LENGTH);
                            sprintf(value,"%d_%d_%d_%d",base_Tfw,FW,base_Tr,BW);
                            Value_Hash_Thres.insert(std::make_pair(key,value)); 

                            //now store the extra germline statistics we need....
                            if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
                            {
                                    if(AF_T<=0.05 && FW>=100 && BW>=100)
                                    {
                                        char key_for_max[50];
                                        memset(key_for_max,0,50);
                                        sprintf(key_for_max,"%s_%s_T",chrom,position);

                                        got_Germline_M_Hash=Germline_M_Hash.find(key_for_max);
                                        if(got_Germline_M_Hash==Germline_M_Hash.end())
                                        {
                                            double value=0;
                                            Germline_M_Hash[key_for_max]=value;
                                        }
                                        else
                                        {
                                            double value=got_Germline_M_Hash->second;
                                            if(value<=AF_T)
                                            {
                                                Germline_M_Hash[key_for_max]=AF_T;
                                            }
                                        }
                                    }
                            }
                            else
                            {
                                    if(AF_T<=0.05 && FW>=100 && BW>=100)
                                    {
                                        char key_for_max[50];
                                        memset(key_for_max,0,50);
                                        sprintf(key_for_max,"%s_%s_T",chrom,position);

                                        got_Germline_M_Hash=Germline_M_Hash.find(key_for_max);
                                        if(got_Germline_M_Hash==Germline_M_Hash.end())
                                        {
                                            double value=0;
                                            Germline_M_Hash[key_for_max]=value;
                                        }
                                        else
                                        {
                                            double value=got_Germline_M_Hash->second;
                                            if(value<=AF_T)
                                            {
                                                Germline_M_Hash[key_for_max]=AF_T;
                                            }
                                        }
                                    }
                            }      
                    }
                }
            }
        }
        fclose(input);
        //std::cout<<"\t\tParsing: "<<current_file<<" Storing:"<<Value_Hash_1.size()<<std::endl;
        //write a message every 50 samples
        if(count % 50==0)
        {
          std::cout<<" \tParsed successfully "<<ANSI_COLOR_GREEN<<count<<"/"<<FILE_Hash.size()<<ANSI_COLOR_RESET<<" germline samples and stored "<<ANSI_COLOR_GREEN<<Value_Hash_Thres.size()<<ANSI_COLOR_RESET<<" germline records"<<std::endl;
        }

    }
}

//function that computes thresholds given germline statistics
void estimateThresholds(float norm_factor,std::unordered_map<std::string,std::string> &Position_Hash,std::unordered_multimap<std::string,std::string> &Value_Hash,std::unordered_map<std::string,std::string> &Thresholds,std::unordered_map<std::string,std::string> &CountHash,std::unordered_map<std::string,std::string> &RatioHash)
{

    char chrom_chrX[]="chrX";
    
    char chrom_X[]="X";

    char chrom_chrY[]="chrY";
    char chrom_Y[]="Y";

    char chrom_chrM[]="chrM";
    char chrom_M[]="M";
    

  char insert_key[50];
  char insert_value[50];
  char count_value[50];
  char coverage_value[50];
  char prompt_key[50];
  char lineCharArray[MINIMUM_READ_LENGTH];
  std::pair <std::unordered_multimap<std::string,std::string>::iterator, std::unordered_multimap<std::string,std::string>::iterator> ret1;
  //read counts in the forward
  char base_fw_char[50];
  //read counts in the reverse
  char base_bw_char[50];
  //total RD forward
  char RD_fw_char[50];
  //total RD reverse
  char RD_bw_char[50];
  //store the integer values here
  int base_fw=-1;
  int base_bw=-1;
  int RD_fw=-1;
  int RD_bw=-1;
  //save the sum of counts per nucleotide and the total RD
  double sum_fw_nt=-1;
  double sum_fw_RD=-1;
  double sum_bw_nt=-1;
  double sum_bw_RD=-1;

  //add two variables for estimating the ratio
  int sum_Ri_fw=0;
  int sum_Ri_bw=0;

  double AF_fw=0;
  double AF_bw=0;

  int count=0;
  char chrom[50];
  char position[50];
  char base[50];
  //scan all reference positions
    for(std::unordered_map<std::string,std::string>::iterator a=Position_Hash.begin(); a!=Position_Hash.end();++a)
    {

        
      //in Positions table, the prompt key is chr_position BUT in the Germline statistics table all records have chr_position_A or C or G or T
      //so this is what we need to generate and promot
      
        //work with A
        memset(prompt_key,0,50);
        sprintf(prompt_key,"%s_A",a->first.c_str());

        //retrieve the chromosome name
        
        memset(chrom,0,50);
        memset(position,0,50);
        memset(base,0,50);
        sscanf(prompt_key,"%[^_]_%[^_]_%[^_]",chrom,position,base);


        ret1 = Value_Hash.equal_range(prompt_key);
        sum_fw_nt=0;
        sum_fw_RD=0;
        sum_bw_nt=0;
        sum_bw_RD=0;
        count=0;
        //add this to print the Ri/RD ratio
        sum_Ri_fw=0;
        sum_Ri_bw=0;

        for(std::unordered_multimap<std::string,std::string>::iterator b=ret1.first; b!=ret1.second; ++b)
        {
          char tmp[50];
          memset(tmp,0,50);
          //store the value
          sprintf(tmp,"%s",b->second.c_str());
          //and analyze it
          memset(base_fw_char,0,50);
          memset(base_bw_char,0,50);
          memset(RD_fw_char,0,50);
          memset(RD_bw_char,0,50);
          base_fw=-1;
          base_bw=-1;
          RD_fw=-1;
          RD_bw=-1;
          sscanf(tmp,"%[^_]_%[^_]_%[^_]_%[^_]",base_fw_char,RD_fw_char,base_bw_char,RD_bw_char);
          base_fw=std::stoi(base_fw_char);
          base_bw=std::stoi(base_bw_char);
          RD_fw=std::stoi(RD_fw_char);
          RD_bw=std::stoi(RD_bw_char);

          //this part of the code is very important because it changes the way we estimate the thresholds
          //for chromosomes with one copy or two copies
          if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
          {
                //we assume that everything below 5% is noise
               //double AF_limit=0.78;
               double AF_limit=0.05;
               double AF_fw=float(base_fw)/float(RD_fw);
               double AF_bw=float(base_bw)/float(RD_bw);
               if(AF_fw<=AF_limit && AF_bw<=AF_limit && RD_fw>=100 && RD_bw>=100)
               {
                sum_fw_nt=sum_fw_nt + base_fw + (float(RD_fw)*float(norm_factor));
                sum_fw_RD=sum_fw_RD + RD_fw;
                sum_bw_nt=sum_bw_nt + base_bw + (float(RD_bw)*float(norm_factor));
                sum_bw_RD=sum_bw_RD + RD_bw;

                //compute the values we need
                sum_Ri_fw=sum_Ri_fw+base_fw;
                sum_Ri_bw=sum_Ri_bw+base_bw;

                count++;
               }
          }
          else
          {
               //we assume that everything below 5% is noise
               double AF_limit=0.05;
               double AF_fw=float(base_fw)/float(RD_fw);
               double AF_bw=float(base_bw)/float(RD_bw);
               if(AF_fw<=AF_limit && AF_bw<=AF_limit && RD_fw>=100 && RD_bw>=100)
               {
                    sum_fw_nt=sum_fw_nt + base_fw + (float(RD_fw)*float(norm_factor));
                    sum_fw_RD=sum_fw_RD + RD_fw;
                    sum_bw_nt=sum_bw_nt + base_bw + (float(RD_bw)*float(norm_factor));
                    sum_bw_RD=sum_bw_RD + RD_bw;

                    //compute the values we need
                    sum_Ri_fw=sum_Ri_fw+base_fw;
                    sum_Ri_bw=sum_Ri_bw+base_bw;

                    count++;
               }
          }

          
        }

        if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
        {
               
               float ratio_fw=-888;
               float ratio_bw=-888;
               if (sum_fw_RD==0 && sum_bw_RD==0 )
               {
                 ratio_fw=-1;
                 ratio_bw=-1;
               }
               else if(sum_fw_RD==0 && sum_bw_RD!=0 )
               {
                 ratio_fw=-1;
                 ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               else if(sum_fw_RD!=0 && sum_bw_RD==0 )
               {
                 ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                 ratio_bw=-1;
               }
               else
               {
                    ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                    ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               //0.38 was the original value for count
               if(count<0.38*Value_Hash.count(prompt_key))
               {
                    //there is a problem here because there are a lot of sample that do not meet the requirements
                    //generate the record details for the forward strand
                    memset(insert_key,0,50);
                    sprintf(insert_key,"%s",prompt_key);
                    memset(insert_value,0,50);
                    sprintf(insert_value,"-1_-1");
                    //store it
                    Thresholds.insert(std::make_pair(insert_key,insert_value));  
                    memset(count_value,0,50);
                    sprintf(count_value,"%d",count);
                    CountHash.insert(std::make_pair(insert_key,count_value));

                    memset(coverage_value,0,50);
                    sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                    RatioHash.insert(std::make_pair(insert_key,coverage_value));
               }
               else
               {
                    AF_fw=float(sum_fw_nt)/float(sum_fw_RD);
                    AF_bw=float(sum_bw_nt)/float(sum_bw_RD);

                    if(isnan(AF_fw) || isnan(AF_bw))
                    {
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"-1_-1");
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"-1_-1");
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));
                    }
                    else
                    {
                         //generate the record details for the forward strand
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"%f_%f",AF_fw,AF_bw);
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value)); 

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));   
                    }      
               }
        }
        else
        {
               float ratio_fw=-888;
               float ratio_bw=-888;
               if (sum_fw_RD==0 && sum_bw_RD==0 )
               {
                 ratio_fw=-1;
                 ratio_bw=-1;
               }
               else if(sum_fw_RD==0 && sum_bw_RD!=0 )
               {
                 ratio_fw=-1;
                 ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               else if(sum_fw_RD!=0 && sum_bw_RD==0 )
               {
                 ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                 ratio_bw=-1;
               }
               else
               {
                    ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                    ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               //0.38 was the original value for count
               if(count<0.38*Value_Hash.count(prompt_key))
               {
                    //there is a problem here because there are a lot of sample that do not meet the requirements
                    //generate the record details for the forward strand
                    memset(insert_key,0,50);
                    sprintf(insert_key,"%s",prompt_key);
                    memset(insert_value,0,50);
                    sprintf(insert_value,"-1_-1");
                    //store it
                    Thresholds.insert(std::make_pair(insert_key,insert_value));  
                    memset(count_value,0,50);
                    sprintf(count_value,"%d",count);
                    CountHash.insert(std::make_pair(insert_key,count_value));

                    memset(coverage_value,0,50);
                    sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                    RatioHash.insert(std::make_pair(insert_key,coverage_value));
               }
               else
               {
                    AF_fw=float(sum_fw_nt)/float(sum_fw_RD);
                    AF_bw=float(sum_bw_nt)/float(sum_bw_RD);

                    if(isnan(AF_fw) || isnan(AF_bw))
                    {
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"-1_-1");
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"-1_-1");
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));
                    }
                    else
                    {
                         //generate the record details for the forward strand
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"%f_%f",AF_fw,AF_bw);
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value)); 

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));   
                    }      
               }
        }

        
        //work with C
        memset(prompt_key,0,50);
        sprintf(prompt_key,"%s_C",a->first.c_str());
        ret1 = Value_Hash.equal_range(prompt_key);
        sum_fw_nt=0;
        sum_fw_RD=0;
        sum_bw_nt=0;
        sum_bw_RD=0;
        count=0;

        //add this to print the Ri/RD ratio
        sum_Ri_fw=0;
        sum_Ri_bw=0;

        for(std::unordered_multimap<std::string,std::string>::iterator b=ret1.first; b!=ret1.second; ++b)
        {

          char tmp[50];
          memset(tmp,0,50);
          //store the value
          sprintf(tmp,"%s",b->second.c_str());
          //and analyze it
          memset(base_fw_char,0,50);
          memset(base_bw_char,0,50);
          memset(RD_fw_char,0,50);
          memset(RD_bw_char,0,50);
          base_fw=-1;
          base_bw=-1;
          RD_fw=-1;
          RD_bw=-1;
          sscanf(tmp,"%[^_]_%[^_]_%[^_]_%[^_]",base_fw_char,RD_fw_char,base_bw_char,RD_bw_char);
          base_fw=std::stoi(base_fw_char);
          base_bw=std::stoi(base_bw_char);
          RD_fw=std::stoi(RD_fw_char);
          RD_bw=std::stoi(RD_bw_char);
          //this part of the code is very important because it changes the way we estimate the thresholds
          //for chromosomes with one copy or two copies
          if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
          {
               double AF_limit=0.05;
               double AF_fw=float(base_fw)/float(RD_fw);
               double AF_bw=float(base_bw)/float(RD_bw);
               if(AF_fw<=AF_limit && AF_bw<=AF_limit && RD_fw>=100 && RD_bw>=100)
               {
                sum_fw_nt=sum_fw_nt + base_fw + (float(RD_fw)*float(norm_factor));
                sum_fw_RD=sum_fw_RD + RD_fw;
                sum_bw_nt=sum_bw_nt + base_bw + (float(RD_bw)*float(norm_factor));
                sum_bw_RD=sum_bw_RD + RD_bw;

                //compute the values we need
                sum_Ri_fw=sum_Ri_fw+base_fw;
                sum_Ri_bw=sum_Ri_bw+base_bw;

                count++;
               }
          }
          else
          {
               //we assume that everything below 5% is noise
               double AF_limit=0.05;
               double AF_fw=float(base_fw)/float(RD_fw);
               double AF_bw=float(base_bw)/float(RD_bw);
               if(AF_fw<=AF_limit && AF_bw<=AF_limit && RD_fw>=100 && RD_bw>=100)
               {
                    sum_fw_nt=sum_fw_nt + base_fw + (float(RD_fw)*float(norm_factor));
                    sum_fw_RD=sum_fw_RD + RD_fw;
                    sum_bw_nt=sum_bw_nt + base_bw + (float(RD_bw)*float(norm_factor));
                    sum_bw_RD=sum_bw_RD + RD_bw;

                    //compute the values we need
                    sum_Ri_fw=sum_Ri_fw+base_fw;
                    sum_Ri_bw=sum_Ri_bw+base_bw;

                    count++;
               }
          }
        }
        if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
        {
              
               float ratio_fw=-888;
               float ratio_bw=-888;
               if (sum_fw_RD==0 && sum_bw_RD==0 )
               {
                 ratio_fw=-1;
                 ratio_bw=-1;
               }
               else if(sum_fw_RD==0 && sum_bw_RD!=0 )
               {
                 ratio_fw=-1;
                 ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               else if(sum_fw_RD!=0 && sum_bw_RD==0 )
               {
                 ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                 ratio_bw=-1;
               }
               else
               {
                    ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                    ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               //0.38 was the original value for count
               if(count<0.38*Value_Hash.count(prompt_key))
               {
                    //there is a problem here because there are a lot of sample that do not meet the requirements
                    //generate the record details for the forward strand
                    memset(insert_key,0,50);
                    sprintf(insert_key,"%s",prompt_key);
                    memset(insert_value,0,50);
                    sprintf(insert_value,"-1_-1");
                    //store it
                    Thresholds.insert(std::make_pair(insert_key,insert_value));  
                    memset(count_value,0,50);
                    sprintf(count_value,"%d",count);
                    CountHash.insert(std::make_pair(insert_key,count_value));

                    memset(coverage_value,0,50);
                    sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                    RatioHash.insert(std::make_pair(insert_key,coverage_value));
               }
               else
               {
                    AF_fw=float(sum_fw_nt)/float(sum_fw_RD);
                    AF_bw=float(sum_bw_nt)/float(sum_bw_RD);

                    if(isnan(AF_fw) || isnan(AF_bw))
                    {
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"-1_-1");
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"-1_-1");
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));
                    }
                    else
                    {
                         //generate the record details for the forward strand
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"%f_%f",AF_fw,AF_bw);
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value)); 

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));   
                    }      
               }
        }
        else
        {
               float ratio_fw=-888;
               float ratio_bw=-888;
               if (sum_fw_RD==0 && sum_bw_RD==0 )
               {
                 ratio_fw=-1;
                 ratio_bw=-1;
               }
               else if(sum_fw_RD==0 && sum_bw_RD!=0 )
               {
                 ratio_fw=-1;
                 ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               else if(sum_fw_RD!=0 && sum_bw_RD==0 )
               {
                 ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                 ratio_bw=-1;
               }
               else
               {
                    ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                    ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               //0.38 was the original value for count
               if(count<0.38*Value_Hash.count(prompt_key))
               {
                    //there is a problem here because there are a lot of sample that do not meet the requirements
                    //generate the record details for the forward strand
                    memset(insert_key,0,50);
                    sprintf(insert_key,"%s",prompt_key);
                    memset(insert_value,0,50);
                    sprintf(insert_value,"-1_-1");
                    //store it
                    Thresholds.insert(std::make_pair(insert_key,insert_value));  
                    memset(count_value,0,50);
                    sprintf(count_value,"%d",count);
                    CountHash.insert(std::make_pair(insert_key,count_value));

                    memset(coverage_value,0,50);
                    sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                    RatioHash.insert(std::make_pair(insert_key,coverage_value)); 
               }
               else
               {
                    AF_fw=float(sum_fw_nt)/float(sum_fw_RD);
                    AF_bw=float(sum_bw_nt)/float(sum_bw_RD);
                    
                    if(isnan(AF_fw) || isnan(AF_bw))
                    {
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"-1_-1");
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"-1_-1");
                         RatioHash.insert(std::make_pair(insert_key,coverage_value)); 
                    }
                    else
                    {
                         //generate the record details for the forward strand
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"%f_%f",AF_fw,AF_bw);
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));  
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));  

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                         RatioHash.insert(std::make_pair(insert_key,coverage_value)); 
                    }      
               }
        }
        

        //work with G
        memset(prompt_key,0,50);
        sprintf(prompt_key,"%s_G",a->first.c_str());
        ret1 = Value_Hash.equal_range(prompt_key);
        sum_fw_nt=0;
        sum_fw_RD=0;
        sum_bw_nt=0;
        sum_bw_RD=0;

        //add this to print the Ri/RD ratio
        sum_Ri_fw=0;
        sum_Ri_bw=0;

        count=0;
        for(std::unordered_multimap<std::string,std::string>::iterator b=ret1.first; b!=ret1.second; ++b)
        {

          char tmp[50];
          memset(tmp,0,50);
          //store the value
          sprintf(tmp,"%s",b->second.c_str());
          //and analyze it
          memset(base_fw_char,0,50);
          memset(base_bw_char,0,50);
          memset(RD_fw_char,0,50);
          memset(RD_bw_char,0,50);
          base_fw=-1;
          base_bw=-1;
          RD_fw=-1;
          RD_bw=-1;
          sscanf(tmp,"%[^_]_%[^_]_%[^_]_%[^_]",base_fw_char,RD_fw_char,base_bw_char,RD_bw_char);
          base_fw=std::stoi(base_fw_char);
          base_bw=std::stoi(base_bw_char);
          RD_fw=std::stoi(RD_fw_char);
          RD_bw=std::stoi(RD_bw_char);
          //this part of the code is very important because it changes the way we estimate the thresholds
          //for chromosomes with one copy or two copies
          if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
          {
               //double AF_limit=0.78;
               double AF_limit=0.05;
               double AF_fw=float(base_fw)/float(RD_fw);
               double AF_bw=float(base_bw)/float(RD_bw);
               if(AF_fw<=AF_limit && AF_bw<=AF_limit && RD_fw>=100 && RD_bw>=100)
               {
                sum_fw_nt=sum_fw_nt + base_fw + (float(RD_fw)*float(norm_factor));
                sum_fw_RD=sum_fw_RD + RD_fw;
                sum_bw_nt=sum_bw_nt + base_bw + (float(RD_bw)*float(norm_factor));
                sum_bw_RD=sum_bw_RD + RD_bw;

                //compute the values we need
                sum_Ri_fw=sum_Ri_fw+base_fw;
                sum_Ri_bw=sum_Ri_bw+base_bw;

                count++;
               }
          }
          else
          {
               //we assume that everything below 5% is noise
               double AF_limit=0.05;
               double AF_fw=float(base_fw)/float(RD_fw);
               double AF_bw=float(base_bw)/float(RD_bw);
               if(AF_fw<=AF_limit && AF_bw<=AF_limit && RD_fw>=100 && RD_bw>=100)
               {
                    sum_fw_nt=sum_fw_nt + base_fw + (float(RD_fw)*float(norm_factor));
                    sum_fw_RD=sum_fw_RD + RD_fw;
                    sum_bw_nt=sum_bw_nt + base_bw + (float(RD_bw)*float(norm_factor));
                    sum_bw_RD=sum_bw_RD + RD_bw;

                    //compute the values we need
                    sum_Ri_fw=sum_Ri_fw+base_fw;
                    sum_Ri_bw=sum_Ri_bw+base_bw;

                    count++;
               }
          }
        }
        
        if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
        {
               
               float ratio_fw=-888;
               float ratio_bw=-888;
               if (sum_fw_RD==0 && sum_bw_RD==0 )
               {
                 ratio_fw=-1;
                 ratio_bw=-1;
               }
               else if(sum_fw_RD==0 && sum_bw_RD!=0 )
               {
                 ratio_fw=-1;
                 ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               else if(sum_fw_RD!=0 && sum_bw_RD==0 )
               {
                 ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                 ratio_bw=-1;
               }
               else
               {
                    ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                    ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               //0.38 was the original value for count
               if(count<0.38*Value_Hash.count(prompt_key))
               {
                    //there is a problem here because there are a lot of sample that do not meet the requirements
                    //generate the record details for the forward strand
                    memset(insert_key,0,50);
                    sprintf(insert_key,"%s",prompt_key);
                    memset(insert_value,0,50);
                    sprintf(insert_value,"-1_-1");
                    //store it
                    Thresholds.insert(std::make_pair(insert_key,insert_value));  
                    memset(count_value,0,50);
                    sprintf(count_value,"%d",count);
                    CountHash.insert(std::make_pair(insert_key,count_value));

                    memset(coverage_value,0,50);
                    sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                    RatioHash.insert(std::make_pair(insert_key,coverage_value));
               }
               else
               {
                    AF_fw=float(sum_fw_nt)/float(sum_fw_RD);
                    AF_bw=float(sum_bw_nt)/float(sum_bw_RD);

                    if(isnan(AF_fw) || isnan(AF_bw))
                    {
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"-1_-1");
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"-1_-1");
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));
                    }
                    else
                    {
                         //generate the record details for the forward strand
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"%f_%f",AF_fw,AF_bw);
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value)); 

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));   
                    }      
               }
        }
        else
        {   
               float ratio_fw=-888;
               float ratio_bw=-888;
               if (sum_fw_RD==0 && sum_bw_RD==0 )
               {
                 ratio_fw=-1;
                 ratio_bw=-1;
               }
               else if(sum_fw_RD==0 && sum_bw_RD!=0 )
               {
                 ratio_fw=-1;
                 ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               else if(sum_fw_RD!=0 && sum_bw_RD==0 )
               {
                 ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                 ratio_bw=-1;
               }
               else
               {
                    ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                    ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               //0.38 was the original value for count
               if(count<0.38*Value_Hash.count(prompt_key))
               {
                    //there is a problem here because there are a lot of sample that do not meet the requirements
                    //generate the record details for the forward strand
                    memset(insert_key,0,50);
                    sprintf(insert_key,"%s",prompt_key);
                    memset(insert_value,0,50);
                    sprintf(insert_value,"-1_-1");
                    //store it
                    Thresholds.insert(std::make_pair(insert_key,insert_value));  
                    memset(count_value,0,50);
                    sprintf(count_value,"%d",count);
                    CountHash.insert(std::make_pair(insert_key,count_value));

                    memset(coverage_value,0,50);
                    sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                    RatioHash.insert(std::make_pair(insert_key,coverage_value)); 
               }
               else
               {
                    AF_fw=float(sum_fw_nt)/float(sum_fw_RD);
                    AF_bw=float(sum_bw_nt)/float(sum_bw_RD);
                    if(isnan(AF_fw) || isnan(AF_bw))
                    {
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"-1_-1");
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value)); 
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"-1_-1");
                         RatioHash.insert(std::make_pair(insert_key,coverage_value)); 
                    }
                    else
                    {
                         //generate the record details for the forward strand
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"%f_%f",AF_fw,AF_bw);
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));  
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));  

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                         RatioHash.insert(std::make_pair(insert_key,coverage_value)); 
                    }      
               }
        }
        

        //work with T
        memset(prompt_key,0,50);
        sprintf(prompt_key,"%s_T",a->first.c_str());
        ret1 = Value_Hash.equal_range(prompt_key);
        sum_fw_nt=0;
        sum_fw_RD=0;
        sum_bw_nt=0;
        sum_bw_RD=0;

        //add this to print the Ri/RD ratio
        sum_Ri_fw=0;
        sum_Ri_bw=0;

        count=0;
        for(std::unordered_multimap<std::string,std::string>::iterator b=ret1.first; b!=ret1.second; ++b)
        {

          char tmp[50];
          memset(tmp,0,50);
          //store the value
          sprintf(tmp,"%s",b->second.c_str());
          //and analyze it
          memset(base_fw_char,0,50);
          memset(base_bw_char,0,50);
          memset(RD_fw_char,0,50);
          memset(RD_bw_char,0,50);
          base_fw=-1;
          base_bw=-1;
          RD_fw=-1;
          RD_bw=-1;
          sscanf(tmp,"%[^_]_%[^_]_%[^_]_%[^_]",base_fw_char,RD_fw_char,base_bw_char,RD_bw_char);
          base_fw=std::stoi(base_fw_char);
          base_bw=std::stoi(base_bw_char);
          RD_fw=std::stoi(RD_fw_char);
          RD_bw=std::stoi(RD_bw_char);
          //this part of the code is very important because it changes the way we estimate the thresholds
          //for chromosomes with one copy or two copies
          if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
          {
               //double AF_limit=0.78;
               double AF_limit=0.05;
               double AF_fw=float(base_fw)/float(RD_fw);
               double AF_bw=float(base_bw)/float(RD_bw);
               if(AF_fw<=AF_limit && AF_bw<=AF_limit && RD_fw>=100 && RD_bw>=100)
               {
                sum_fw_nt=sum_fw_nt + base_fw + (float(RD_fw)*float(norm_factor));
                sum_fw_RD=sum_fw_RD + RD_fw;
                sum_bw_nt=sum_bw_nt + base_bw + (float(RD_bw)*float(norm_factor));
                sum_bw_RD=sum_bw_RD + RD_bw;

                //compute the values we need
                sum_Ri_fw=sum_Ri_fw+base_fw;
                sum_Ri_bw=sum_Ri_bw+base_bw;

                count++;
               }
          }
          else
          {
               //we assume that everything below 5% is noise
               double AF_limit=0.05;
               double AF_fw=float(base_fw)/float(RD_fw);
               double AF_bw=float(base_bw)/float(RD_bw);
               if(AF_fw<=AF_limit && AF_bw<=AF_limit && RD_fw>=100 && RD_bw>=100)
               {
                    sum_fw_nt=sum_fw_nt + base_fw + (float(RD_fw)*float(norm_factor));
                    sum_fw_RD=sum_fw_RD + RD_fw;
                    sum_bw_nt=sum_bw_nt + base_bw + (float(RD_bw)*float(norm_factor));
                    sum_bw_RD=sum_bw_RD + RD_bw;

                    //compute the values we need
                    sum_Ri_fw=sum_Ri_fw+base_fw;
                    sum_Ri_bw=sum_Ri_bw+base_bw;

                    count++;
               }
          }
        }
        if(strcmp(chrom,chrom_X)==0 || strcmp(chrom,chrom_Y)==0 || strcmp(chrom,chrom_M)==0 || strcmp(chrom,chrom_chrX)==0 || strcmp(chrom,chrom_chrY)==0 || strcmp(chrom,chrom_chrM)==0)
        {
               
               float ratio_fw=-888;
               float ratio_bw=-888;
               if (sum_fw_RD==0 && sum_bw_RD==0 )
               {
                 ratio_fw=-1;
                 ratio_bw=-1;
               }
               else if(sum_fw_RD==0 && sum_bw_RD!=0 )
               {
                 ratio_fw=-1;
                 ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               else if(sum_fw_RD!=0 && sum_bw_RD==0 )
               {
                 ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                 ratio_bw=-1;
               }
               else
               {
                    ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                    ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               //0.38 was the original value for count
               if(count<0.38*Value_Hash.count(prompt_key))
               {
                    //there is a problem here because there are a lot of sample that do not meet the requirements
                    //generate the record details for the forward strand
                    memset(insert_key,0,50);
                    sprintf(insert_key,"%s",prompt_key);
                    memset(insert_value,0,50);
                    sprintf(insert_value,"-1_-1");
                    //store it
                    Thresholds.insert(std::make_pair(insert_key,insert_value));  
                    memset(count_value,0,50);
                    sprintf(count_value,"%d",count);
                    CountHash.insert(std::make_pair(insert_key,count_value));

                    memset(coverage_value,0,50);
                    sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                    RatioHash.insert(std::make_pair(insert_key,coverage_value));
               }
               else
               {
                    AF_fw=float(sum_fw_nt)/float(sum_fw_RD);
                    AF_bw=float(sum_bw_nt)/float(sum_bw_RD);

                    if(isnan(AF_fw) || isnan(AF_bw))
                    {
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"-1_-1");
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"-1_-1");
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));
                    }
                    else
                    {
                         //generate the record details for the forward strand
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"%f_%f",AF_fw,AF_bw);
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value)); 

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                         RatioHash.insert(std::make_pair(insert_key,coverage_value));   
                    }      
               }
        }
        else
        { 
               float ratio_fw=-888;
               float ratio_bw=-888;
               if (sum_fw_RD==0 && sum_bw_RD==0 )
               {
                 ratio_fw=-1;
                 ratio_bw=-1;
               }
               else if(sum_fw_RD==0 && sum_bw_RD!=0 )
               {
                 ratio_fw=-1;
                 ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               else if(sum_fw_RD!=0 && sum_bw_RD==0 )
               {
                 ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                 ratio_bw=-1;
               }
               else
               {
                    ratio_fw=float(sum_Ri_fw)/float(sum_fw_RD);
                    ratio_bw=float(sum_Ri_bw)/float(sum_bw_RD);
               }
               //0.38 was the original value for count
               if(count<0.38*Value_Hash.count(prompt_key))
               {
                    //there is a problem here because there are a lot of sample that do not meet the requirements
                    //generate the record details for the forward strand
                    memset(insert_key,0,50);
                    sprintf(insert_key,"%s",prompt_key);
                    memset(insert_value,0,50);
                    sprintf(insert_value,"-1_-1");
                    //store it
                    Thresholds.insert(std::make_pair(insert_key,insert_value));  
                    memset(count_value,0,50);
                    sprintf(count_value,"%d",count);
                    CountHash.insert(std::make_pair(insert_key,count_value));

                    memset(coverage_value,0,50);
                    sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                    RatioHash.insert(std::make_pair(insert_key,coverage_value)); 
               }
               else
               {
                    AF_fw=float(sum_fw_nt)/float(sum_fw_RD);
                    AF_bw=float(sum_bw_nt)/float(sum_bw_RD);
                    if(isnan(AF_fw) || isnan(AF_bw))
                    {
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"-1_-1");
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value));

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                         RatioHash.insert(std::make_pair(insert_key,coverage_value)); 
                    }
                    else
                    {
                         //generate the record details for the forward strand
                         memset(insert_key,0,50);
                         sprintf(insert_key,"%s",prompt_key);
                         memset(insert_value,0,50);
                         sprintf(insert_value,"%f_%f",AF_fw,AF_bw);
                         //store it
                         Thresholds.insert(std::make_pair(insert_key,insert_value));   
                         memset(count_value,0,50);
                         sprintf(count_value,"%d",count);
                         CountHash.insert(std::make_pair(insert_key,count_value)); 

                         memset(coverage_value,0,50);
                         sprintf(coverage_value,"%f_%f",ratio_fw,ratio_bw);
                         RatioHash.insert(std::make_pair(insert_key,coverage_value)); 
                    }      
               }
        }
        
        
    }
    std::cout<<ANSI_COLOR_GREEN<< Thresholds.size() <<ANSI_COLOR_RESET<<" thresholds computed successfully" <<std::endl;
    //printThresholds(Thresholds);
}

void generateFinalOutput(float C_value_float,char *panelDesign,std::unordered_map<std::string,std::string> &Reference_Hash,std::unordered_map<std::string,std::string> &Duplicate_Hash,std::unordered_map<std::string,std::string> &Thresholds,std::unordered_map<std::string,double> &Germline_M_Hash,char *output_dir,std::unordered_map<std::string,std::string> &CountHash,std::unordered_map<std::string,std::string> &RatioHash)
{

    char lineCharArray[MINIMUM_READ_LENGTH];
    FILE *input;
    int ret=-1;

    //write the output
    char outputFile[MINIMUM_READ_LENGTH];
    memset(outputFile,0,MINIMUM_READ_LENGTH);
    sprintf(outputFile,"%s/positionSpecificNoise_%.4f.txt",output_dir,C_value_float);
    std::ofstream output;
    output.open(outputFile);

    //output<<"chrom\tposition\treference\tduplicate\tThres_A\tThres_C\tThres_G\tThres_T\tGerm_Max_A\tGerm_Max_C\tGerm_Max_G\tGerm_Max_T\tCount_A\tCount_C\tCount_G\tCount_T\tRatio_A\tRatio_C\tRatio_G\tRatio_T"<<std::endl;
    output<<"chrom\tposition\treference\tduplicate\tThres_A\tThres_C\tThres_G\tThres_T\tGerm_Max_A\tGerm_Max_C\tGerm_Max_G\tGerm_Max_T"<<std::endl;


    //check if file exists
    if((input=fopen(panelDesign,"r"))==NULL)
    {
        printf("Error from function generateFinalOutput: Cannot open %s\n",panelDesign);
        exit(0);
    }
    else
    {
        //open the file
        input=fopen(panelDesign,"r");
        
        while(!feof(input))
        {
            //read line by line
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            //then get the ranges and extract the positions
            char chrom[50];
            char start[50];
            char end[50];
            char ampli_name[50];
            char rs_number[50];
            char gene_name[50];
            if(ret!=EOF)
            {
                memset(chrom,0,50);
                memset(start,0,50);
                memset(end,0,50);
                memset(ampli_name,0,50);
                memset(rs_number,0,50);
                memset(gene_name,0,50);
                sscanf(lineCharArray,"%s\t%s\t%s\t%s\t%s\t%s",chrom,start,end,ampli_name,rs_number,gene_name);
                //debug
                //cout<<chrom<<":"<<start<<":"<<end<<endl;
                
                //at this point we have the ranges and we need to extract all single positions in the ranges
                //this will be done as follows
                int start_num=-1;
                int end_num=-1;
                int idx=-1;
                start_num=std::stoi(start);
                end_num=std::stoi(end);
                for(idx=start_num;idx<=end_num;idx++)
                {
                    

                    char ref_prompt_key[50];
                    memset(ref_prompt_key,0,50);
                    sprintf(ref_prompt_key,"%s_%d",chrom,idx);

                    //use this key to find reference and duplicate

                    got_ReferenceBase_Hash=Reference_Hash.find(ref_prompt_key);
                    if(got_ReferenceBase_Hash==Reference_Hash.end())
                    {
                        std::cout<<"Problem found..: "<<std::endl;
                    }
                    else
                    {

                        output<<chrom<<"\t"<<idx<<"\t"<<got_ReferenceBase_Hash->second;
                    }

                    got_DuplicatePosition_Hash=Duplicate_Hash.find(ref_prompt_key);
                    if(got_DuplicatePosition_Hash==Duplicate_Hash.end())
                    {
                        output<<"\t"<<"NO";
                    }
                    else
                    {
                        output<<"\t"<<"YES";
                    }

                    char A_prompt_key[50];
                    memset(A_prompt_key,0,50);
                    sprintf(A_prompt_key,"%s_%d_A",chrom,idx);

                    char C_prompt_key[50];
                    memset(C_prompt_key,0,50);
                    sprintf(C_prompt_key,"%s_%d_C",chrom,idx);

                    char G_prompt_key[50];
                    memset(G_prompt_key,0,50);
                    sprintf(G_prompt_key,"%s_%d_G",chrom,idx);

                    char T_prompt_key[50];
                    memset(T_prompt_key,0,50);
                    sprintf(T_prompt_key,"%s_%d_T",chrom,idx);

                    //threshold A
                    got_Thresholds_Hash_Analytic=Thresholds.find(A_prompt_key);
                    if(got_Thresholds_Hash_Analytic==Thresholds.end())
                    {
                        std::cout<<"Problem with base 	 A"<<std::endl;
                    }
                    else
                    {

                        //this part of the code fixes the issues with the error estimation

                        char myRef[50];
                        memset(myRef,0,50);
                        sprintf(myRef,"%s",got_ReferenceBase_Hash->second.c_str());

                        char myRefFlag[]="A";
                        //the base is A so write -2_-2
                        if(strcmp(myRefFlag,myRef)==0)
                        {
                          output<<"\t-2_-2";
                        }
                        else
                        {
                          char myValue[50];
                          memset(myValue,0,50);
                          sprintf(myValue,"%s",got_Thresholds_Hash_Analytic->second.c_str());

                          char myValueFlag[]="-1_-1";
                          if(strcmp(myValueFlag,myValue)==0)
                          {
                            output<<"\t0.01_0.01";
                          }
                          else
                          {
                            output<<"\t"<<got_Thresholds_Hash_Analytic->second;
                          }

                        }

                        
                    }

                    //threshold C
                    got_Thresholds_Hash_Analytic=Thresholds.find(C_prompt_key);
                    if(got_Thresholds_Hash_Analytic==Thresholds.end())
                    {
                        std::cout<<"Problem with base 	 C"<<std::endl;
                    }
                    else
                    {
                        char myRef[50];
                        memset(myRef,0,50);
                        sprintf(myRef,"%s",got_ReferenceBase_Hash->second.c_str());

                        char myRefFlag[]="C";
                        //the base is A so write -2_-2
                        if(strcmp(myRefFlag,myRef)==0)
                        {
                          output<<"\t-2_-2";
                        }
                        else
                        {
                          char myValue[50];
                          memset(myValue,0,50);
                          sprintf(myValue,"%s",got_Thresholds_Hash_Analytic->second.c_str());

                          char myValueFlag[]="-1_-1";
                          if(strcmp(myValueFlag,myValue)==0)
                          {
                            output<<"\t0.01_0.01";
                          }
                          else
                          {
                            output<<"\t"<<got_Thresholds_Hash_Analytic->second;
                          }

                        }
                    }

                    //threshold G
                    got_Thresholds_Hash_Analytic=Thresholds.find(G_prompt_key);
                    if(got_Thresholds_Hash_Analytic==Thresholds.end())
                    {
                        std::cout<<"Problem with base 	 G"<<std::endl;
                    }
                    else
                    {
                        char myRef[50];
                        memset(myRef,0,50);
                        sprintf(myRef,"%s",got_ReferenceBase_Hash->second.c_str());

                        char myRefFlag[]="G";
                        //the base is A so write -2_-2
                        if(strcmp(myRefFlag,myRef)==0)
                        {
                          output<<"\t-2_-2";
                        }
                        else
                        {
                          char myValue[50];
                          memset(myValue,0,50);
                          sprintf(myValue,"%s",got_Thresholds_Hash_Analytic->second.c_str());

                          char myValueFlag[]="-1_-1";
                          if(strcmp(myValueFlag,myValue)==0)
                          {
                            output<<"\t0.01_0.01";
                          }
                          else
                          {
                            output<<"\t"<<got_Thresholds_Hash_Analytic->second;
                          }

                        }
                    }

                    //threshold T
                    got_Thresholds_Hash_Analytic=Thresholds.find(T_prompt_key);
                    if(got_Thresholds_Hash_Analytic==Thresholds.end())
                    {
                        std::cout<<"Problem with base 	 T"<<std::endl;
                    }
                    else
                    {
                        char myRef[50];
                        memset(myRef,0,50);
                        sprintf(myRef,"%s",got_ReferenceBase_Hash->second.c_str());

                        char myRefFlag[]="T";
                        //the base is A so write -2_-2
                        if(strcmp(myRefFlag,myRef)==0)
                        {
                          output<<"\t-2_-2";
                        }
                        else
                        {
                          char myValue[50];
                          memset(myValue,0,50);
                          sprintf(myValue,"%s",got_Thresholds_Hash_Analytic->second.c_str());

                          char myValueFlag[]="-1_-1";
                          if(strcmp(myValueFlag,myValue)==0)
                          {
                            output<<"\t0.01_0.01";
                          }
                          else
                          {
                            output<<"\t"<<got_Thresholds_Hash_Analytic->second;
                          }

                        }
                    }

                    //max A
                    got_Germline_Max_Hash=Germline_M_Hash.find(A_prompt_key);
                    if(got_Germline_Max_Hash==Germline_M_Hash.end())
                    {
                        //std::cout<<"Problem with base 	 max A"<<std::endl;
                        output<<"\t-";
                    }
                    else
                    {
                        output<<"\t"<<got_Germline_Max_Hash->second;
                    }
                    //max C
                    got_Germline_Max_Hash=Germline_M_Hash.find(C_prompt_key);
                    if(got_Germline_Max_Hash==Germline_M_Hash.end())
                    {
                        //std::cout<<"Problem with base 	 max C"<<std::endl;
                        output<<"\t-";
                    }
                    else
                    {
                        output<<"\t"<<got_Germline_Max_Hash->second;
                    }
                    //max G
                    got_Germline_Max_Hash=Germline_M_Hash.find(G_prompt_key);
                    if(got_Germline_Max_Hash==Germline_M_Hash.end())
                    {
                        //std::cout<<"Problem with base 	 max G"<<std::endl;
                        output<<"\t-";
                    }
                    else
                    {
                        output<<"\t"<<got_Germline_Max_Hash->second;
                    }
                    //max T
                    got_Germline_Max_Hash=Germline_M_Hash.find(T_prompt_key);
                    if(got_Germline_Max_Hash==Germline_M_Hash.end())
                    {
                        //std::cout<<"Problem with base 	 max T"<<std::endl;
                        output<<"\t-";
                    }
                    else
                    {
                        output<<"\t"<<got_Germline_Max_Hash->second;
                    }
                    //count A
                    //got_Count_Hash=CountHash.find(A_prompt_key);
                    //if(got_Count_Hash==CountHash.end())
                    //{
                        //std::cout<<"Problem with base 	 count A"<<std::endl;
                        //output<<"\t-";
                    //}
                    //else
                    //{
                        //output<<"\t"<<got_Count_Hash->second;
                    //}
                    //count C
                    //got_Count_Hash=CountHash.find(C_prompt_key);
                    //if(got_Count_Hash==CountHash.end())
                    //{
                        //std::cout<<"Problem with base 	 count C"<<std::endl;
                        //output<<"\t-";
                    //}
                    //else
                    //{
                        //output<<"\t"<<got_Count_Hash->second;
                    //}
                    //count G
                    //got_Count_Hash=CountHash.find(G_prompt_key);
                    //if(got_Count_Hash==CountHash.end())
                    //{
                        //std::cout<<"Problem with base 	 count G"<<std::endl;
                        //output<<"\t-";
                    //}
                    //else
                    //{
                        //output<<"\t"<<got_Count_Hash->second;
                    //}
                    //count T
                    //got_Count_Hash=CountHash.find(T_prompt_key);
                    //if(got_Count_Hash==CountHash.end())
                    //{
                        //std::cout<<"Problem with base 	 count T"<<std::endl;
                        //output<<"\t-";
                    //}
                    //else
                    //{
                        //output<<"\t"<<got_Count_Hash->second;
                    //}
                    //ratio A
                    //got_Ratio_Hash=RatioHash.find(A_prompt_key);
                    //if(got_Ratio_Hash==RatioHash.end())
                    //{
                        //std::cout<<"Problem with base 	 ratio A"<<std::endl;
                        //output<<"\t-";
                    //}
                    //else
                    //{
                        //output<<"\t"<<got_Ratio_Hash->second;
                    //}
                    //ratio C
                    //got_Ratio_Hash=RatioHash.find(C_prompt_key);
                    //if(got_Ratio_Hash==RatioHash.end())
                    //{
                        //std::cout<<"Problem with base 	 ratio C"<<std::endl;
                        //output<<"\t-";
                    //}
                    //else
                    //{
                        //output<<"\t"<<got_Ratio_Hash->second;
                    //}
                    //ratio G
                    //got_Ratio_Hash=RatioHash.find(G_prompt_key);
                    //if(got_Ratio_Hash==RatioHash.end())
                    //{
                        //std::cout<<"Problem with base 	 ratio G"<<std::endl;
                        //output<<"\t-";
                    //}
                    //else
                    //{
                        //output<<"\t"<<got_Ratio_Hash->second;
                    //}
                    //ratio T
                    //got_Ratio_Hash=RatioHash.find(T_prompt_key);
                    //if(got_Ratio_Hash==RatioHash.end())
                    //{
                        //std::cout<<"Problem with base 	 ratio T"<<std::endl;
                        //output<<"\t-";
                    //}
                    //else
                    //{
                        //output<<"\t"<<got_Ratio_Hash->second;
                    //}
                    output<<std::endl;
                }
            }
        }
    }

}

//this is a default function that assigns the same error per position

void generateFinalOutput_default(float C_value_float,char *panelDesign,std::unordered_map<std::string,std::string> &Reference_Hash,std::unordered_map<std::string,std::string> &Duplicate_Hash,char *output_dir)
{

    char lineCharArray[MINIMUM_READ_LENGTH];
    FILE *input;
    int ret=-1;

    //write the output
    char outputFile[MINIMUM_READ_LENGTH];
    memset(outputFile,0,MINIMUM_READ_LENGTH);
    sprintf(outputFile,"%s/positionSpecificNoise_default.txt",output_dir);
    std::ofstream output;
    output.open(outputFile);

    output<<"chrom\tposition\treference\tduplicate\tThres_A\tThres_C\tThres_G\tThres_T\tGerm_Max_A\tGerm_Max_C\tGerm_Max_G\tGerm_Max_T"<<std::endl;

    //check if file exists
    if((input=fopen(panelDesign,"r"))==NULL)
    {
        printf("Error from function generateFinalOutput: Cannot open %s\n",panelDesign);
        exit(0);
    }
    else
    {
        //open the file
        input=fopen(panelDesign,"r");
        
        while(!feof(input))
        {
            //read line by line
            memset(lineCharArray,0,MINIMUM_READ_LENGTH);
            ret=fscanf(input,"%1000[^\n]\n",lineCharArray);
            //then get the ranges and extract the positions
            char chrom[50];
            char start[50];
            char end[50];
            char ampli_name[50];
            char rs_number[50];
            char gene_name[50];
            if(ret!=EOF)
            {
                memset(chrom,0,50);
                memset(start,0,50);
                memset(end,0,50);
                memset(ampli_name,0,50);
                memset(rs_number,0,50);
                memset(gene_name,0,50);
                sscanf(lineCharArray,"%s\t%s\t%s\t%s\t%s\t%s",chrom,start,end,ampli_name,rs_number,gene_name);
                //debug
                //cout<<chrom<<":"<<start<<":"<<end<<endl;
                
                //at this point we have the ranges and we need to extract all single positions in the ranges
                //this will be done as follows
                int start_num=-1;
                int end_num=-1;
                int idx=-1;
                start_num=std::stoi(start);
                end_num=std::stoi(end);
                for(idx=start_num;idx<=end_num;idx++)
                {
                    char ref_prompt_key[50];
                    memset(ref_prompt_key,0,50);
                    sprintf(ref_prompt_key,"%s_%d",chrom,idx);

                    //use this key to find reference and duplicate

                    got_ReferenceBase_Hash=Reference_Hash.find(ref_prompt_key);
                    if(got_ReferenceBase_Hash==Reference_Hash.end())
                    {
                        std::cout<<"Problem found.."<<std::endl;
                    }
                    else
                    {

                        output<<chrom<<"\t"<<idx<<"\t"<<got_ReferenceBase_Hash->second;
                    }

                    got_DuplicatePosition_Hash=Duplicate_Hash.find(ref_prompt_key);
                    if(got_DuplicatePosition_Hash==Duplicate_Hash.end())
                    {
                        output<<"\t"<<"NO";
                    }
                    else
                    {
                        output<<"\t"<<"YES";
                    }
                    char value[50];
                    memset(value,0,50);
                    //sprintf(value,"%.4f_%.4f",C_value_float,C_value_float);
                    sprintf(value,"0.01_0.01");
                    output<<"\t"<<value<<"\t"<<value<<"\t"<<value<<"\t"<<value<<"\t-\t-\t-\t-"<<std::endl;
                }
            }
        }
    }
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

//generate a folder
void generateFolder(char *dir_path)
{
    char command[MINIMUM_READ_LENGTH];
    memset(command,0,MINIMUM_READ_LENGTH);
    sprintf(command, "mkdir -p %s",dir_path);
    myExec(command);
        
}

