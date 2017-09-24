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
     Email to dimitrios.kleftogiannis@kicr.ac.uk 
     I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
     You are free to modify, extend or distribute this code, as long 
     as this copyright notice is included whole and unchanged. 

END COPYRIGHT NOTICE
 
UTILITY
  This program estimates the CN status given tumour-normal matched files.
  For the tumour sample the output contains the ratio X/Y per amplicon of
  interest where X is the coverage and Y is the total number of reads in the sample.
  This information can be used with the output of AmplisolvePreProCN to visualize
  the amplification/deletions with respect to the pool of normal samples.


INPUT ARGUMENTS

    1. panel design                 : a tab limited bed-like file with info about the gene panel. Typically this file has 6 columns
                                    as follows: 
                                    chrom \t amplicon_start \t amplicon_end \t info \t info \t gene_name
                              
    2. tumour bam file              : a bam file for the tumour. Has to be indexed

    3. germline bam file            : a germline bam file that matches tumour. The file has to be indexed 
                                    If germline is not available type not_available 

    4. mapping quality Q            : mapping quality Q. Default is 20 for Ion Torrent platforms

    5. chromosome                   : the chromosome you want to estimate the CN status
                              
    6. output dir                   : directory for storing all results of execution. Becareful with the full paths
    

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

    cc -o  AmpliSolvePreProCN  AmpliSolvePreProCN.cpp -std=c++0x -lstdc++

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

//color the output
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[36m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"

//definitions of read buffers and other variables
#define MAXIMUM_READ_LENGTH 10000
#define MINIMUM_READ_LENGTH 1000

//functions
void printUsage();
void generateFolder(char *dir_path);
void myExec(char *command);
void storeBothPanels(char *amplicon_coordinates,std::map<int,std::string> &Hash1,std::map<int,std::string> &Hash2, char *CHROM,char *filtered_panel);
void computeNumberOfMappedReads(char *panel_design,std::map<std::string,std::string> &BAM_Hash,int quality,std::unordered_map<std::string,int> &RDS_Hash,char *result_dir);
void computeAmpliconCoverage(std::map<std::string,std::string> &BAM_Hash,int quality, std::unordered_map<std::string,int> &Cov_Hash, std::map<int,std::string> &Amplicon_Hash,char *result_dir);
void outputTumourRatios(char *output_name,std::map<std::string,std::string> &BAM_Hash,std::map<int,std::string> &Amplicon_Hash,std::unordered_map<std::string,int> &Cov_Hash,std::unordered_map<std::string,int> &RDS_Hash);
void outTumourCNstatus(char *output_name,std::unordered_map<std::string,int> &Tumour_Ampl_Hash,std::unordered_map<std::string,int> &Tumour_RDS_Hash,std::unordered_map<std::string,int> &Germline_Ampl_Hash,std::unordered_map<std::string,int> &Germline_RDS_Hash);
int checkFileExistance(char *name);

//data structures
//store the file name and the actual paths of tumour and normal
std::map<std::string,std::string> BamFileList_Tumour_Hash;
std::map<std::string,std::string>::iterator got_BamFileList_Tumour_Hash;
std::map<std::string,std::string> BamFileList_Germline_Hash;
std::map<std::string,std::string>::iterator got_BamFileList_Germline_Hash;
//store the panel design
std::map<int,std::string> PanelDesign_Hash;
std::map<int,std::string>::iterator got_PanelDesign_Hash;
//store the filtered panel design
std::map<int,std::string> FilteredPanelDesign_Hash;
std::map<int,std::string>::iterator got_FilteredPanelDesign_Hash;
//store the number of reads mapped on all amplicons
std::unordered_map<std::string,int> TumourNumReadsAllPanel_Hash;
std::unordered_map<std::string,int>::iterator got_TumourNumReadsAllPanel_Hash;
//store the number of reads mapped on filtered amplicons
std::unordered_map<std::string,int> TumourNumReadsFilteredPanel_Hash;
std::unordered_map<std::string,int>::iterator got_TumourNumReadsFilteredPanel_Hash;
//store the number of reads mapped on filtered amplicons for germline
std::unordered_map<std::string,int> GermlineNumReadsFilteredPanel_Hash;
std::unordered_map<std::string,int>::iterator got_GermlineNumReadsFilteredPanel_Hash;
//store the number of reads mapped on specific amplicon of interest for tumour
std::unordered_map<std::string,int> TumourAmpliconNumReads_Hash;
std::unordered_map<std::string,int>::iterator got_TumourAmpliconNumReads_Hash;
//store the number of reads mapped on specific amplicon of interest for germline
std::unordered_map<std::string,int> GermlineAmpliconNumReads_Hash;
std::unordered_map<std::string,int>::iterator got_GermlineAmpliconNumReads_Hash;

int main(int argc, char **argv)
{
        //initialize the input arguments, read the inputs from the command line and check for correctness
        char *user_panel_design;
        user_panel_design=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        char *user_tumour_bam;
        user_tumour_bam=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        char *user_germline_bam;
        user_germline_bam=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        char *user_Q;
        user_Q=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        char *user_chromosome;
        user_chromosome=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
        char *user_output_dir;
        user_output_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

        if(argc!=7)
        {
            std::cout<<"************************************************************************************************************************************"<<std::endl;
            std::cout<<ANSI_COLOR_RED<<"                                        Your input arguments are not correct"<<ANSI_COLOR_RESET<<std::endl;
            std::cout<<"                         Copyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n"<<std::endl;
            printUsage();
            return 0;
        }
        else
        {
            //execute commands
            char command[MINIMUM_READ_LENGTH];

            //store the input variable only if the number of input variables is correct   
            char *panel_design;
            panel_design=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            char *tumour_bam;
            tumour_bam=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            char *germline_bam;
            germline_bam=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            char *Q;
            Q=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            char *chromosome;
            chromosome=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));
            char *output_dir;
            output_dir=(char*)malloc(sizeof(char) * (MINIMUM_READ_LENGTH));

            strcpy(user_panel_design,argv[1]);
            memset(panel_design,0,MINIMUM_READ_LENGTH);
            sscanf(user_panel_design,"panel_design=%s",panel_design);

            strcpy(user_tumour_bam,argv[2]);
            memset(tumour_bam,0,MINIMUM_READ_LENGTH);
            sscanf(user_tumour_bam,"tumour_bam=%s",tumour_bam);

            strcpy(user_germline_bam,argv[3]);
            memset(germline_bam,0,MINIMUM_READ_LENGTH);
            sscanf(user_germline_bam,"germline_bam=%s",germline_bam);

            strcpy(user_Q,argv[4]);
            memset(Q,0,MINIMUM_READ_LENGTH);
            sscanf(user_Q,"Q=%s",Q);

            strcpy(user_chromosome,argv[5]);
            memset(chromosome,0,MINIMUM_READ_LENGTH);
            sscanf(user_chromosome,"chromosome=%s",chromosome);

            strcpy(user_output_dir,argv[6]);
            memset(output_dir,0,MINIMUM_READ_LENGTH);
            sscanf(user_output_dir,"output_dir=%s",output_dir);

            //flag the germline dir value
            char flag_germline_bam[]="not_available";
            int no_germline=-1;
            //save the Q in integer format
            int Q_int=std::atoi(Q);

            
            std::cout<<"************************************************************************************************************************************\n"<<std::endl;
            std::cout<<"                                AmpliSolveCN: Copy Number status and ratio estimation per chromosome\n"<<std::endl;
            std::cout<<"                        Copyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n"<<std::endl;
            std::cout<<"Execution started under the following parameters:"<<std::endl;
            std::cout<<"\t1. panel design                           : "<<ANSI_COLOR_GREEN<<panel_design<<ANSI_COLOR_RESET<<std::endl;
            
            int flag=checkFileExistance(tumour_bam);
            if(flag==1)
            {
                std::cout<<"\t2. tumour bam                             : "<<ANSI_COLOR_GREEN<<tumour_bam<<ANSI_COLOR_RESET<<std::endl;
                memset(command,0,MINIMUM_READ_LENGTH);
                //generate teh bam index. it helps...
                sprintf(command,"samtools index %s",tumour_bam);
                myExec(command);

            }
            else
            {
                    std::cout<<ANSI_COLOR_RED<<"                                         THE PROGRAM WILL TERMINATE: Tumour file does not exist    "<<ANSI_COLOR_RESET<<std::endl;
                    std::cout<<"************************************************************************************************************************************"<<std::endl;
                    exit(0);   
            }
            
            if(strcmp(germline_bam,flag_germline_bam)==0)
            {
                no_germline=-1;
                std::cout<<"\t3. germline bam                           : "<<ANSI_COLOR_RED<<"NO matched germline BAM available. Only the tumour ratio can be estimated."<<ANSI_COLOR_RESET<<std::endl;
            }
            else
            {
                flag=checkFileExistance(germline_bam);
                if(flag==1)
                {
                    std::cout<<"\t3. germline bam                           : "<<ANSI_COLOR_GREEN<<germline_bam<<ANSI_COLOR_RESET<<std::endl;
                    no_germline=1;
                    memset(command,0,MINIMUM_READ_LENGTH);
                    //generate the bam index. helps...
                    sprintf(command,"samtools index %s",germline_bam);
                    myExec(command);
                }
                else
                {
                    std::cout<<ANSI_COLOR_RED<<"                                         THE PROGRAM WILL TERMINATE: Germline file does not exist    "<<ANSI_COLOR_RESET<<std::endl;
                    std::cout<<"************************************************************************************************************************************"<<std::endl;
                    exit(0); 
                }
                
            }
            if(Q_int<=0)
            {
                Q_int=20;
                std::cout<<"\t4. Q                                    : "<<ANSI_COLOR_RED<<"User gave: "<<Q<<ANSI_COLOR_RESET<<" The value has changed to 20 as default"<<std::endl;
            }
            else
            {
                std::cout<<"\t4. Q                                      : "<<ANSI_COLOR_GREEN<<Q_int<<ANSI_COLOR_RESET<<std::endl;
            }
            std::cout<<"\t5. chromosome                             : "<<ANSI_COLOR_GREEN<<chromosome<<ANSI_COLOR_RESET<<std::endl;
            std::cout<<"\t6. output                                 : "<<ANSI_COLOR_GREEN<<output_dir<<ANSI_COLOR_RESET<<std::endl;

             //check if mathced germline bam is available
            if(no_germline==1)
            {
                //matched germline file is available so the output contains both
                //the X/Y ratio for the tumour and the CN status estimation given the ratio
                //X1/Y1 from the normal file

                //we have germlines so we continue...
                std::cout<<"Store matched tumour-normal files";
                //make a trick and get the actual file names of tumour and germline
                char tmp[200];
              
                //germline bam
                std::string str = germline_bam;
                std::string res = str.substr( str.find_last_of("/") + 1 );
                char filename_germline[200];
                memset(filename_germline,0,200);
                memset(tmp,0,200);
                sprintf(tmp,"%s",res.c_str());
                strncpy(filename_germline,tmp,strlen(tmp)-4);
                //tumour bam
                str = tumour_bam;
                res = str.substr( str.find_last_of("/") + 1 );
                char filename_tumour[200];
                memset(filename_tumour,0,200);
                memset(tmp,0,200);
                sprintf(tmp,"%s",res.c_str());
                strncpy(filename_tumour,tmp,strlen(tmp)-4);
                std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;

                //store filenames in Hash to be consistent with other implementations in Amplisolve
                //store the bam list in a Hash
                BamFileList_Tumour_Hash.insert(std::make_pair<std::string,std::string>(tumour_bam,filename_tumour));
                BamFileList_Germline_Hash.insert(std::make_pair<std::string,std::string>(germline_bam,filename_germline));

                //store the results of the program and the intermediate files
                char interm_results_dir[200];
                memset(interm_results_dir,0,200);
                //and the intermediate results
                sprintf(interm_results_dir,"%s/%s_CN_interm_files",output_dir,chromosome);
                generateFolder(interm_results_dir);

                //we store the full panel and the panel without the selected chromosome.
                //this is because we want to estimate the tumour ratio as we did in AmpliSolvePreProCN
                //and the CN status.

                //first lets compute the filtered panel filename
                char filtered_panel_design[MINIMUM_READ_LENGTH];
                memset(filtered_panel_design,0,MINIMUM_READ_LENGTH);
                sprintf(filtered_panel_design,"%s/no%s_panel.bed",interm_results_dir,chromosome);
                //generate it and flush it to drive
                //store the panel design and the filtered panel
                storeBothPanels(panel_design,PanelDesign_Hash,FilteredPanelDesign_Hash,chromosome,filtered_panel_design);
                //check if user's chromosome exist 
                int F=PanelDesign_Hash.size();
                if(F<=0)
                {
                    std::cout<<ANSI_COLOR_RED<<"                                         THE PROGRAM WILL TERMINATE: chromosome provided by user is not in the panel    "<<ANSI_COLOR_RESET<<std::endl;
                    std::cout<<"************************************************************************************************************************************"<<std::endl;
                    exit(0);
                }
                else
                {
                    //continue the execution from here

                    //store the number of reads on all amplicons for tumour
                    std::cout<<"Running function computeNumberOfMappedReads     \t";
                    computeNumberOfMappedReads(panel_design,BamFileList_Tumour_Hash,Q_int,TumourNumReadsAllPanel_Hash,interm_results_dir);
                    //store the number of reads on filtered amplicons for tumour
                    computeNumberOfMappedReads(filtered_panel_design,BamFileList_Tumour_Hash,Q_int,TumourNumReadsFilteredPanel_Hash,interm_results_dir);
                    //store the number of reads on filtered amplicons for germline
                    computeNumberOfMappedReads(filtered_panel_design,BamFileList_Germline_Hash,Q_int,GermlineNumReadsFilteredPanel_Hash,interm_results_dir);
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;

                    std::cout<<"Running function computeAmpliconCoverage     \t";
                    computeAmpliconCoverage(BamFileList_Tumour_Hash,Q_int,TumourAmpliconNumReads_Hash,PanelDesign_Hash,interm_results_dir);
                    computeAmpliconCoverage(BamFileList_Germline_Hash,Q_int,GermlineAmpliconNumReads_Hash,PanelDesign_Hash,interm_results_dir);
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;

                    //compute the simple ratio for tumour; write it in this file
                    char output_name1[MINIMUM_READ_LENGTH];
                    memset(output_name1,0,MINIMUM_READ_LENGTH);
                    sprintf(output_name1,"%s/%s_%s_ratios.txt",output_dir,filename_tumour,chromosome);
                    std::cout<<"Running function outputTumourRatios     \t";
                    outputTumourRatios(output_name1,BamFileList_Tumour_Hash,PanelDesign_Hash,TumourAmpliconNumReads_Hash,TumourNumReadsAllPanel_Hash);
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;

                    //here we compute the CN status usign the matched germline sample
                    char output_name2[MINIMUM_READ_LENGTH];
                    memset(output_name2,0,MINIMUM_READ_LENGTH);
                    sprintf(output_name2,"%s/%s_%s_CNstatus.txt",output_dir,filename_tumour,chromosome);
                     std::cout<<"Running function outTumourCNstatus     \t";
                    outTumourCNstatus(output_name2,TumourAmpliconNumReads_Hash,TumourNumReadsFilteredPanel_Hash,GermlineAmpliconNumReads_Hash,GermlineNumReadsFilteredPanel_Hash);
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;
                    std::cout<<"Execution was successful. Results can be found at: "<<ANSI_COLOR_YELLOW<<output_name1<<ANSI_COLOR_RESET<<" and "<<ANSI_COLOR_YELLOW<<output_name2<<ANSI_COLOR_RESET<<std::endl;
                }
            }
            else
            {
                //we dont have germline..so the output contains only the ratio X/Y for the tumour file
                //we have germlines so we continue...
                std::cout<<"Store the tumour file";
                //make a trick and get the actual file names of tumour and germline
                char tmp[200];

                //tumour bam
                std::string str = tumour_bam;
                std::string res = str.substr( str.find_last_of("/") + 1 );
                char filename_tumour[200];
                memset(filename_tumour,0,200);
                memset(tmp,0,200);
                sprintf(tmp,"%s",res.c_str());
                strncpy(filename_tumour,tmp,strlen(tmp)-4);
                std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;

                //store filenames in Hash to be consistent with other implementations in Amplisolve
                //store the tumour bam in a Hash
                BamFileList_Tumour_Hash.insert(std::make_pair<std::string,std::string>(tumour_bam,filename_tumour));

                //store the results of the program and the intermediate files
                char interm_results_dir[200];
                memset(interm_results_dir,0,200);
                //and the intermediate results
                sprintf(interm_results_dir,"%s/%s_CN_interm_files",output_dir,chromosome);
                generateFolder(interm_results_dir);

                //we store the full panel and the panel without the selected chromosome.
                //this is because we want to estimate the tumour ratio as we did in AmpliSolvePreProCN
                //and the CN status.

                //first lets compute the filtered panel filename
                char filtered_panel_design[MINIMUM_READ_LENGTH];
                memset(filtered_panel_design,0,MINIMUM_READ_LENGTH);
                sprintf(filtered_panel_design,"%s/no%s_panel.bed",interm_results_dir,chromosome);
                //generate it and flush it to drive
                //store the panel design and the filtered panel
                storeBothPanels(panel_design,PanelDesign_Hash,FilteredPanelDesign_Hash,chromosome,filtered_panel_design);
                //check if user's chromosome exist 
                int F=PanelDesign_Hash.size();
                if(F<=0)
                {
                    std::cout<<ANSI_COLOR_RED<<"                                         THE PROGRAM WILL TERMINATE: chromosome provided by the user is not in the panel    "<<ANSI_COLOR_RESET<<std::endl;
                    std::cout<<"************************************************************************************************************************************"<<std::endl;
                    exit(0);
                }
                else
                {
                    //continue the execution from here

                    //store the number of reads on all amplicons for tumour
                    std::cout<<"Running function computeNumberOfMappedReads     \t";
                    computeNumberOfMappedReads(panel_design,BamFileList_Tumour_Hash,Q_int,TumourNumReadsAllPanel_Hash,interm_results_dir);
                    //store the number of reads on filtered amplicons for tumour
                    computeNumberOfMappedReads(filtered_panel_design,BamFileList_Tumour_Hash,Q_int,TumourNumReadsFilteredPanel_Hash,interm_results_dir);
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;

                    std::cout<<"Running function computeAmpliconCoverage     \t";
                    computeAmpliconCoverage(BamFileList_Tumour_Hash,Q_int,TumourAmpliconNumReads_Hash,PanelDesign_Hash,interm_results_dir);
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;

                    //compute the simple ratio for tumour; write it in this file
                    char output_name1[MINIMUM_READ_LENGTH];
                    memset(output_name1,0,MINIMUM_READ_LENGTH);
                    sprintf(output_name1,"%s/%s_%s_ratios.txt",output_dir,filename_tumour,chromosome);
                    std::cout<<"Running function outputTumourRatios     \t";
                    outputTumourRatios(output_name1,BamFileList_Tumour_Hash,PanelDesign_Hash,TumourAmpliconNumReads_Hash,TumourNumReadsAllPanel_Hash);
                    std::cout<<ANSI_COLOR_GREEN<<"  OK"<<ANSI_COLOR_RESET<<std::endl;

                    std::cout<<"Execution was successful. Results can be found at: "<<ANSI_COLOR_YELLOW<<output_name1<<ANSI_COLOR_RESET<<std::endl;
                }
            }
                std::cout<<"************************************************************************************************************************************\n"<<std::endl;
                return 0;
        }
        //clear the memory before exiting
        BamFileList_Tumour_Hash.clear();
        BamFileList_Germline_Hash.clear();
        PanelDesign_Hash.clear();
        FilteredPanelDesign_Hash.clear();
        TumourNumReadsAllPanel_Hash.clear();
        TumourNumReadsFilteredPanel_Hash.clear();
        GermlineNumReadsFilteredPanel_Hash.clear();
        TumourAmpliconNumReads_Hash.clear();
        GermlineAmpliconNumReads_Hash.clear();
}

//print program's usage
void printUsage()
{
  std::cout<<"Please type the following: "<<std::endl;
    std::cout<<"\n./AmpliSolveCN"<<ANSI_COLOR_GREEN<<" panel_design="<<ANSI_COLOR_RESET<<"/your/panel/design/bed/file"<<ANSI_COLOR_GREEN<<" tumour_bam="<<ANSI_COLOR_RESET<<"/full/path/of/your/tumour/bam/file"<<ANSI_COLOR_GREEN<<" germline_bam="<<ANSI_COLOR_RESET<<"/full/path/of/your/germline/bam/file"<<ANSI_COLOR_GREEN<<" Q="<<ANSI_COLOR_RESET<<"mapping_quality_score"<<ANSI_COLOR_GREEN<<" chromosome="<<ANSI_COLOR_RESET<<"chr_to_estimate_CN"<<ANSI_COLOR_GREEN<<" output_dir="<<ANSI_COLOR_RESET<<"/dir/name/to/store/the/results/and/intermediate/results"<<std::endl;
    std::cout<<"\nRemember that: "<<std::endl;
    std::cout<<"This code internally uses SAMtools and Bedtools; make sure that are installed and configured (in your PATH)"<<std::endl;
    std::cout<<"\nIf no matched germline bam file is available just type germline_bam=not_available ; If germline files are available should be named bam and not BAM (case-sensitive) "<<std::endl;
    std::cout<<"\nExecution example:"<<std::endl;
    std::cout<<"The program takes 6 input arguments, so please fill them similar to the following command:"<<std::endl;
    std::cout<<"./AmpliSolveCN panel_design=/Users/dkleftog/Desktop/AmpliSolve_Execution_Example/myBED/example_panel_design.bed tumour_bam=/Users/dkleftog/Desktop/AmpliSolve_Execution_Example/Plasma_samples/Plasma_S1.bam germline_bam=/Users/dkleftog/Desktop/AmpliSolve_Execution_Example/Germline_samples/Germline_S1.bam Q=20 chromosome=chrX output_dir=/Users/dkleftog/Desktop/AmpliSolveTesting"<<std::endl;
    std::cout<<"\n\tIt is important to give the arguments in this order. Otherwise the program will crach!"<<std::endl;
    std::cout<<"\nMore info about the input data and detailed execution examples can be found at our web-repository"<<std::endl;
    std::cout<<"************************************************************************************************************************************"<<std::endl;
}

//generate a folder
void generateFolder(char *dir_path)
{
    char command[MINIMUM_READ_LENGTH];
    memset(command,0,MINIMUM_READ_LENGTH);
    sprintf(command, "mkdir -p %s",dir_path);
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

//store the chromosome of interest for the panel design and the rest for the filtered panel
//the function also flushes the filtered panel into the disk for further use
void storeBothPanels(char *amplicon_coordinates,std::map<int,std::string> &Hash1,std::map<int,std::string> &Hash2, char *CHROM,char *filtered_panel)
{

    std::ofstream output;
    output.open(filtered_panel);

    char lineCharArray[MAXIMUM_READ_LENGTH];
    FILE *input;
    int ret=-1;
    int count1=1;
    int count2=1;
    char complete_name[MAXIMUM_READ_LENGTH];
    sprintf(complete_name,"%s",amplicon_coordinates);
    if((input=fopen(complete_name,"r"))==NULL)
    {
        printf("Error from function storePanel: Cannot open %s\n",complete_name);
        exit(0);
    }
    else
    {
        //read the file line by line
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
                
                if(strcmp(CHROM,chrom)==0)
                {
                    //store only the ones you are interested
                    char value[50];
                    memset(value,0,50);
                    sprintf(value,"%s:%s-%s",chrom,pos_start,pos_end);
                    Hash1[count1]=value;
                    count1++;
                }
                else
                {
                    //store the filtered panel on disk
                    output<<lineCharArray<<std::endl;
                    //store the rest
                    char value[50];
                    memset(value,0,50);
                    sprintf(value,"%s:%s-%s",chrom,pos_start,pos_end);
                    Hash2[count2]=value;
                    count2++;

                }
            }
        }
    }
    fclose(input);
    output.close();
}

//this function runs internally samtools and bedtools and computes the total number of mapped reads on amplicons
void computeNumberOfMappedReads(char *panel_design,std::map<std::string,std::string> &BAM_Hash,int quality,std::unordered_map<std::string,int> &RDS_Hash,char *result_dir)
{
    //gerenate a random number unique for each execution
    time_t t;
    int seed;
    srand((unsigned) time(&t));
    seed=floor(rand()%1000);

    //save the command
    char command[MINIMUM_READ_LENGTH];
    //parse the files and generate a file with two columns: file name and number of reads mapped on all amplicons
    for(std::map<std::string,std::string>::iterator it=BAM_Hash.begin(); it!=BAM_Hash.end();++it)
    {
        
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
void computeAmpliconCoverage(std::map<std::string,std::string> &BAM_Hash,int quality, std::unordered_map<std::string,int> &Cov_Hash, std::map<int,std::string> &Amplicon_Hash,char *result_dir)
{
    //gerenate a random number unique for each execution
    time_t t;
    int seed;
    srand((unsigned) time(&t));
    seed=floor(rand()%1000);

    //save the command
    char command[MINIMUM_READ_LENGTH];
    //parse the files and compute for every amplicon the coverage per amplicon
    for(std::map<std::string,std::string>::iterator it=BAM_Hash.begin(); it!=BAM_Hash.end();++it)
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
                        //here the key is just the ID number of the amplicon for simplicity
                        sprintf(key,"%d",iit->first);
                        Cov_Hash.insert(std::pair<std::string,int>(key,rds));
                        count++;
                    }
                }
            }
            fclose(input);   
        }
    }

    //debuging: check the results
    //for(std::unordered_map<std::string,int>::iterator it=Cov_Hash.begin(); it!=Cov_Hash.end();++it)
    //{
        //std::cout<<it->first<<" and "<<it->second<<std::endl;
    //}

}

//this function prints the ratios for tumour file
void outputTumourRatios(char *output_name,std::map<std::string,std::string> &BAM_Hash,std::map<int,std::string> &Amplicon_Hash,std::unordered_map<std::string,int> &Cov_Hash,std::unordered_map<std::string,int> &RDS_Hash)
{

    //write the header of the file
    std::ofstream output;
    output.open(output_name);
    //write the header of the file
    output<<"Filename/SimpleRatio";
    for(std::map<int,std::string>::iterator it=Amplicon_Hash.begin(); it!=Amplicon_Hash.end();++it)
    {
        output<<"\t"<<it->second;
    }
    output<<std::endl;
    for(std::map<std::string,std::string>::iterator a=BAM_Hash.begin(); a!=BAM_Hash.end();++a)
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
            sprintf(prompt_key,"%d",b->first);
            //retrieve the number of reads on the amplicon
            got_TumourAmpliconNumReads_Hash=Cov_Hash.find(prompt_key);
            if(got_TumourAmpliconNumReads_Hash==Cov_Hash.end())
            {
                //something is wrong so write NaN
                output<<"\t-1";
            }
            else
            {
                //no obtain the number of reads mapped on all amplicons to generate the ratio
                got_TumourNumReadsAllPanel_Hash=RDS_Hash.find(filename);
                if(got_TumourNumReadsAllPanel_Hash==RDS_Hash.end())
                {
                    //something is wrong here so write NaN
                    output<<"\t-2";
                }
                else
                {
                    int rds=got_TumourNumReadsAllPanel_Hash->second;
                    if(rds==0)
                    {
                        //we cannot devide with 0 so write NaN
                        output<<"\t-3";
                    }
                    else
                    {
                        //we hae everything so we simply write the ratio
                        output<<"\t"<<float(got_TumourAmpliconNumReads_Hash->second)/float(rds);
                    }
                }
            }

        }
        output<<std::endl;
    }

    //here we compute the CN status which is similar to the previous computation BUT, we normalize the 
    //amplicon coverage by the number of reads mapped to all amplicons except for the reads in chromosome
    //of interest and the tumour ratio is further normalized by the germline ratio computed in the same way.

    output.close();
}


void outTumourCNstatus(char *output_name,std::unordered_map<std::string,int> &Tumour_Ampl_Hash,std::unordered_map<std::string,int> &Tumour_RDS_Hash,std::unordered_map<std::string,int> &Germline_Ampl_Hash,std::unordered_map<std::string,int> &Germline_RDS_Hash)
{
    //write the header of the file
    std::ofstream output;
    output.open(output_name);
    output<<"ID\tInfered_CN_status"<<std::endl;
    //get the number of reads for tumour
    int T=0;
    for(std::unordered_map<std::string,int>::iterator a=Tumour_RDS_Hash.begin(); a!=Tumour_RDS_Hash.end();++a)
    {
        T=a->second;
    }
    //get the number of reads for germline
    int N=0;
    for(std::unordered_map<std::string,int>::iterator a=Germline_RDS_Hash.begin(); a!=Germline_RDS_Hash.end();++a)
    {
        N=a->second;
    }
    //store the ratios
    std::multimap<float,std::string> Ratios;
    std::multimap<float,std::string>::iterator got_Ratios;
    //generate the ratio per amplicon
    char prompt_key[100];
    float ratio=0;
    float t=0;
    float n=0;
    for(std::unordered_map<std::string,int>::iterator it=Tumour_Ampl_Hash.begin(); it!=Tumour_Ampl_Hash.end();++it)
    {
        //key is the amplicon ID; it can be also the coordinate....
        memset(prompt_key,0,100);
        sprintf(prompt_key,"%s",it->first.c_str());
        //get the same value from the germline file...
        got_GermlineAmpliconNumReads_Hash=Germline_Ampl_Hash.find(prompt_key);
        if(got_GermlineAmpliconNumReads_Hash==Germline_Ampl_Hash.end())
        {
            std::cout<<"Kati malakia paizei edo"<<std::endl;
        }
        else
        {
            t=it->second;
            n=got_GermlineAmpliconNumReads_Hash->second;
        }

        if(n==0 || N==0 || T==0)
        {
            //skip this something is wrong give ratio -1
            ratio=-1;
        }
        else
        {
            //remember that the formula is
            //ratio=(t/T)/(n/N)
            ratio=float(t*N)/float(n*T);
            Ratios.insert(std::pair<float,std::string>(ratio,prompt_key));
        }
    }

     //estimate the quartiles. I use the sample formulas that are suitable for normal distribution
    int hash_size=Ratios.size();
    int Q1=floor(float(25*hash_size)/float(100)+1);
    int Q3=floor(float(75*hash_size)/float(100)+1);
    int k=1;
    float Q1_value=0;
    float Q3_value=0;
    for(std::multimap<float,std::string>::iterator it=Ratios.begin(); it!=Ratios.end();++it)
    {
        if(k==Q1)
        {
            Q1_value=it->first;
        }

        if(k==Q3)
        {
            Q3_value=it->first;
        }
        k++;
    }
    //estimate the interquartile range and then compute the cutoff
    float IQR=Q3_value-Q1_value;
    float cutoff=Q1_value-(1.5*IQR);

    float my_sum=0;
    int count=0;
    //output the ratios and estimate the average ratio
    for(std::multimap<float,std::string>::iterator it=Ratios.begin(); it!=Ratios.end();++it)
    {
        if(cutoff>=it->first)
        {
            //it is an outlier so dont consider it
        }
        else
        {
            my_sum=my_sum+it->first;
            count++;
            output<<it->second<<"\t"<<it->first<<std::endl;

        }
    }
    output<<"AVERAGE\t"<<float(my_sum/count)<<std::endl;
    output.close();
}

int checkFileExistance(char *name)
{
    int flag=0;
    if (FILE *file = fopen(name, "r")) 
    {
        fclose(file);
        flag=1;
    } 
    else 
    {
        flag=-1;
    } 
    return flag;
}




