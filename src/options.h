#ifndef OPTIONS_H
#define OPTIONS_H

#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include "json.hpp"
#include "util.h"

/** struct to store input/output diecttory Options */
class IODirectoryOptions{
    public:
        std::string out_dir; ///< output directory of analysis
        std::string bin_dir; ///< binary directory of pipeline
        std::string db_dir ; ///< database directory of pipeline
        std::string sjm_dir; ///< directory to store sjm control files
        std::string cut_dir; ///< directory to store library splitted
        std::string fil_dir; ///< directory to store library filtered by databases such as rRNA
        std::string dfq_dir; ///< directory to store library downsampled to predefined volume
        std::string aln_dir; ///< directory to store alignment results of library
        std::string mkd_dir; ///< directory to store markdup results of library
        std::string bqc_dir; ///< directory to store alignment qc results of library
        std::string exp_dir; ///< directory to store expression quantification of library
        std::string rep_dir; ///< directory to store final analysis report of pipeline
        std::string log_dir; ///< directory to store analysis success/fail status of pipeline
    
    IODirectoryOptions(){
        out_dir = "./";           
        bin_dir = "./";           
        db_dir = "./";            
        sjm_dir = "00.sjm";       
        cut_dir = "01.fqtool";
        fil_dir = "02.filter";    
        dfq_dir = "03.seqtk";     
        aln_dir = "04.align"; 
        mkd_dir = "05.mkdup";   
        bqc_dir = "06.bamqc";     
        exp_dir = "07.express";   
        rep_dir = "report";       
        log_dir = "log";
    }    
};

/** class to store pipeline execution control options */
class PipeControlOptions{
    public:
        int minstage;                ///< minimum stages this pipeline supports
        int maxstage;                ///< maximum stages this pipeline supports
        std::string sample_list;     ///< sample list path
        std::string ref;             ///< reference fasta file path
        std::string reg;             ///< region bed file path
        std::string gset;            ///< gene set tsv file path
        std::string queue;           ///< sge queue name
        std::vector<int> ana_marker; ///< analysis marker number
        int ini_marker;              ///< starting analysis stage marker
        int end_marker;              ///< ending analysis stage marker
        std::string dfq_vol;         ///< fastq reads number to downsample to
        bool rerun;                  ///< resume running from last failing if true
        bool local;                  ///< running pipeline in local server if true
        bool gensjm;                 ///< only generate sjm control files without running pipeline if true
        bool update;                 ///< update sjm control file before running pipeline if true
        bool noclean;                ///< do not cleanup intermediate files after running pipeline if true

    PipeControlOptions(){
        minstage = 1;
        maxstage = 9;
        sample_list = "";
        ref = "";
        reg = "";
        gset = "";
        queue = "";
        ana_marker = {};
        ini_marker = 1;
        end_marker = 9;
        dfq_vol = "VOL";
        rerun = false;
        local = false;
        gensjm = false;
        update = false;
        noclean = false;
    }
};

/** struct to store various  options of pipeline */
class Options{
    public:
        IODirectoryOptions ioOpt;                  ///< IODirectoryOptions object
        PipeControlOptions clOpt;                  ///< PipeControlOptions object
        std::string version = "0.0.0";             ///< pipeline version
        size_t minFqVol = 0;                       ///< minimum fastq volume in sample list
        int nSubPipe = 1;                          ///< total sub pipeline needed
        int nSamples = 0;                          ///< total samples to execute
        std::string goodMarkFile;                  ///< pipeline success markfile
        std::string failMarkFile;                  ///< pipeline failure markfile
    public:
        /** construct a Options object */
        Options();
        /** destroy a Options object */
        ~Options();

        /** update Options after commandline arguments parse */
        void updateOptions();

        /** update minFqVolMap for down sample fastq to samve volume in each pooled library after filter db */
        void updateMinFqVolMap();

        /** generate subdirectories */
        void genDirectory();

        /** show marker help info of each stage */
        static void showMark();
};

#endif
