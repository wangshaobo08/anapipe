#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <climits>
#include <mutex>
#include "CLI.hpp"
#include <libgen.h>
#include "pipeline.h"
#include "genpipe.h"
#include "options.h"
#include "util.h"

int main(int argc, char** argv)
{
    Options* opt = new Options();
    if(argc == 1){
    std::string sys_cmd = std::string(argv[0]) + " -h";
        std::system(sys_cmd.c_str());
        return 0;
    }
    if(argc == 2 && strcmp(argv[1], "-L") == 0){
        opt->showMark();
        return 0;
    }
    std::mutex logmtx;
    std::string cmp_time = std::string(__TIME__) + " " + std::string(__DATE__);
    CLI::App app("program: " + std::string(basename(argv[0])) + "\nversion: " + opt->version + "\nupdated: " + cmp_time);
    app.get_formatter()->column_width(36);
    app.add_option("-s,--slist", opt->clOpt.sample_list, "sample list file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-r,--ref", opt->clOpt.ref, "reference file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-b,--bed", opt->clOpt.reg, "bed region file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-t,--gset", opt->clOpt.gset, "gene list")->check(CLI::ExistingFile);
    app.add_option("-v,--vread", opt->clOpt.dfq_vol, "fastq subset reads number");
    app.add_option("-o,--out", opt->ioOpt.out_dir, "output directory");
    app.add_option("-a,--amark", opt->clOpt.ana_marker, "analysis marker range")->check(CLI::Range(opt->clOpt.minstage, opt->clOpt.maxstage));
    app.add_option("-i,--imark", opt->clOpt.ini_marker, "initial analysis marker")->check(CLI::Range(opt->clOpt.minstage, opt->clOpt.maxstage));
    app.add_option("-e,--emark", opt->clOpt.end_marker, "end analysis marker")->check(CLI::Range(opt->clOpt.minstage, opt->clOpt.maxstage));
    app.add_option("-q,--queue", opt->clOpt.queue, "queue to run tasks");
    CLI::Option* prerun = app.add_flag("-c,--ctd", opt->clOpt.rerun, "continue from last failure");
    app.add_flag("-l,--loc", opt->clOpt.local, "run in localhost");
    app.add_flag("-g,--gen", opt->clOpt.gensjm, "generate sjms, not run tasks");
    app.add_flag("-u,--update", opt->clOpt.update, "update command to execute")->needs(prerun);
    app.add_flag("-n,--noclean", opt->clOpt.noclean, "not cleanup intermediate files");
    CLI_PARSE(app, argc, argv);
    util::loginfo("parsing arguments finished.", logmtx);
    std::string cmd = "";
    for(int i = 0; i < argc; ++i){
        cmd.append(argv[i] + std::string(1, ' '));
    }
    util::loginfo("command: " + cmd, logmtx);
    util::loginfo("preparing output parent directory.", logmtx);
    util::makedir(opt->ioOpt.out_dir);
    util::loginfo("output parent directory prepared.", logmtx);
    util::loginfo("updating arguments.", logmtx);
    opt->updateOptions();
    util::loginfo("arguments updated.", logmtx);
    Pipeline* p = new Pipeline(opt->nSubPipe, opt->failMarkFile, opt->goodMarkFile);
    util::loginfo("pipeline object construced.", logmtx);
    util::loginfo("generate subdirectories.", logmtx);
    opt->genDirectory();
    util::loginfo("subdirectories generated.", logmtx);
    GenPipe* g = new GenPipe(opt, p);
    util::loginfo("generate library analysis pipeline.", logmtx);
    g->genAnalibTask();
    util::loginfo("library analysis pipeline generated.", logmtx);
    if(opt->clOpt.rerun){
        util::loginfo("resume last running of pipeline started.", logmtx);
        p->prepareRerun();
        util::loginfo("finished resume running.", logmtx);
    }
    if(!opt->clOpt.gensjm){
        util::loginfo("running pipeline now.", logmtx);
        p->runPipeline();
        util::loginfo("pipeline finished.", logmtx);
    }
}
