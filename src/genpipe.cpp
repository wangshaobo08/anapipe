#include "genpipe.h"

GenPipe::GenPipe(Options* opt, Pipeline* p){
    mOpt = opt;
    mPipe = p;
}

GenPipe::~GenPipe(){
}

void GenPipe::genAnalibTask(){
    std::ifstream fr1(mOpt->clOpt.sample_list);
    std::string sampleNo, flowCell, libName, read1, read2, tmpStr;
    std::getline(fr1, tmpStr);
    std::istringstream iss;
    int count = 0;
    while(std::getline(fr1, tmpStr)){
        iss.clear();
        iss.str(tmpStr);
        iss >> sampleNo >> flowCell >> libName >> read1 >> read2;
        GenJob* genJob = new GenJob(mOpt);
        Task* task = new Task(7);
        // fqtool
        Job* jFqtool = new Job("fqtool", mOpt->ioOpt.cut_dir, libName, 1);
        genJob->setLib(read1, read2);
        genJob->genFqtoolJob(jFqtool);
        task->addJob(jFqtool, 0);
        // filter
        Job* jFiltdb = new Job("filter", mOpt->ioOpt.fil_dir, libName, 2);
        genJob->setLib(jFqtool->o1, jFqtool->o2);
        genJob->genFilterJob(jFiltdb);
        task->addJob(jFiltdb, 1); 
        // seqtk
        Job* jSeqtk = new Job("seqtk", mOpt->ioOpt.dfq_dir, libName, 3);
        jSeqtk->optPre = sampleNo;
        genJob->setLib(jFiltdb->o1, jFiltdb->o2);
        genJob->genSeqtkJob(jSeqtk);
        task->addJob(jSeqtk, 2);
        // bwa mem
        Job* jAln = new Job("bwa", mOpt->ioOpt.aln_dir, libName, 4);
        genJob->setLib(jSeqtk->o1, jSeqtk->o2);
        genJob->genAlignJob(jAln);
        task->addJob(jAln, 3);
        // mkdup
        Job* jMkdup = new Job("duplexer", mOpt->ioOpt.mkd_dir, libName, 5);
        genJob->setBam(jAln->o1);
        genJob->genMkdupJob(jMkdup);
        task->addJob(jMkdup, 4);
        // bamqc
        Job* jBamqc = new Job("bamqc", mOpt->ioOpt.bqc_dir, libName, 6);
        genJob->setBam(jMkdup->o1);
        genJob->genBamqcJob(jBamqc);
        task->addJob(jBamqc, 5);
        // express
        Job* jExpress = new Job("kallisto", mOpt->ioOpt.exp_dir, libName, 7);
        genJob->setLib(jSeqtk->o1, jSeqtk->o2);
        genJob->genExpressJob(jExpress);
        task->addJob(jExpress, 3);
        // report
        Job* jReport = new Job("anarpt", mOpt->ioOpt.rep_dir, libName, 8);
        genJob->genReportJob(jReport);
        task->addJob(jReport, 6);
        // cleanup
        if(!mOpt->clOpt.noclean){
            Job* jClean = new Job("cleanup", mOpt->ioOpt.log_dir, libName, 9);
            genJob->genCleanupJob(jClean);
            task->addJob(jClean, 7);
        }
        // updating status
        std::string sjmStatus = mOpt->ioOpt.sjm_dir + "/" + libName + "_analib.sjm.status";
        std::map<std::string, std::string> jMap;
        Job::getStatus(jMap, sjmStatus);
        for(auto& e: task->joblist){
            for(auto& f: e){
                if(std::find(mOpt->clOpt.ana_marker.cbegin(), mOpt->clOpt.ana_marker.cend(), f->stage_marker) == mOpt->clOpt.ana_marker.cend()){
                    f->status.second = "done";
                }
                if(mOpt->clOpt.update && jMap.find(f->name.second) != jMap.end()){
                    f->status.second = jMap[f->name.second];
                }
                if(!mOpt->clOpt.queue.empty()){
                    f->queue.second = mOpt->clOpt.queue;
                }
            }
        }
        std::string sjmFile(mOpt->ioOpt.sjm_dir + "/" + libName + "_analib.sjm");
        std::string sjmLog(mOpt->ioOpt.sjm_dir + "/" + libName + "_analib.log");
        RunTask* runTask = new RunTask();
        runTask->sjmCMD = mOpt->ioOpt.bin_dir + "/sjm -i -l " + sjmLog + " " + sjmFile;
        runTask->goodMarkFile = mOpt->ioOpt.log_dir + "/" + libName + ".Analysis.SUCCESS";
        runTask->failMarkFile = mOpt->ioOpt.log_dir + "/" + libName + ".Analysis.FAIL";
        runTask->logFile = sjmLog;
        mPipe->addRunFile(runTask, count, 0);
        std::ofstream fw(sjmFile);
        fw << task;
        fw.close();
    }
    ++count;
}
