#include "genjob.h"
#include "util.h"

GenJob::GenJob(Options* opt){
    mOpt = opt;
}

GenJob::~GenJob(){
}

void GenJob::setLib(const std::string& l1, const std::string& l2){
    lib1 = l1;
    lib2 = l2;
}

void GenJob::setBam(const std::string& b){
    bam = b;
}

void GenJob::genFqtoolJob(Job* j){
    std::string ofq1 = j->workdir.second + j->pre + ".R1.fq";
    std::string ofq2 = j->workdir.second + j->pre + ".R2.fq";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/fqtool";
    j->cmd.second += " -q -a -x";
    j->cmd.second += " -i " + lib1 + " -I " + lib2;
    j->cmd.second += " -o " + ofq1;
    j->cmd.second += " -O " + ofq2;
    j->cmd.second += " -J " + j->workdir.second + j->pre + ".fqtool.json";
    j->cmd.second += " -H " + j->workdir.second + j->pre + ".fqtool.html";
    j->memory.second = "1g";
    j->slots.second = "4";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = ofq1;
    j->o2 = ofq2;
}

void GenJob::genSeqtkJob(Job* j){
    j->cmd.second = "rm -f " + j->workdir.second + util::basename(lib1);
    j->cmd.second += " && rm -f " + j->workdir.second + util::basename(lib2);
    j->cmd.second += " && ln -sf " + lib1 + " " + j->workdir.second;
    j->cmd.second += " && ln -sf " + lib2 + " " + j->workdir.second;
    if(mOpt->clOpt.dfq_vol != "VOL"){
        size_t vol = std::atoi(mOpt->clOpt.dfq_vol.c_str());
        if(vol == 0){
            vol = mOpt->minFqVol;
        }
        if(vol != 0){
            j->cmd.second = "rm -f " + j->workdir.second + util::basename(lib1);
            j->cmd.second += " && rm -f " + j->workdir.second + util::basename(lib2);
            j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/seqtk sample -2 -s 100 " + lib1 + " " + std::to_string(vol);
            j->cmd.second += " > " + j->workdir.second + util::basename(lib1);
            j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/seqtk sample -2 -s 100 " + lib2 + " " + std::to_string(vol);
            j->cmd.second += " > " + j->workdir.second + util::basename(lib2);
        }
    }
    j->memory.second = "1g";
    j->slots.second = "1";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = j->workdir.second + util::basename(lib1);
    j->o2 = j->workdir.second + util::basename(lib2);
}

void GenJob::genFilterJob(Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/filter -d";
    j->cmd.second += " -i " + lib1;
    j->cmd.second += " -I " + lib2;
    j->cmd.second += " -r " + mOpt->ioOpt.db_dir + "/ncrna/Homo_sapiens.GRCh37.ncrna.class.fa";
    j->cmd.second += " -o " + j->workdir.second + "/" + util::basename(lib1);
    j->cmd.second += " -O " + j->workdir.second + "/" + util::basename(lib2);
    j->cmd.second += " -l " + j->workdir.second + "/" + j->pre + ".filter.json";
    j->memory.second = "1g";
    j->slots.second = "6";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
    j->o1 = j->workdir.second + "/" + util::basename(lib1);
    j->o2 = j->workdir.second + "/" + util::basename(lib2);
}

void GenJob::genAlignJob(Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/bwa mem";
    j->cmd.second += " -t 8 " + mOpt->clOpt.ref + " " + lib1 + " " + lib2;
    j->cmd.second += " | " + mOpt->ioOpt.bin_dir + "/samtools sort -@ 8";
    j->cmd.second += " -o " + j->workdir.second + j->pre + ".aln.sort.bam";
    j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/samtools index ";
    j->cmd.second += j->workdir.second + j->pre + ".aln.sort.bam";
    j->memory.second = "8g";
    j->slots.second = "8";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = j->workdir.second + j->pre + ".aln.sort.bam";
}

void GenJob::genMkdupJob(Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/duplexer markdup";
    j->cmd.second += " -i " + bam;
    j->cmd.second += " -l " + j->workdir.second + "/" + j->pre + ".mkdup.json";
    j->cmd.second += " | " + mOpt->ioOpt.bin_dir + "/samtools sort -@ 8 -o " + j->workdir.second + j->pre + ".mkdup.sort.bam";
    j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/samtools index ";
    j->cmd.second += j->workdir.second + j->pre + ".mkdup.sort.bam";
    j->memory.second = "3g";
    j->slots.second = "3";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = j->workdir.second + j->pre + ".mkdup.sort.bam";
}

void GenJob::genBamqcJob(Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/bamqc";
    j->cmd.second += " -i " + bam;
    j->cmd.second += " -b " + mOpt->clOpt.reg;
    j->cmd.second += " -r " + mOpt->clOpt.ref;
    j->cmd.second += " -o " + j->workdir.second + "/" + j->pre + ".bamqc.json";
    j->cmd.second += " -I --reg " + mOpt->clOpt.reg;
    j->cmd.second += " -U -l " + mOpt->ioOpt.db_dir + "/refMrna/refGene.3utr.len";
    j->memory.second = "1g";
    j->slots.second = "4";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
}

void GenJob::genExpressJob(Job* j){
    j->cmd.second += "mkdir -p " + j->workdir.second + j->pre + " && ";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/kallisto";
    j->cmd.second += " quant -t 8 -i " + mOpt->ioOpt.db_dir + "/ensembl/Homo_sapiens.GRCh37.cdna.all.fa.idx";
    j->cmd.second += " -o " + j->workdir.second + j->pre + " ";
    j->cmd.second += lib1 + " " + lib2;
    j->memory.second = "4g";
    j->slots.second = "8";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
}

void GenJob::genCleanupJob(Job* j){
    j->cmd.second = "rm -f ";
    j->cmd.second += mOpt->ioOpt.cut_dir + "/" + j->pre + "*.fq ";
    if(mOpt->clOpt.dfq_vol != "VOL"){
        j->cmd.second += mOpt->ioOpt.fil_dir + "/" + j->pre + "*.fq ";
    }
    j->cmd.second += mOpt->ioOpt.aln_dir + "/" + j->pre + "*.aln.bam ";
    j->cmd.second += mOpt->ioOpt.mkd_dir + "/" + j->pre + "*.mkdup.bam ";
    j->memory.second = "1g";
    j->slots.second = "1";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
}
