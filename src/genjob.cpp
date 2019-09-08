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
    std::string ofq1 = j->workdir.second + j->pre + ".R1.fq.gz";
    std::string ofq2 = j->workdir.second + j->pre + ".R2.fq.gz";
    std::string jrpt = j->workdir.second + j->pre + ".fqtool.json";
    std::string hrpt = j->workdir.second + j->pre + ".fqtool.html"; 
    j->cmd.second += mOpt->ioOpt.bin_dir + "/fqtool";
    j->cmd.second += " -q -a -x -d";
    j->cmd.second += " -i " + lib1 + " -I " + lib2;
    j->cmd.second += " -o " + ofq1;
    j->cmd.second += " -O " + ofq2;
    j->cmd.second += " -J " + jrpt;
    j->cmd.second += " -H " + hrpt;
    j->memory.second = "1g";
    j->slots.second = "4";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = ofq1;
    j->o2 = ofq2;
}

void GenJob::genSeqtkJob(Job* j){
    std::string bnlib1 = util::basename(lib1);
    std::string bnlib2 = util::basename(lib2);
    std::string ofq1 = j->workdir.second + bnlib1;
    std::string ofq2 = j->workdir.second + bnlib2;
    j->cmd.second = "rm -f " + ofq1;
    j->cmd.second += " && rm -f " + ofq2;
    j->cmd.second += " && ln -sf " + lib1 + " " + j->workdir.second;
    j->cmd.second += " && ln -sf " + lib2 + " " + j->workdir.second;
    if(mOpt->clOpt.dfq_vol != "VOL"){
        size_t vol = std::atoi(mOpt->clOpt.dfq_vol.c_str());
        if(vol == 0) vol = mOpt->minFqVol;
        if(vol != 0){
            j->cmd.second = "rm -f " + ofq1;
            j->cmd.second += " && rm -f " + ofq2;
            j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/seqtk sample -2 -s 100";
            j->cmd.second += " -o " + ofq1 + " " + lib1 + " " + std::to_string(vol);
            j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/seqtk sample -2 -s 100";
            j->cmd.second += " -o " + ofq2 + " " + lib2 + " " + std::to_string(vol);
        }
    }
    j->memory.second = "1g";
    j->slots.second = "1";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = ofq1;
    j->o2 = ofq2;
}

void GenJob::genFilterJob(Job* j){
    std::string ncrf = mOpt->ioOpt.db_dir + "/ncrna/Homo_sapiens.GRCh37.ncrna.class.fa";
    std::string ofq1 = j->workdir.second + "/" + util::basename(lib1);
    std::string ofq2 = j->workdir.second + "/" + util::basename(lib2);
    std::string jrpt = j->workdir.second + "/" + j->pre + ".filter.json";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/filter -d";
    j->cmd.second += " -i " + lib1;
    j->cmd.second += " -I " + lib2;
    j->cmd.second += " -r " + ncrf;
    j->cmd.second += " -o " + ofq1;
    j->cmd.second += " -O " + ofq2;
    j->cmd.second += " -l " + jrpt;
    j->memory.second = "1g";
    j->slots.second = "6";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
    j->o1 = ofq1;
    j->o2 = ofq2;
}

void GenJob::genAlignJob(Job* j){
    std::string obam = j->workdir.second + j->pre + ".aln.sort.bam";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/bwa mem";
    j->cmd.second += " -t 8 " + mOpt->clOpt.ref + " " + lib1 + " " + lib2;
    j->cmd.second += " | " + mOpt->ioOpt.bin_dir + "/samtools sort -@ 8 -o " + obam;
    j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/samtools index " + obam;
    j->memory.second = "8g";
    j->slots.second = "8";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = obam;
}

void GenJob::genMkdupJob(Job* j){
    std::string jrpt = j->workdir.second + j->pre + ".mkdup.json";
    std::string obam = j->workdir.second + j->pre + ".mkdup.sort.bam";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/duplexer markdup";
    j->cmd.second += " -i " + bam;
    j->cmd.second += " -l " + jrpt;
    j->cmd.second += " | " + mOpt->ioOpt.bin_dir + "/samtools sort -@ 8 -o " + obam;
    j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/samtools index " + obam;
    j->memory.second = "3g";
    j->slots.second = "3";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = obam;
}

void GenJob::genBamqcJob(Job* j){
    std::string jrpt = j->workdir.second + j->pre + ".bamqc.json";
    std::string creg = mOpt->ioOpt.db_dir + "/regfile/refMrna.CDS.bed";
    std::string utr3 = mOpt->ioOpt.db_dir + "/refMrna/refGene.3utr.len";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/bamqc";
    j->cmd.second += " -i " + bam;
    j->cmd.second += " -b " + mOpt->clOpt.reg;
    j->cmd.second += " -r " + mOpt->clOpt.ref;
    j->cmd.second += " -o " + jrpt;
    j->cmd.second += " -I --reg " + creg;
    j->cmd.second += " -U -l " + utr3;
    j->memory.second = "1g";
    j->slots.second = "4";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
}

void GenJob::genExpressJob(Job* j){
    std::string cidx = mOpt->ioOpt.db_dir + "/ensembl/Homo_sapiens.GRCh37.cdna.all.fa.idx";
    j->cmd.second += "mkdir -p " + j->workdir.second + j->pre + " && ";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/kallisto";
    j->cmd.second += " quant -t 8 -i " + cidx;
    j->cmd.second += " -o " + j->workdir.second + j->pre + " ";
    j->cmd.second += lib1 + " " + lib2;
    j->memory.second = "4g";
    j->slots.second = "8";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
}

void GenJob::genReportJob(Job* j){
    std::string filtlog = mOpt->ioOpt.fil_dir + "/" + j->pre + ".filter.json";
    std::string bamqc = mOpt->ioOpt.bqc_dir + "/" + j->pre + ".bamqc.json";
    std::string abundance = mOpt->ioOpt.exp_dir + "/" + j->pre + "/abundance.tsv";
    std::string ens2gen = mOpt->ioOpt.db_dir + "/NCBI/ensebml2genename";
    std::string outf = mOpt->ioOpt.rep_dir + "/" + j->pre + ".report.xlsx";
    j->cmd.second = mOpt->ioOpt.bin_dir + "/anarpt";
    j->cmd.second += " -f " + filtlog + " -b " + bamqc + " -k " + abundance + " -e " + ens2gen + " -o " + outf;
    j->memory.second = "1g";
    j->slots.second = "1";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
}

void GenJob::genCleanupJob(Job* j){
    std::string cutfq = mOpt->ioOpt.cut_dir + "/" + j->pre + "*.fq.gz";
    std::string filfq = mOpt->ioOpt.fil_dir + "/" + j->pre + "*.fq.gz";
    std::string alnbam = mOpt->ioOpt.aln_dir + "/" + j->pre + "*.aln.bam";
    std::string mkdbam = mOpt->ioOpt.mkd_dir + "/" + j->pre + "*.mkdup.bam";
    j->cmd.second = "rm -f ";
    j->cmd.second += cutfq + " " + alnbam + " " + mkdbam;
    if(mOpt->clOpt.dfq_vol != "VOL") j->cmd.second += " " + filfq;
    j->memory.second = "1g";
    j->slots.second = "1";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
}
