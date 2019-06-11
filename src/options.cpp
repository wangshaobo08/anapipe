#include "util.h"
#include "dirutil.h"
#include "options.h"

Options::Options(){
}

Options::~Options(){
}

void Options::updateOptions(){
    clOpt.sample_list = util::abspath(clOpt.sample_list);
    if(clOpt.ana_marker.empty()){
        for(int i = clOpt.ini_marker; i <= clOpt.end_marker; ++i){
            clOpt.ana_marker.push_back(i);
        }
    }else{
        std::remove_if(clOpt.ana_marker.begin(), clOpt.ana_marker.end(), 
                       [&](int& e){return e < clOpt.ini_marker || e > clOpt.end_marker;}
                       );
    }
    if(clOpt.gset.empty()){
        clOpt.gset = ioOpt.db_dir + "gset/kgset.tsv";
    }

    ioOpt.bin_dir = util::dirname(dirutil::getExecutablePath());
    ioOpt.out_dir = util::abspath(ioOpt.out_dir);
    ioOpt.db_dir = util::dirname(ioOpt.bin_dir) + "/db/";
    
    std::string sep = "/";
    ioOpt.sjm_dir = ioOpt.out_dir + sep + ioOpt.sjm_dir;
    ioOpt.cut_dir = ioOpt.out_dir + sep + ioOpt.cut_dir;
    ioOpt.fil_dir = ioOpt.out_dir + sep + ioOpt.fil_dir;
    ioOpt.dfq_dir = ioOpt.out_dir + sep + ioOpt.dfq_dir;
    ioOpt.aln_dir = ioOpt.out_dir + sep + ioOpt.aln_dir;
    ioOpt.mkd_dir = ioOpt.out_dir + sep + ioOpt.mkd_dir;
    ioOpt.bqc_dir = ioOpt.out_dir + sep + ioOpt.bqc_dir;
    ioOpt.exp_dir = ioOpt.out_dir + sep + ioOpt.exp_dir;
    ioOpt.rep_dir = ioOpt.out_dir + sep + ioOpt.rep_dir;
    ioOpt.log_dir = ioOpt.out_dir + sep + ioOpt.log_dir;

    std::ifstream fr(clOpt.sample_list);
    std::string line;
    int count = 0;
    while(std::getline(fr, line)){
        ++count;
    }
    nSamples = count;
    goodMarkFile = ioOpt.log_dir + "/SUCCESS";
    failMarkFile = ioOpt.log_dir + "/FAIL";
    if(clOpt.dfq_vol == "0"){
        updateMinFqVolMap();
    }
}

void Options::updateMinFqVolMap(){
    if(!util::exists(clOpt.sample_list)){
        util::error_exit("configure file \"" + clOpt.sample_list + "\" does not existd!\n");
    }
    std::vector<size_t> readsNum;
    std::ifstream fr(clOpt.sample_list);
    std::string line, logfile;
    std::vector<std::string> vstr;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        logfile = util::joinpath(ioOpt.fil_dir, vstr[2] + ".filter.json");
        std::ifstream fr1(logfile);
        jsn::json j;
        fr1 >> j;
        readsNum.push_back(j["FilterResult"]["Summary"]["ReadsGot"]);
        fr1.close();
    }
    minFqVol = *std::min_element(readsNum.begin(), readsNum.end());
}

void Options::genDirectory(){
    util::makedir(ioOpt.sjm_dir);
    util::makedir(ioOpt.cut_dir);
    util::makedir(ioOpt.fil_dir);
    util::makedir(ioOpt.dfq_dir);
    util::makedir(ioOpt.aln_dir);
    util::makedir(ioOpt.mkd_dir);
    util::makedir(ioOpt.bqc_dir);
    util::makedir(ioOpt.exp_dir);
    util::makedir(ioOpt.rep_dir);
    util::makedir(ioOpt.log_dir);
}

void Options::showMark(){
    std::vector<std::pair<std::string, std::string>> stg;
    stg.push_back({"cutadapter and qc", "fqtool"});
    stg.push_back({"filter ncrna",      "filter"});
    stg.push_back({"downsample fastq",  "seqtk"});
    stg.push_back({"genome slignment",  "bwa"});
    stg.push_back({"markdup",           "duplexer"});
    stg.push_back({"bam QC",            "bamqc"});
    stg.push_back({"express quant",     "kallisto"});
    stg.push_back({"report",            "anarpt"});
    stg.push_back({"cleanup",           "rm"});
    
    std::map<std::string, std::string> sver;
    sver["fqtool"]     = "0.0.0";
    sver["filter"]     = "0.0.0";
    sver["seqtk"]      = "1.3-r106";
    sver["bwa"]        = "0.7.17-r1188";
    sver["duplexer"]   = "0.0.0";
    sver["bamqc"]      = "0.0.0";
    sver["kallisto"]   = "0.45.1";
    sver["anarpt"]     = "0.0.0";
    sver["rm"]         = "8.4";
    
    std::cout << std::left;
    std::cout << "  ┌------┬-----------------┬---------┬------------┐" << std::endl;
    std::cout << "  |" << std::setw(6) << "Marker" << "|" << std::setw(17) << "Analysis" << "|";
    std::cout << std::setw(9) << "Software" << "|" << std::setw(12) << "Version" << "|" << std::endl;
    for(size_t i = 0; i < stg.size(); ++i){
    std::cout << "  |------┼-----------------┼---------┼------------|" << std::endl;
        std::cout << "  |" << std::setw(6) << i + 1 << "|" << std::setw(17) << stg[i].first << "|";
        std::cout << std::setw(9) << stg[i].second << "|" << std::setw(12) << sver[stg[i].second] << "|" << std::endl;
    }
    std::cout << "  └------┴-----------------┴---------┴------------┘" << std::endl;
}
