#ifndef GENJOB_H
#define GENJOB_H

#include <string>
#include "job.h"
#include "options.h"

/** Class to generate a job */
class GenJob{
    public:
        Options* mOpt;      ///< pointer to Options object
        std::string bam;    ///< bam file
        std::string lib1;   ///< fastq read1 file path
        std::string lib2;   ///< fastq read2 file path
    public:
        /** Construct a GenJob object
         * @param opt pointer to Options object
         */
        GenJob(Options* opt);
        
        /** Destroy a GenJob object */
        ~GenJob();

        /** set library of GenJob 
         * @param l1 read1 file path
         * @param l2 read2 file path
         */
        void setLib(const std::string& l1, const std::string& l2);
        
        /** set bam of GenJob
         * @param b bam file path
         */
        void setBam(const std::string& b);

        /** generate spliter Job
         * @param conf barcode configure file path
         * @param j pointer to Job
         */
        void genSpliterJob(const std::string& conf, Job* j);
        
        /** generate fqtool Job
         * @param j pointer to Job
         */
        void genFqtoolJob(Job* j);
        
        /** generate seqtk Job
         * @param j pointer to Job
         */
        void genSeqtkJob(Job* j);
        
        /** generate filter Job
         * @param j pointer to Job
         */
        void genFilterJob(Job* j);
        
        /** generate alignment Job
         * @param j pointer to Job
         */ 
        void genAlignJob(Job* j);
        
        /** generate markdup Job
         * @param j pointer to Job
         */
        void genMkdupJob(Job* j);
        
        /** generate bamqc Job
         * @param j pointer to Job
         */
        void genBamqcJob(Job* j);
        
        /** generate kallisto Job
         * @param j pointer to Job
         */
        void genExpressJob(Job* j);
        
        /** generate cleanup Job
         * @param j pointer to Job
         */
        void genCleanupJob(Job* j);
};

#endif
