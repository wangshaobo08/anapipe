#ifndef GENPIPE_H
#define GENPIPE_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "job.h"
#include "task.h"
#include "pipeline.h"
#include "genjob.h"
#include "options.h"

/** Class to generate pipeline */
class GenPipe{
    public:
        Options* mOpt;   ///< pointer to Options object
        Pipeline* mPipe; ///< pointer to Pipeline object

        /** construct a GenPipe object
         * @param opt pointer to Options object
         * @param p pointer to Pipeline object
         */
        GenPipe(Options* opt, Pipeline* p);
        
        /** destroy a GenPipe object */
        ~GenPipe();

        /** generate library analysisi sub pipeline */
        void genAnalibTask();
};

#endif
