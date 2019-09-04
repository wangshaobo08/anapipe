program: anapipe  
version: 0.0.0  
updated: 09:39:19 Sep  4 2019  
Usage: anapipe [OPTIONS]  

|  Options                    | Explanations
|-----------------------------|---------------------------------
|  -h,--help                  | Print this help message and exit
|  -s,--slist FILE REQUIRED   | sample list file
|  -r,--ref FILE REQUIRED     | reference file
|  -b,--bed FILE REQUIRED     | bed region file
|  -t,--gset FILE             | gene list
|  -v,--vread TEXT            | fastq  subset reads number
|  -o,--out TEXT              | output directory
|  -a,--amark INT in [1 - 9]  | analysis marker range
|  -i,--imark INT in [1 - 9]  | initial analysis marker
|  -e,--emark INT in [1 - 9]  | end analysis marker
|  -q,--queue TEXT            | queue to run tasks
|  -c,--ctd                   | continue from last failure
|  -l,--loc                   | run in localhost
|  -g,--gen                   | generate sjms, not run tasks
|  -u,--update Needs: --ctd   | update command to execute
|  -n,--noclean               | not cleanup intermediate files
|  --retryn INT in [0 - 5]=1  | # of retries after failure of one job
Installation

1. clone repo  
`git clone https://github.com/vanNul/anapipe.git`

2. compile  
`cd anapipe`  
`./autogen.sh`  
`./configure --prefix=/path/to/install/dir/`  
`make`  
`make install`  

3. execute  
`/path/to/install/dir/bin/anapipe`  

PS  
this is a sjm based analysis pipeline  

|Marker|Analysis         |Software |Version
|------|-----------------|---------|----------- 
|1     |cutadapter and qc|fqtool   |0.0.0       
|2     |filter ncrna     |filter   |0.0.0       
|3     |downsample fastq |seqtk    |1.3-r106    
|4     |genome slignment |bwa      |0.7.17-r988
|5     |markdup          |duplexer |0.0.0       
|6     |bam QC           |bamqc    |0.0.0       
|7     |express quant    |kallisto |0.45.1 
|8     |report           |anarpt   |0.0.0     
|9     |cleanup          |rm       |8.4         
