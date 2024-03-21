# Installation

Please install HiCPro(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0831-x) as v3.0.0 and build the supporting environment. Conda is recommended for the process.

# HiC-Pro mapping

Run the following code in terminal

"HiC-Pro_3.0.0/bin/HiC-Pro -i /lustre/home/fhliu/HiCPro_data -o /lustre/home/fhliu/HiCPro_data_results_fat -c /HiCPro_data/configure/config-hicpro.txt"

The configuration file used in the process needs to be carefully edited, here is an example of our work:

##configuration of HiCPro##

#Please change the variable settings below if necessary


##Paths and Settings  - Do not edit !


TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = HiCPro_data


##SYSTEM AND SCHEDULER - Start Editing Here !!

N_CPU = 64
LOGFILE = hicpro.log

JOB_NAME =  hictest
JOB_MEM = 3TB
JOB_WALLTIME = 1000:00:00 
JOB_QUEUE = cpuPartition
JOB_MAIL = 2911040196@qq.com


##Data


PAIR1_EXT = _raw_1
PAIR2_EXT = _raw_2


##Alignment options


FORMAT = phred33
MIN_MAPQ = 0

BOWTIE2_IDX_PATH = /lustre/home/fhliu/HiCPro_data/bowtie2-build/mm9_noYM/mm9
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder


##Annotation files


REFERENCE_GENOME = mm9_noYM
GENOME_SIZE = /lustre/home/fhliu/HiCPro_data/genome_reference/mm9/mm9_noYM.sizes.txt
CAPTURE_TARGET =


##Allele specific analysis


ALLELE_SPECIFIC_SNP = 


##Digestion Hi-C


GENOME_FRAGMENT = /lustre/home/fhliu/HiCPro_data/restriction_sites/mm9_noYM_DpnII.bed
LIGATION_SITE = GATCGATC
MIN_FRAG_SIZE = 
MAX_FRAG_SIZE =
MIN_INSERT_SIZE =
MAX_INSERT_SIZE =


##Hi-C processing


MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 1
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

##Contact Maps


BIN_SIZE = 1000000
MATRIX_FORMAT = upper


##Normalization

MAX_ITER = 10
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1

# results visulization

The output of HiCPro is a matrix file containing interaction pairs. Using juicer tools, it can be converted into hic file and viewed with juicebox. After installation of juicer(https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start), run the following code in the terminal to commit the convertion.

cd HiCPro_data_results_fat/hic_results


find . -name "*allValidPairs" | xargs -i juicertools/hicpro2juicebox.sh -i {} -g genome_reference/mm9/mm9_noYM.sizes.txt -j pkgs/juicer_tools_1.19.02.jar


# call compartments

To call compartments, cool files must firstly be converted from hic files. Using the code bellow to generate cool files with required resolutions. hicexplorer v3.7 should be installed before the process.

"
for sample in *.hic; do
hicConvertFormat --matrices $sample  --outFileName $sample.cool_10k --inputFormat hic --outputFormat cool --resolutions 10000
hicConvertFormat --matrices $sample  --outFileName $sample.cool_25k --inputFormat hic --outputFormat cool --resolutions 25000
hicConvertFormat --matrices $sample  --outFileName $sample.cool_50k --inputFormat hic --outputFormat cool --resolutions 50000
hicConvertFormat --matrices $sample  --outFileName $sample.cool_100k --inputFormat hic --outputFormat cool --resolutions 100000

done

"

After acquiring cool files, we can call compartments with cooltools(v0.4.0).  

"
for sample in *.cool; do
cooltools eigs-cis --n-eigs 3 -o $sample.callcompartment  $sample;
done

"






