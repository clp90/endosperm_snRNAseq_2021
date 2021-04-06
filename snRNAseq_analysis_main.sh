#!/bin/bash

# ------------------------------------------------------------------------------------
# v1.0 by Colette L. Picard
# 02/20/2021
# ------------------------------------------------------------------------------------

# -------------------------
# Version history:
# v.1.0: final 1.0 - 02/20/2021
# -------------------------

# Description printed when "help" option specified:
read -d '' usage <<"EOF"
v.1.0 by Colette L Picard, 02/20/2021

Usage:
(1) To start from beginning: 
	./snRNAseq_analysis_main.sh -i indir -o outdir -d required_files/snRNAseq_samplist.txt

(2) To start from a specific step, for example step 3:
	./snRNAseq_analysis_main.sh -o outdir -d required_files/snRNAseq_samplist.txt -3

(3) To test that all dependencies are properly installed/can be located without running the script itself:
	./snRNAseq_analysis_main.sh -0

This is the master script for project looking at expression in single nuclei sorted from
A. thaliana endosperm, published 2021 (Picard et al.). This script is designed to 
reproduce every analysis in the paper. These analyses were run on the Whitehead
Institute cluster, which is an LSF cluster. This script farms out various jobs
to the cluster using -bsub-, and will need to be adjusted to run on clusters using
different schedulers (e.g. -qsub-). Dependencies are listed below.

Analysis is divided into these main steps:

******************************************************
Step 1) Align reads to metagenome, get alignment stats (Fig. S1)
Step 2) Get per-gene read counts and CPMs for each library (Fig. S1)
Step 3) Identify nuclei passing basic QC, filter out neg/pos ctrls and examine spike-ins (Fig. S1, S2)
Step 4) Clustering by tSNE and SC3, assign tissue of origin based on clustering (Fig. 1, Ext. data Figs. 1, 2, 3)
Step 5) Cell cycle analysis (Ext. data Fig. 7, Figs. S6, S7)
Step 6) Identify DE genes between 4 DAP endosperm and seed coat clusters
Step 7) Identify imprinted genes from 4 DAP endosperm nuclei (Fig. 3, Ext. data Fig. 8, Figs. S8-S12)
Step 8) Characterize clusters by examining patterns of overall and allele-specific expression over cell cycle & clusters (Figs. 1, 3, Ext. data Figs. 3-5, 9, Figs. S4,S5,S13-S16)
Step 9) Additional analyses/plots (Fig. 3, Ext. data Fig. 3, 9, Figs. S2, S3, S13, S16)
******************************************************

Required arguments:

(1) Required file snRNAseq_samplist.txt is a tab-del file containing a description of each library, currently with the following
	fields (in this order!):
	sampID	project	sort_date	well	libname	seqtype	cross	dap	n_nuc	peak	stage	stagecode	batch	platform	final ERCC dil.
	D17-187001	171003Geh	8/10/17	2A	P1_2A	single	CxV	3	1	6N	early heart	2	1	HiSeq 2000	0
	D17-187002	171003Geh	8/10/17	2B	P1_2B	single	CxV	3	1	6N	early heart	2	1	HiSeq 2000	1:10M
	D17-187003	171003Geh	8/10/17	2C	P1_2C	single	CxV	3	1	6N	early heart	2	1	HiSeq 2000	1:50M
	D17-187004	171003Geh	8/10/17	2D	P1_2D	single	CxV	3	1	6N	early heart	2	1	HiSeq 2000	1:100M
	D17-187005	171003Geh	8/10/17	2E	P1_2E	single	CxV	3	1	6N	early heart	2	1	HiSeq 2000	1:200M

	This file is provided in the github repo here: required_files/snRNAseq_samplist.txt

(2) Point -i to the directory containing all the raw sequencing files named according to the 'libname' field in required_files/snRNAseq_samplist.txt. 
	For example, the fastq file for the first sample in the file (P1_2A) should be named P1_2A.fastq. Paired-end files (run on 
	Nextseq) should be named e.g. P11_1A_R1.fastq and P11_1A_R2.fastq for forward and reverse reads, respectively.

(3) Point -o to a suitable output directory. Can use an existing directory, although
	files will be re-generated and overwritten if they existed already in that directory so use with
	caution.

******************************************************

Notes: 

(1) Although this script should rerun every analysis in the paper, it is very long, has a lot of dependencies and is designed 
	to run on an LSF cluster and uses bsub a lot to speed up analyses by running things in parallel, so it probably won't work 
	out of the box on most setups. To facilitate running, it's separated into 9 main parts and can be restarted from the beginning 
	of each part without re-running the previous parts. For example,

	./snRNAseq_analysis_main.sh -o myoutdir -d required_files/snRNAseq_samplist.txt -3

	will restart the analysis from step 3, assuming steps 1 and 2 already successfully completed and are stored
	in myoutdir/ already. Dependencies are provided in the github where possible, and otherwise are noted when possible.

(2) Used STAR-indexed metagenomes were created of TAIR10+Ler+ERCC and TAIR10+Cvi+ERCC and are included in the Github.
	Mapping and metagenome generation was carried out as described in this paper and associated github:
	https://dx.doi.org/10.1007/978-1-0716-0179-2_13
	https://github.com/clp90/imprinting_analysis

(3) Relevant files from that github are also included in this repo for simplicity.
	SNP files, annotation GTF files etc. are also included in github in required_files/
	
(4) All custom helper scripts, etc. are provided in either the same directory as this script, or
	the scripts/ subdirectory. By default, this script will look for those files in those locations,
	and should all be visible without the user needing to specify anything. However, if any of those
	files are moved relative to this script (e.g. the scripts/ subdir is moved somewhere else, without
	moving this script too), then you can re-point to the required scripts and files using the -p, -s
	and -f flags.

******************************************************

Usage:
snRNAseq_analysis_main.sh [options] -d descfile.txt -o outdir

User-specified options:
Required arguments:
	-d descfile : file containing list of samples (see above for expected fields)
	-o outdir : where to put all output files from this project
Additional arguments:
	-i indir : folder containing all the raw data from all libraries (see above for expected format); only required if starting from step 1
	-p path_to_pipelines : path to folder containing all required pipeline scripts (see list below) [default assumes in same directory as this script]
	-s path_to_scripts : path to folder containing all required helper scripts (see list below) [default assumes in a scripts/ subdirectory relative to this script]
	-f path_to_files : path to folder containing all required other files (see list below) [default assumes in a required_files/ subdirectory relative to this script]
Flag options:
	-0 : checks that all required programs installed on PATH and all required helper scripts and other files can be located, then exits without running
	-[step] : skip to step [step] ([step] can be values 2-9, for example -8 will skip to step 8)
	-h : prints this version and usage information

******************************************************

Dependencies:
(1) Required installed on user PATH:
	- `bedtools` (by Aaron R Quinlan and Ira M Hall) [v2.23.0] [link](https://bedtools.readthedocs.io/en/latest/)
	- `samtools` ((c) Genome Research Ltd.) [v.1.11] [link](http://www.htslib.org)
	- `python` v.2.x.x (Python Software Foundation) [v2.7.17] [link](https://www.python.org/downloads/release/python-2717/)
	- `R` (R Core Team) [v3.6.3] [link](https://www.r-project.org)
  - `fastqc` (Andrews 2010) [v0.11.8] [link](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
  - `trim_galore` (Krueger 2012) [v0.5.0] [link](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
  - `samtools` (Li *et al.* 2009) [v1.9] [link](http://www.htslib.org/)
  - `STAR` (Dobin *et al.* 2012) [v2.6.1d] [link](https://github.com/alexdobin/STAR)
  - `Java 1.8+` [v1.8.0_191] [link](https://www.java.com/en/download/manual.jsp)
	
(2) Must be in path_to_pipelines - set path to this folder using -p option:
	- rna_seq_map.sh (by CLP)
	
(3) Must be in path_to_scripts - set path to this folder using -s option:
	- MarkDuplicates.jar (part of the picard-tools suite, Broad Institute)
	- assign_to_allele.py (by CLP)
	- barchart.R (by CLP)
	- cluster_gene_expression.R (by CLP)
	- gene_overlaps_dotplot.R (by CLP)
	- merge_by_column.R (by CLP)
	- merge_many_files.sh (by CLP)
	- piechart.R (by CLP)
	- plot_SC3_heatmap.py (by CLP)
	- run_DEsingle.R (by CLP)
	- run_topGO.R (by CLP)
	- scatterplot.R (by CLP)
	- single_cell_ASE_analysis.R (by CLP)
	- single_cell_ASE_src.R (by CLP)
	- single_cell_RNAseq_plots.R (by CLP)
	- single_cell_cluster_PCA_tSNE.R (by CLP)
	- single_cell_cluster_SC3.R (by CLP)
	- single_cell_simulate_reads.R (by CLP)
	- single_cell_trajectory_analysis.R (by CLP)
	- subset_large_file.R (by CLP)
	- check_dependencies.R (by CLP)

(4) Must be in path_to_files - set path to this folder using -f option:
	- TAIR10_plus_araport11_nonoverlapping.gtf
	- araport11_TEs_not_in_genes.gtf
	- ERCC_spikeins/ERCC92.gtf
	- genomes/TAIR10_Cvi_ERCC_OH50
	- genomes/TAIR10_Cvi_metachrom.txt
	- Col_Cvi_SNPs_clean.bed
	- genomes/TAIR10_Ler_ERCC_OH50
	- genomes/TAIR10_Ler_metachrom.txt
	- Col_Ler_SNPs_clean.bed
	- ERCC_chrlist.txt 
	
(5) Required R packages:
	- See scripts/check_dependencies.R	
	
(6) Required python packages:
	- argparse (Thomas Waldmann) [1.1]
	- numpy (Travis E. Oliphant et al.) [1.15.1]
	- HTSeq (Anders et al. 2014) [0.11.0]
	- cutadapt (Martin 2011) [1.18]

------------------------------------------------------------------------------------
EOF

# ----------------------
# MAIN
# ----------------------
[[ $# -eq 0 ]] && { printf "%s\n" "$usage"; exit 0; } 		# if no user-supplied arguments, print usage and exit

# ----------------------
# Get user-specified arguments
# ----------------------

# Initiate environment
scriptDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )		# location of this script 
workdir=$( pwd )												# current working directory

# Required arguments:
# ----------------------
indir=""										# folder containing all the raw data from all libraries
descfile=""							# file containing list of samples (see above for expected fields)
outdir=""							# where to put all output files from this project

# Additional arguments:
# ----------------------
path_to_pipelines="$scriptDir"							# path to folder containing all required pipeline scripts (see list below)
path_to_scripts="$scriptDir/scripts"					# path to folder containing all required helper scripts (see list below)
path_to_files="$scriptDir/required_files"				# path to folder containing all required other files (see list below)

# Flag options:
# ----------------------
skip1=false							# skip step 1
skip2=false							# skip step 2
skip3=false							# skip step 3
skip4=false							# skip step 4
skip5=false							# skip step 5
skip6=false							# skip step 6
skip7=false							# skip step 7
skip8=false							# skip step 8
skip9=false

# ----------------------
while getopts "i:d:o:p:s:f:234567890h" opt; do
	case $opt in
		i)	# folder containing all the raw data from all libraries (see above for expected format)
			indir="$OPTARG"
			;;
		d)	# file containing list of samples (see above for expected fields)
			descfile="$OPTARG"
			;;
		o)	# where to put all output files from this project
			outdir="$OPTARG"
			;;
		p)	# path to folder containing all required pipeline scripts (see list below)
			path_to_pipelines="$OPTARG"
			;;
		s)	# path to folder containing all required helper scripts (see list below)
			path_to_scripts="$OPTARG"
			;;
		f)	# path to folder containing all required other files (see list below)
			path_to_files="$OPTARG"
			;;
		2)	# skip to step 2
			skip1=true
			;;
		3)	# skip to step 3
			skip1=true
			skip2=true
			;;
		4)	# skip to step 4
			skip1=true
			skip2=true
			skip3=true
			;;
		5)	# skip to step 5
			skip1=true
			skip2=true
			skip3=true
			skip4=true
			;;
		6)	# skip to step 6
			skip1=true
			skip2=true
			skip3=true
			skip4=true
			skip5=true
			;;
		7)	# skip to step 7
			skip1=true
			skip2=true
			skip3=true
			skip4=true
			skip5=true
			skip6=true
			;;
		8)	# skip to step 8
			skip1=true
			skip2=true
			skip3=true
			skip4=true
			skip5=true
			skip6=true
			skip7=true
			;;
		9)	# skip to step 9
			skip1=true
			skip2=true
			skip3=true
			skip4=true
			skip5=true
			skip6=true
			skip7=true
			skip8=true
			;;
		0)	# check dependencies ok then exit
			checkdep=true
			;;
		h)	# print usage and version information to stdout and exit
			echo "$usage"
			exit 0
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

# Check all dependencies can be located
# ----------------------

# Check that all other required programs on $PATH are installed
# ----------------------
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools is required on PATH but was not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools is required on PATH but was not found"; exit 1; }
command -v python >/dev/null 2>&1 || { echo "Error: python is required on PATH but was not found"; exit 1; }
command -v R >/dev/null 2>&1 || { echo "Error: R is required on PATH but was not found"; exit 1; }

# Python must be version 2.x.x
pver=$( python --version | cut -f2 -d ' ' | cut -f1 -d'.' )
pverfull=$( python --version )
[ "$pver" -eq 2 ] || { echo "Error: this script requires python v.2.x.x, installed version is $pver"; exit 1; }

# Check that required pipelines are in path_to_pipelines (set custom path to this location with option -p)
# ----------------------
[ ! -f "$path_to_pipelines/rna_seq_map.sh" ] && { echo "Error: could not find required script rna_seq_se_map.sh in provided folder (${path_to_pipelines}), change this location with -p"; exit 1; }

# Check all pipeline dependencies
# ----------------------
$path_to_pipelines/rna_seq_map.sh -0

# Check that required scripts are in path_to_scripts (set custom path to this location with option -s)
# ----------------------
required_scripts_list=( 
"MarkDuplicates.jar"
"assign_to_allele.py"
"barchart.R"
"cluster_gene_expression.R"
"gene_overlaps_dotplot.R"
"merge_by_column.R"
"merge_many_files.sh"
"piechart.R"
"plot_SC3_heatmap.py"
"run_DEsingle.R"
"run_topGO.R"
"scatterplot.R"
"single_cell_ASE_analysis.R"
"single_cell_ASE_src.R"
"single_cell_RNAseq_plots.R"
"single_cell_cluster_PCA_tSNE.R"
"single_cell_cluster_SC3.R"
"single_cell_simulate_reads.R"
"single_cell_trajectory_analysis.R"
"subset_large_file.R"
"check_dependencies.R"
)

for ((i=0;i<${#required_scripts_list[@]};++i)); do
[ ! -f "$path_to_scripts/${required_scripts_list[i]}" ] && { echo "Error: could not find required script ${required_scripts_list[i]} in provided folder (${path_to_scripts}), change this location with -s"; exit 1; }
done

# Check that all other required files are in path_to_files (set custom path to this location with option -f)
# ----------------------
required_files_list=( 
"TAIR10_plus_araport11_nonoverlapping.gtf"
"araport11_TEs_not_in_genes.gtf"
"ERCC_spikeins/ERCC92.gtf"
"genomes/TAIR10_Cvi_metachrom.txt"
"Col_Cvi_SNPs_clean.bed"
"genomes/TAIR10_Ler_metachrom.txt"
"Col_Ler_SNPs_clean.bed"
"ERCC_chrlist.txt" 
)

for ((i=0;i<${#required_files_list[@]};++i)); do
[ ! -f "$path_to_files/${required_files_list[i]}" ] && { echo "Error: could not find required file ${required_files_list[i]} in provided folder (${path_to_files}), change this location with -f"; exit 1; }
done

# Check all required R packages are installed
# ----------------------
$path_to_scripts/check_dependencies.R > /dev/null || { echo "Not all required R packages are installed, exiting."; exit 1; }

# Check all required python libraries are installed
# ----------------------
pip show argparse -q 2> /dev/null || { echo "Required python package argparse not installed, exiting."; exit 1; }
pip show numpy -q 2> /dev/null || { echo "Required python package numpy not installed, exiting."; exit 1; }
pip show cutadapt -q 2> /dev/null || { echo "Required python package cutadapt not installed, exiting."; exit 1; }
pip show HTSeq -q 2> /dev/null || { echo "Required python package HTSeq not installed, exiting."; exit 1; }

# Done checking all dependencies. Stop here if -0 flagged.
# ----------------------
"$checkdep" && exit 0


# Check all required user inputs are provided
# ----------------------
[ "$skip1" = "false" ] && [ -z "$indir" ] && { echo "Error: -i indir is a required argument (folder containing all the raw data from all libraries (see above for expected format))"; exit 1; }
[ -z "$descfile" ] && { echo "Error: -d descfile is a required argument (file containing list of samples (see above for expected fields))"; exit 1; }
[ -z "$outdir" ] && { echo "Error: -o outdir is a required argument (where to put all output files from this project)"; exit 1; }

# ----------------------
# Parameters to use:
# ----------------------
# Annotation files:
gtf_genes="$path_to_files/TAIR10_plus_araport11_nonoverlapping.gtf"
gtf_TEs="$path_to_files/araport11_TEs_not_in_genes.gtf"
gtf_ercc="$path_to_files/ERCC_spikeins/ERCC92.gtf"

# FOR STEP 1 (preprocessing and mapping):
ColCvimetagenome="$path_to_files/genomes/TAIR10_Cvi_ERCC_OH50"				# TAIR10 + Cvi (SNP subs only) + ERCC, with --sjdbOverhang 50
ColCvimetachrom="$path_to_files/genomes/TAIR10_Cvi_metachrom.txt"
ColCvisnpfile="$path_to_files/Col_Cvi_SNPs_clean.bed"
ColLermetagenome="$path_to_files/genomes/TAIR10_Ler_ERCC_OH50"			# TAIR10 + Ler (SNP subs only) + ERCC, with --sjdbOverhang 50
ColLermetachrom="$path_to_files/genomes/TAIR10_Ler_metachrom.txt"
ColLersnpfile="$path_to_files/Col_Ler_SNPs_clean.bed"
spikeinchrs="$path_to_files/ERCC_chrlist.txt"
trim=0															# trim this many bases of 5' end of all reads
frac_mismatch=0.05												# max # mismatches (tophat -N option)
minIntron=20													# min intron size (see tophat -i option)
maxIntron=5000													# max intron size (see tophat -I option)
maxdist=5000													# max distance between two mates in order to report
adapter="CTGTCTCTTATA"											# nextera xt adapter

# FOR STEP 2: (htseq-count)

# FOR STEP 3: (quality filtering)
mingenescov=1500		# min # of genes that must be covered
mingenescov5=1000		# min # of genes that must have at least 5 coverage

# FOR STEP 6 (identifying DE genes)
pvalcutoff=0.0001		# p-value cutoff
log2fccutoff=2			# log2 fold change cutoff
pvalcutoff_sc=0.001		# less stringent p-value cutoff for identifying DE genes specifically between SC and endosperm
log2fccutoff_sc=1		# less stringent log2 fold change cutoff for identifying DE genes specifically between SC and endosperm


# ----------------------
# HELPER FUNCTIONS 
# ----------------------
err_msg () {
  	printf "Error: $1 \n" | tee -a "$2"
  	echo "Exiting."
  	exit 1
}

displaytime () {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  [[ $D > 0 ]] && printf '%d days ' $D
  [[ $H > 0 ]] && printf '%d hours ' $H
  [[ $M > 0 ]] && printf '%d minutes ' $M
  [[ $D > 0 || $H > 0 || $M > 0 ]] && printf 'and '
  printf '%d seconds\n' $S
}

rna_seq_map_summary () {
	info="$3"
	strain2="$4"
	stub=`awk '/Stubname/ {print $2}' $1`
	readlen=`awk '/Reads are/ {print $3}' $1`
	libtype=`awk '/Library is stranded/ {print $4}' $1`
	encoding=`awk '/encoding/ {print $5}' $1`
	trim=`awk '/bases to trim at/ {print $11}' $1`
	qual=`awk '/Quality cutoff/ {print $8}' $1`
	fmismatch=`awk '/Max number of mismatches/ {print $12}' $1`

	# get info about mapping
	n_tot=$( awk '/Input file contained/ {print $4}' $1 )
	n_qc=$( awk '/passed quality filtering steps/ {print $1}' $1 )

	uniq_Col=$(awk '/mapped uniquely to Col/ {print $1}' $1 )
	[ "$strain2" = "Cvi" ] && uniq_alt=$(awk '$0 ~ /mapped uniquely to Cvi/ {print $1}' $1 )
	[ "$strain2" = "Ler" ] && uniq_alt=$(awk '$0 ~ /mapped uniquely to Ler/ {print $1}' $1 )
	uniq_none=$(awk '/mapped uniquely to both strains/ {print $1}' $1 )
	overall_uniq=$(( $uniq_Col + $uniq_alt + $uniq_none ))
	
	uniq_dedup_Col=$(awk '/mapped uniquely after PCR dedup to Col/ {print $1}' $1 )
	[ "$strain2" = "Cvi" ] && uniq_dedup_alt=$(awk '$0 ~ /mapped uniquely after PCR dedup to Cvi/ {print $1}' $1 )
	[ "$strain2" = "Ler" ] && uniq_dedup_alt=$(awk '$0 ~ /mapped uniquely after PCR dedup to Ler/ {print $1}' $1 )
	uniq_dedup_none=$(awk '/uniquely after PCR dedup to both strains/ {print $1}' $1 )
	overall_uniq_dedup=$(( $uniq_dedup_Col + $uniq_dedup_alt + $uniq_dedup_none ))
	
	tot_spikein=$(awk '/mapped uniquely to the spike-in genome/ {print $1}' $1 )
	tot_spikein_dedup=$(awk '/uniquely after PCR dedup to the spike-in genome/ {print $1}' $1 )
				
	printf "$stub\t${info}\t$readlen\t$libtype\t$encoding\t$trim\t$qual\t$fmismatch\t" >> $2
	printf "$n_tot\t$n_qc\t$overall_uniq\t$uniq_Col\t$uniq_alt\t$uniq_none\t" >> $2
	printf "$overall_uniq_dedup\t$uniq_dedup_Col\t$uniq_dedup_alt\t$uniq_dedup_none\t" >> $2
	printf "$tot_spikein\t$tot_spikein_dedup\n" >> $2
}


# ----------------------
# RUN PIPELINE
# ----------------------
mkdir -p "$outdir/LSF_logs"; lsf="$outdir/LSF_logs"
mkdir -p "$outdir/_summary_stats"; summ="$outdir/_summary_stats"
log="$outdir/logfile.txt"
time_start=$(date)	# time run was started
time_ss=$(date +%s)	# time run was started (in seconds)

# Output user-derived options to stdout and to log file
# ----------------------
echo "Running snRNAseq_analysis.sh v1.0 (03/31/2021):" | tee "$log"
echo "Run start on: $time_start" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Working directory: $( pwd )" | tee -a "$log"
echo "Output directory: $outdir" | tee -a "$log"
echo "Sample descriptions file: $descfile" | tee -a "$log"
[ "$skip1" = "false" ] && { echo "Raw fastq files are in: $indir" | tee -a "$log"; }
echo "Log file: $log" | tee -a "$log"
echo "LSF job logs are in: $outdir/LSF_logs" | tee -a "$log"
echo "Pipelines are in: $path_to_pipelines" | tee -a "$log"
echo "Helper scripts are in: $path_to_scripts" | tee -a "$log"
echo "Required files are in: $path_to_files" | tee -a "$log"
echo "-------------------------" | tee -a "$log"

# Info on libraries is stored in descfile; extract relevant info and save
# ----------------------
i=0
libID=()
stubname=()
project=()
seqtype=()
cross=()
dap=()
num_nuc=()
peak=()
stage=()
stagecode=()
sbatch=()
strain2=()
platform=()
ercc=()

while IFS=$'\t' read -r -a aa; do
	libID[i]="${aa[0]}"
	stubname[i]="${aa[4]}"
	project[i]="${aa[1]}"
	seqtype[i]="${aa[5]}"
	cross[i]="${aa[6]}"
	dap[i]="${aa[7]}"
	num_nuc[i]="${aa[8]}"
	peak[i]="${aa[9]}"
	sbatch[i]="${aa[12]}"
	stage[i]="${aa[10]}"
	stagecode[i]="${aa[11]}"
	platform[i]="${aa[13]}"
	ercc[i]="${aa[14]}"
	[[ ${cross[i]} == *"L"* ]] && strain2[i]="Ler"
	[[ ${cross[i]} == *"V"* ]] && strain2[i]="Cvi"
	i=$(( $i + 1 ))
done < <( tail -n+2 "$descfile" )

# Assert all raw data files can be found (if running step 1)
if [ "$skip1" = "false" ]; then
	for ((i=0;i<${#stubname[@]};++i)); do
		if [ "${seqtype[i]}" = "paired" ]; then
			# data are paired-end
			[ -f "$indir/${stubname[i]}_R1.fastq" ] || { echo "Could not find fastq file $indir/${stubname[i]}_R1.fastq"; exit 1; }
			[ -f "$indir/${stubname[i]}_R2.fastq" ] || { echo "Could not find fastq file $indir/${stubname[i]}_R2.fastq"; exit 1; }
		else
			# data are single-end
			[ -f "$indir/${stubname[i]}.fastq" ] || { echo "Could not find fastq file $indir/${stubname[i]}.fastq"; exit 1; }
		fi
	done
fi
		
# Print out summary of libraries in dataset, from descfile
echo "" | tee -a "$log"
echo "Summary of all libraries in this project:" | tee -a "$log"
echo "Total libraries: ${#libID[@]}" | tee -a "$log"
numnocellctrl=0; numblankctrl=0; numtwonuc=0; numsingle=0
numCxV=0; numVxC=0; numCxL=0; numLxC=0
num3N=0; num6N=0
num2dap=0; num3dap=0; num4dap=0; num5dap=0

for ((i=0;i<${#libID[@]};++i)); do
	[ "${dap[i]}" == "blank" ] && numblankctrl=$(( $numblankctrl + 1 ))
	[ "${dap[i]}" == "nocell" ] && numnocellctrl=$(( $numnocellctrl + 1 ))
	[ "${num_nuc[i]}" == "2" ] && numtwonuc=$(( $numtwonuc + 1 ))
	[ "${num_nuc[i]}" == "1" ] && numsingle=$(( $numsingle + 1 ))

	[ "${num_nuc[i]}" == "1" ] && [ "${cross[i]}" == "CxV" ] && numCxV=$(( $numCxV + 1 ))
	[ "${num_nuc[i]}" == "1" ] && [ "${cross[i]}" == "VxC" ] && numVxC=$(( $numVxC + 1 ))
	[ "${num_nuc[i]}" == "1" ] && [ "${cross[i]}" == "CxL" ] && numCxL=$(( $numCxL + 1 ))
	[ "${num_nuc[i]}" == "1" ] && [ "${cross[i]}" == "LxC" ] && numLxC=$(( $numLxC + 1 ))
	
	[ "${num_nuc[i]}" == "1" ] && [ "${peak[i]}" == "3N" ] && num3N=$(( $num3N + 1 ))
	[ "${num_nuc[i]}" == "1" ] && [ "${peak[i]}" == "6N" ] && num6N=$(( $num6N + 1 ))
	
	[ "${num_nuc[i]}" == "1" ] && [ "${dap[i]}" == "2" ] && num2dap=$(( $num2dap + 1 ))
	[ "${num_nuc[i]}" == "1" ] && [ "${dap[i]}" == "3" ] && num3dap=$(( $num3dap + 1 ))
	[ "${num_nuc[i]}" == "1" ] && [ "${dap[i]}" == "4" ] && num4dap=$(( $num4dap + 1 ))
	[ "${num_nuc[i]}" == "1" ] && [ "${dap[i]}" == "5" ] && num5dap=$(( $num5dap + 1 ))
done

echo " -> $numblankctrl are \"blank\" controls (water used for Nextera prep)" | tee -a "$log"
echo " -> $numnocellctrl are \"nocell\" controls (no nucleus sorted into well)" | tee -a "$log"
echo " -> $numtwonuc are \"two nuclei\" controls (two 3N nuclei sorted into well)" | tee -a "$log"
echo " -> $numsingle are \"single nuclei\" samples" | tee -a "$log"
echo "" | tee -a "$log"
echo "Of the $numsingle single nuclei samples:" | tee -a "$log"
echo " -> $numCxV are from a Col x Cvi cross" | tee -a "$log"
echo " -> $numVxC are from a Cvi x Col cross" | tee -a "$log"
echo " -> $numCxL are from a Col x Ler cross" | tee -a "$log"
echo " -> $numLxC are from a Ler x Col cross" | tee -a "$log"
echo "" | tee -a "$log"
echo "Of the $numsingle single nuclei samples:" | tee -a "$log"
echo " -> $num3N are from the 3N ploidy peak" | tee -a "$log"
echo " -> $num6N are from the 6N ploidy peak" | tee -a "$log"
echo "" | tee -a "$log"
echo "Of the $numsingle single nuclei samples:" | tee -a "$log"
echo " -> $num2dap are from seeds harvested at 2 DAP" | tee -a "$log"
echo " -> $num3dap are from seeds harvested at 3 DAP" | tee -a "$log"
echo " -> $num4dap are from seeds harvested at 4 DAP" | tee -a "$log"
echo " -> $num5dap are from seeds harvested at 5 DAP" | tee -a "$log"
echo "" | tee -a "$log"


# ----------------------
# Step 1: run rna_seq_se_map or rna_seq_pe_map on each library
# ----------------------
if "$skip1"; then printf "\nSkipping step 1 (mapping reads) by user request...\n" | tee -a "$log"
else printf "\nStep 1: mapping reads\n" | tee -a "$log"; ts=$(date +%s)
			
	# filter and map reads for each library
	# ---------------
	# note that all libraries are unstranded and encoded with Sanger/Illumina 1.9, which is Phred+33 not Phred+64		
	pid=(); mkdir -p "$lsf/rna_seq_map"
	for ((i=0;i<${#stubname[@]};++i)); do
		mkdir -p "$outdir/${project[i]}_map"
		
		# determine if cross is Col to Ler or to Cvi and get corresponding metagenome
		if [[ ${strain2[i]} == "Ler" ]]; then
			metagenometouse="$ColLermetagenome"
			metachromtouse="$ColLermetachrom"
			snpfiletouse="$ColLersnpfile"
		else
			metagenometouse="$ColCvimetagenome"
			metachromtouse="$ColCvimetachrom"
			snpfiletouse="$ColCvisnpfile"
		fi	
		
		# get input files for this library
		[ "${seqtype[i]}" = "paired" ] && infilestr="-1 $indir/${stubname[i]}_R1.fastq -2 $indir/${stubname[i]}_R2.fastq" || pestr="-1 $indir/${stubname[i]}.fastq"
					
		# make and run mapping command
		cmd="$path_to_pipelines/rna_seq_map.sh ${infilestr} -o $outdir/${project[i]}_map/${stubname[i]}_map -n ${stubname[i]} -g $metagenometouse -C $metachromtouse -A Col -B ${strain2[i]} -a $adapter -i $minIntron -I $maxIntron -t $trim -N $frac_mismatch -d $maxdist -P $spikeinchrs -3r"
		bsub -o "$lsf/rna_seq_map/${stubname[i]}.txt" -K "$cmd" & pid[i]=$!
	done
		
	# wait for all jobs to finish
	for ((i=0;i<${#stubname[@]};++i)); do
		wait "${pid[i]}" || err_msg "read mapping failed for ${stubname[i]}, see $lsf/rna_seq_map/${stubname[i]}.txt" "$log"	
	done
	
	# use assign_to_allele to assign reads to parent-of-origin explicitly using SNPs
	# ---------------
	pid=()
	for ((i=0;i<${#stubname[@]};++i)); do
		mkdir -p "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele"
		[ ${strain2[i]} = "Ler" ] && snpfile="$ColLersnpfile" || snpfile="$ColCvisnpfile"
		samtools view -ho "$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments.sam" "$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments.bam"
		
		cmd="$path_to_scripts/assign_to_allele.py $snpfile $outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments.sam $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]} --refname Col --altname ${strain2[i]}"
		bsub -o "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_log.txt" -K "$cmd" & pid[i]=$!
	done
	
	for ((i=0;i<${#stubname[@]};++i)); do
		wait "${pid[i]}" || err_msg "error running assign_to_allele, see $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_log.txt" "$log"	
		rm "$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments.sam"
	done

	# grab generic SAM header
	i=0
	samtools view -H "$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments.bam" | head -n -2 > "$outdir/SAM_header.txt"
	
	# sort and compress all files; bsub in groups of 100 to speed things up
	mkdir -p "$outdir/commands_tmp"
	pids=(); j=0
	echo '#!/bin/bash' > "$outdir/commands_tmp/sort_compress_${j}.sh"
	echo "" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
	for ((i=0;i<${#stubname[@]};++i)); do
		[ "$(($i % 100))" = "0" ] && [ "$i" != 0 ] && { echo "exit 0" >> "$outdir/commands_tmp/sort_compress_${j}.sh"; chmod a+x "$outdir/commands_tmp/sort_compress_${j}.sh"; rm -f "$lsf/sort_compress_${j}.txt"; bsub -o "$lsf/sort_compress_${j}.txt" -K "$outdir/commands_tmp/sort_compress_${j}.sh" & pids[j]=$!; j=$(( $j + 1 )); echo '#!/bin/bash' > "$outdir/commands_tmp/sort_compress_${j}.sh"; echo "" >> "$outdir/commands_tmp/sort_compress_${j}.sh"; }
		echo "[ -f $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col.sam ] || cat $outdir/SAM_header.txt > $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col.sam" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
		cc1="samtools sort $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col.sam -o $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col.bam -T $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
		echo "$cc1" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
		echo "[ \$? != 0 ] && { echo \"Error in cmd: $cc1\"; exit 1; }" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
		echo "[ -f $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}.sam ] || cat $outdir/SAM_header.txt > $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}.sam" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
		cc2="samtools sort $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}.sam -o $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}.bam -T $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}"
		echo "$cc2" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
		echo "[ \$? != 0 ] && { echo \"Error in cmd: $cc2\"; exit 1; }" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
		echo "[ -f $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none.sam ] || cat $outdir/SAM_header.txt > $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none.sam" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
		cc3="samtools sort $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none.sam -o $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none.bam -T $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none"
		echo "$cc3" >> "$outdir/commands_tmp/sort_compress_${j}.sh"
		echo "[ \$? != 0 ] && { echo \"Error in cmd: $cc3\"; exit 1; }" >> "$outdir/commands_tmp/sort_compress_${j}.sh"	
	done
	
	# submit final job
	echo "exit 0" >> "$outdir/commands_tmp/sort_compress_${j}.sh"; chmod a+x "$outdir/commands_tmp/sort_compress_${j}.sh"; rm -f "$lsf/sort_compress_${j}.txt"; bsub -o "$lsf/sort_compress_${j}.txt" -K "$outdir/commands_tmp/sort_compress_${j}.sh" & pids[j]=$!
	
	for ((j=0;j<${#pids[@]};++j)); do
		wait "${pids[j]}" || err_msg "error with samtools, see $lsf/sort_compress_${j}.txt" "$log"	
	done

	rm -rf "$outdir/commands_tmp/*"
	pids=(); j=0
	echo '#!/bin/bash' > "$outdir/commands_tmp/markduplicates_${j}.sh"
	echo "" >> "$outdir/commands_tmp/markduplicates_${j}.sh"
		 	
	# remove PCR duplicates using MarkDuplicates
	for ((i=0;i<${#stubname[@]};++i)); do
		[ "$(($i % 100))" = "0" ] && [ "$i" != 0 ] && { echo "exit 0" >> "$outdir/commands_tmp/markduplicates_${j}.sh"; chmod a+x "$outdir/commands_tmp/markduplicates_${j}.sh"; rm -f "$lsf/markduplicates_${j}.txt"; bsub -o "$lsf/markduplicates_${j}.txt" -K "$outdir/commands_tmp/markduplicates_${j}.sh" & pids[j]=$!; j=$(( $j + 1 )); echo "#!/bin/bash" > "$outdir/commands_tmp/markduplicates_${j}.sh"; echo "" >> "$outdir/commands_tmp/markduplicates_${j}.sh"; }

		cmd1="java -Xmx10g -jar $path_to_scripts/MarkDuplicates.jar I=$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col.bam O=$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col_dedup.bam METRICS_FILE=$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col_dedup_log.txt REMOVE_DUPLICATES=true"
		cmd2="java -Xmx10g -jar $path_to_scripts/MarkDuplicates.jar I=$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}.bam O=$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}_dedup.bam METRICS_FILE=$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}_dedup_log.txt REMOVE_DUPLICATES=true"
		cmd3="java -Xmx10g -jar $path_to_scripts/MarkDuplicates.jar I=$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none.bam O=$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none_dedup.bam METRICS_FILE=$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none_dedup_log.txt REMOVE_DUPLICATES=true"

		c1=$( samtools view -c "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col.bam" )
		[ "$c1" != 0 ] && { echo "$cmd1" >> "$outdir/commands_tmp/markduplicates_${j}.sh"; echo "[ \$? != 0 ] && { echo \"Error in cmd: $cmd1\"; exit 1; }" >> "$outdir/commands_tmp/markduplicates_${j}.sh"; }
		[ "$c1" != 0 ] || cat "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col.bam" > "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col_dedup.bam"
		c2=$( samtools view -c "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}.bam" )
		[ "$c2" != 0 ] && { echo "$cmd2" >> "$outdir/commands_tmp/markduplicates_${j}.sh"; echo "[ \$? != 0 ] && { echo \"Error in cmd: $cmd2\"; exit 1; }" >> "$outdir/commands_tmp/markduplicates_${j}.sh"; }
		[ "$c2" != 0 ] || cat "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}.bam" > "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}_dedup.bam"
		c3=$( samtools view -c "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none.bam" )
		[ "$c3" != 0 ] && { echo "$cmd3" >> "$outdir/commands_tmp/markduplicates_${j}.sh"; echo "[ \$? != 0 ] && { echo \"Error in cmd: $cmd3\"; exit 1; }" >> "$outdir/commands_tmp/markduplicates_${j}.sh"; }
		[ "$c3" != 0 ] || cat "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none.bam" > "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none_dedup.bam"
	done
	
	echo "exit 0" >> "$outdir/commands_tmp/markduplicates_${j}.sh"; chmod a+x "$outdir/commands_tmp/markduplicates_${j}.sh"; rm -f "$lsf/markduplicates_${j}.txt"; bsub -o "$lsf/markduplicates_${j}.txt" -K "$outdir/commands_tmp/markduplicates_${j}.sh" & pids[j]=$!

	for ((j=0;j<${#pids[@]};++j)); do
		wait "${pids[j]}" || err_msg "error with MarkDuplicates, see $lsf/markduplicates_${j}.txt" "$log"	
	done
	
	# merge the dedup'd files back together
	rm -rf "$outdir/commands_tmp/*"
	pids=(); j=0
	echo "#!/bin/bash" > "$outdir/commands_tmp/merge_${j}.sh"
	echo "" >> "$outdir/commands_tmp/merge_${j}.sh"
		 	
	for ((i=0;i<${#stubname[@]};++i)); do
		[ "$(($i % 100))" = "0" ] && [ "$i" != 0 ] && { echo "exit 0" >> "$outdir/commands_tmp/merge_${j}.sh"; chmod a+x "$outdir/commands_tmp/merge_${j}.sh"; rm -f "$lsf/merge_${j}.txt"; bsub -o "$lsf/merge_${j}.txt" -K "$outdir/commands_tmp/merge_${j}.sh" & pids[j]=$!; j=$(( $j + 1 )); echo "#!/bin/bash" > "$outdir/commands_tmp/merge_${j}.sh"; echo "" >> "$outdir/commands_tmp/merge_${j}.sh"; }
		cmd="samtools merge -f $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_all_dedup.bam $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_Col_dedup.bam $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${strain2[i]}_dedup.bam $outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_none_dedup.bam"
		echo "$cmd" >> "$outdir/commands_tmp/merge_${j}.sh"; echo "[ \$? != 0 ] && { echo \"Error in cmd: $cmd\"; exit 1; }" >> "$outdir/commands_tmp/merge_${j}.sh"
	done
	
	echo "exit 0" >> "$outdir/commands_tmp/merge_${j}.sh"; chmod a+x "$outdir/commands_tmp/merge_${j}.sh"; rm -f "$lsf/merge_${j}.txt"; bsub -o "$lsf/merge_${j}.txt" -K "$outdir/commands_tmp/merge_${j}.sh" & pids[j]=$!

	for ((j=0;j<${#pids[@]};++j)); do
		wait "${pids[j]}" || err_msg "error with samtools merge, see $lsf/merge_${j}.txt" "$log"	
	done
	
	# summarize results	(Table S1)		
	# ---------------
	printf "stubname\tcross\tDAP\tN_nuc\tpeak\treadlen\tlibtype\tencoding\ttrim\tqualencode\tf_mismatch\t" > "$summ/rna_seq_map_summary.txt"
	printf "n_tot\tn_qc\toverall_uniq\tuniq_Col\tuniq_alt\tuniq_none\t" >> "$summ/rna_seq_map_summary.txt"
	printf "overall_uniq_dedup\tuniq_dedup_Col\tuniq_dedup_alt\tuniq_dedup_none\t" >> "$summ/rna_seq_map_summary.txt"
	printf "tot_spikein\ttot_spikein_dedup\n" >> "$summ/rna_seq_map_summary.txt"
	
	for ((i=0;i<${#stubname[@]};++i)); do
		info=$( printf "${cross[i]}\t${dap[i]}\t${num_nuc[i]}\t${peak[i]}" )
		rna_seq_map_summary "$outdir/${project[i]}_map/${stubname[i]}_map/${stubname[i]}_log.txt" "$summ/rna_seq_map_summary.txt" "$info" "${strain2[i]}"
	done		
	
	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	

fi


# ----------------------
# Step 2: get read counts and CPMs for each gene in each library
# ----------------------
if "$skip2"; then printf "\nSkipping step 2 (getting per gene counts and FPKMs) by user request...\n" | tee -a "$log"
else printf "\nStep 2: getting per gene counts and FPKMs\n" | tee -a "$log"; ts=$(date +%s)
		
	# run htseq-count on all libraries (all reads, Col reads only, alt reads only, both with and without PCR duplicates), over genes, TEs
	# this is done in batches of 20 commands together, which are saved together into scripts in $outdir/commands_tmp
	# ---------------
	pids=(); j=0
	for type in "w_dup" "dedup"; do
		echo "#!/bin/bash" > "$outdir/commands_tmp/htseq_${type}_${j}.sh"
		echo "" >> "$outdir/commands_tmp/htseq_${type}_${j}.sh"
			
		for ((i=0;i<${#stubname[@]};++i)); do
			[ "$(($i % 20))" = "0" ] && [ "$i" != 0 ] && { echo "exit 0" >> "$outdir/commands_tmp/htseq_${type}_${j}.sh"; chmod a+x "$outdir/commands_tmp/htseq_${type}_${j}.sh"; rm -f "$lsf/htseq_${type}_${j}.txt"; bsub -o "$lsf/htseq_${type}_${j}.txt" -K "$outdir/commands_tmp/htseq_${type}_${j}.sh" & pids+=( $! ); j=$(( $j + 1 )); echo "#!/bin/bash" > "$outdir/commands_tmp/htseq_${type}_${j}.sh"; echo "" >> "$outdir/commands_tmp/htseq_${type}_${j}.sh"; }
			echo "echo \"Counting $type reads in ${stubname[i]}\"" >> "$outdir/commands_tmp/htseq_${type}_${j}.sh"
			loc="$outdir/${project[i]}_map/${stubname[i]}_map/htseq_count"
			mkdir -p "$loc"
			for allele in "all" "Col" "${strain2[i]}"; do	
				[ "$type" = "w_dup" ] && ff="$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${allele}.bam" || ff="$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_${allele}_dedup.bam"
				[ "$type" = "w_dup" ] && [ "$allele" = "all" ] && ff="$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments.bam"
				if [ $( samtools view -c "$ff" ) -ne 0 ]; then
					cmd="htseq-count -r pos --nonunique none -m intersection-nonempty -f bam -s no $ff $gtf_genes > $loc/${stubname[i]}_counts_${type}_${allele}_genes.txt"
					echo "$cmd" >> "$outdir/commands_tmp/htseq_${type}_${j}.sh"; echo "[ \$? != 0 ] && { echo \"Error in cmd: $cmd\"; exit 1; }" >> "$outdir/commands_tmp/htseq_${type}_${j}.sh"

					cmd="htseq-count -r pos --nonunique none -m intersection-nonempty -t transposable_element -f bam -s no $ff $gtf_TEs > $loc/${stubname[i]}_counts_${type}_${allele}_TEs.txt"
					echo "$cmd" >> "$outdir/commands_tmp/htseq_${type}_${j}.sh"; echo "[ \$? != 0 ] && { echo \"Error in cmd: $cmd\"; exit 1; }" >> "$outdir/commands_tmp/htseq_${type}_${j}.sh"
				fi
			done
		done
		echo "exit 0" >> "$outdir/commands_tmp/htseq_${type}_${j}.sh"; chmod a+x "$outdir/commands_tmp/htseq_${type}_${j}.sh"; rm -f "$lsf/htseq_${type}_${j}.txt"; bsub -o "$lsf/htseq_${type}_${j}.txt" -K "$outdir/commands_tmp/htseq_${type}_${j}.sh" & pids+=( $! )
	done

	for ((j=0;j<${#pids[@]};++j)); do
		wait "${pids[j]}" || { [ $(( $j % 2 )) = 0 ] && { err_msg "error with htseq-count, see $lsf/htseq_w_dup_${j}.txt" "$log"; } || { err_msg "error with htseq-count, see $lsf/htseq_dedup_${j}.txt" "$log"; }; }
	done
	
	# any samples where there were no reads will not have generated files above; generate them now with all zeros for read counts (just to simplify some looping later)
	i=0
	loc="$outdir/${project[i]}_map/${stubname[i]}_map/htseq_count"
	cut -f1 "$loc/${stubname[i]}_counts_w_dup_all_genes.txt" > "$outdir/genelist.txt"
	awk -F$'\t' '{OFS=FS} {print $0,"0"}' "$outdir/genelist.txt" > "$outdir/genelist_emptycounts.txt"

	for type in "w_dup" "dedup"; do
		for ((i=0;i<${#stubname[@]};++i)); do
			loc="$outdir/${project[i]}_map/${stubname[i]}_map/htseq_count"
			for allele in "all" "Col" "${strain2[i]}"; do	
				[ ! -f "$loc/${stubname[i]}_counts_${type}_${allele}_genes.txt" ] && { echo "${stubname[i]} $type $allele missing"; cat "$outdir/genelist_emptycounts.txt" > "$loc/${stubname[i]}_counts_${type}_${allele}_genes.txt"; }
				[ ! -f "$loc/${stubname[i]}_counts_${type}_${allele}_TEs.txt" ] && { cat "$outdir/genelist_emptycounts.txt" > "$loc/${stubname[i]}_counts_${type}_${allele}_TEs.txt"; }
				[ ! -f "$loc/${stubname[i]}_counts_${type}_${allele}_intergenic.txt" ] && { cat "$outdir/genelist_emptycounts.txt" > "$loc/${stubname[i]}_counts_${type}_${allele}_intergenic.txt"; }
			done
		done
	done	
			
	# for each library, make a separate file without mito/chloro counts (genes only) and the __no_feature counts at the end
	pids=(); j=0
	for type in "w_dup" "dedup"; do
		echo "#!/bin/bash" > "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"
		echo "" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"
			
		for ((i=0;i<${#stubname[@]};++i)); do
			loc="$outdir/${project[i]}_map/${stubname[i]}_map/htseq_count"
			[ "$(($i % 100))" = "0" ] && [ "$i" != 0 ] && { echo "exit 0" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"; chmod a+x "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"; rm -f "$lsf/rmChrCM_${type}_${j}.txt"; bsub -o "$lsf/rmChrCM_${type}_${j}.txt" -K "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh" & pids+=( $! ); j=$(( $j + 1 )); echo "#!/bin/bash" > "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"; echo "" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"; }

			for allele in "all" "Col" "${strain2[i]}"; do	
				cmd="printf \"locus_name\t${stubname[i]}\n\" > $loc/${stubname[i]}_counts_${type}_${allele}_genes_cleaned.txt; head -n -5 $loc/${stubname[i]}_counts_${type}_${allele}_genes.txt | awk '\$1 !~ /ATM/ && \$1 !~ /ATC/' >> $loc/${stubname[i]}_counts_${type}_${allele}_genes_cleaned.txt"
				echo "$cmd" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"; echo "[ \$? != 0 ] && { echo \"Error in cmd: $cmd\"; exit 1; }" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"
			
				cmd="printf \"locus_name\t${stubname[i]}\n\" > $loc/${stubname[i]}_counts_${type}_${allele}_TEs_cleaned.txt; head -n -5 $loc/${stubname[i]}_counts_${type}_${allele}_TEs.txt | awk '\$1 !~ /ATM/ && \$1 !~ /ATC/' >> $loc/${stubname[i]}_counts_${type}_${allele}_TEs_cleaned.txt"
				echo "$cmd" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"; echo "[ \$? != 0 ] && { echo \"Error in cmd: $cmd\"; exit 1; }" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"

				cmd="printf \"locus_name\t${stubname[i]}\n\" > $loc/${stubname[i]}_counts_${type}_${allele}_intergenic_cleaned.txt; head -n -5 $loc/${stubname[i]}_counts_${type}_${allele}_intergenic.txt | awk '\$1 !~ /ATM/ && \$1 !~ /ATC/' >> $loc/${stubname[i]}_counts_${type}_${allele}_intergenic_cleaned.txt"
				echo "$cmd" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"; echo "[ \$? != 0 ] && { echo \"Error in cmd: $cmd\"; exit 1; }" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"
			done
		done
	
		echo "exit 0" >> "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"; chmod a+x "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh"; rm -f "$lsf/rmChrCM_${type}_${j}.txt"; bsub -o "$lsf/rmChrCM_${type}_${j}.txt" -K "$outdir/commands_tmp/rmChrCM_${type}_${j}.sh" & pids+=( $! )
	done

	for ((j=0;j<${#pids[@]};++j)); do
		wait "${pids[j]}" || { [ $(( $j % 2 )) = 0 ] && { err_msg "error with cleanup of htseq-count, see $lsf/rmChrCM_w_dup_${j}.txt" "$log"; } || { err_msg "error with htseq-count, see $lsf/rmChrCM_dedup_${j}.txt" "$log"; }; }
	done
	
	# get all filenames of counts files and merge together
	rm -rf "$outdir/count_matrices"
	mkdir -p "$outdir/count_matrices/file_lists"
	for type in "w_dup" "dedup"; do
		for ((i=0;i<${#stubname[@]};++i)); do
			loc="$outdir/${project[i]}_map/${stubname[i]}_map/htseq_count"
			for allele in "all" "Col" "${strain2[i]}"; do	
				echo "$loc/${stubname[i]}_counts_${type}_${allele}_genes_cleaned.txt" >> "$outdir/count_matrices/file_lists/${allele}_counts_in_genes_${type}.txt"
				echo "$loc/${stubname[i]}_counts_${type}_${allele}_TEs_cleaned.txt" >> "$outdir/count_matrices/file_lists/${allele}_counts_in_TEs_${type}.txt"
				echo "$loc/${stubname[i]}_counts_${type}_${allele}_intergenic_cleaned.txt" >> "$outdir/count_matrices/file_lists/${allele}_counts_in_intergenic_${type}.txt"
			done
		done
	done
		
	# merge all files together (note - you would be safe doing this with -paste-, which would
	# also be -much- faster; however, it would assume that all count files contain the same genes in the same
	# order; that should always be true with current code, however early versions were less stable and the
	# approach used here instead will always work even if count files contain different sets of genes or are out of order)
	# ---------------
	pids=()
	rm -f "$lsf/merge_many_files_"*
	for allele in "all" "Col" "Cvi" "Ler"; do	
		for type in "w_dup" "dedup"; do
			cmd="$path_to_scripts/merge_many_files.sh $outdir/count_matrices/file_lists/${allele}_counts_in_genes_${type}.txt $outdir/count_matrices/${allele}_counts_in_genes_${type}.txt locus_name"
			bsub -o "$lsf/merge_many_files_${allele}_counts_in_genes_${type}.txt" -K "$cmd" & pids+=( $! )

			cmd="$path_to_scripts/merge_many_files.sh $outdir/count_matrices/file_lists/${allele}_counts_in_TEs_${type}.txt $outdir/count_matrices/${allele}_counts_in_TEs_${type}.txt locus_name"
			bsub -o "$lsf/merge_many_files_${allele}_counts_in_TEs_${type}.txt" -K "$cmd" & pids+=( $! )

			cmd="$path_to_scripts/merge_many_files.sh $outdir/count_matrices/file_lists/${allele}_counts_in_intergenic_${type}.txt $outdir/count_matrices/${allele}_counts_in_intergenic_${type}.txt locus_name"
			bsub -o "$lsf/merge_many_files_${allele}_counts_in_genes_${type}.txt" -K "$cmd" & pids+=( $! )
		done
		
		for ((j=0;j<${#pids[@]};++j)); do
			wait "${pids[j]}" || { err_msg "error with merge_many_files, see $lsf/merge_many_files_*.txt" "$log"; }
		done
		pids=()
	done
			
	# convert raw counts to CPM (note: decided to only use CPM and not TPM because of heavy 3' bias
	# in read distribution, see Fig. S3)
	# ---------------
	# for each sample, transform counts to CPM by dividing by the total counts and multiplying by 1M
	echo "Calculating CPM for each library..."
	for type in "w_dup" "dedup"; do
		printf "libname\tsize_factor\n" > "$outdir/_summary_stats/cpm_size_factors_${type}.txt"
		for ((i=0;i<${#stubname[@]};++i)); do
			printf "Processing ${type} reads in ${stubname[i]}..."
			loc="$outdir/${project[i]}_map/${stubname[i]}_map/cpm"; mkdir -p "$loc"
			
			if [ "$type" = "w_dup" ]; then
				totmapped=$( awk '$0 ~ /mapped uniquely\./ {print $1}' "$outdir/${project[i]}_map/${stubname[i]}_map/${stubname[i]}_log.txt" )
				totspike=$( awk '$0 ~ /mapped uniquely to the spike-in genome/ {print $1}' "$outdir/${project[i]}_map/${stubname[i]}_map/${stubname[i]}_log.txt" )
				totreads=$(( $totmapped - $totspike ))
			else
				totmapped=$( awk '$0 ~ /mapped uniquely after PCR dedup\./ {print $1}' "$outdir/${project[i]}_map/${stubname[i]}_map/${stubname[i]}_log.txt" )
				totspike=$( awk '$0 ~ /mapped uniquely after PCR dedup to the spike-in genome/ {print $1}' "$outdir/${project[i]}_map/${stubname[i]}_map/${stubname[i]}_log.txt" )
				totreads=$(( $totmapped - $totspike ))
			fi	
			echo "total aligned reads = ${totreads}"
			
			ff="$outdir/${project[i]}_map/${stubname[i]}_map/htseq_count/${stubname[i]}_counts_${type}_all_genes_cleaned.txt"
			
			if [ "$totreads" -gt 0 ]; then
				awk -F$'\t' -v a="$totreads" '{OFS=FS} NR==1{print $0} NR!=1{print $1,($2/a)*1000000}' "$ff" > "$loc/${stubname[i]}_cpm_${type}_all_genes.txt"
				sf=$( echo "$totreads" | awk '{print $0 / 1000000}' )
				printf "${stubname[i]}\t$sf\n" >> "$outdir/_summary_stats/cpm_size_factors_${type}.txt"
			else 
				echo "Library ${stubname[i]} had no counts and was not normalized"
				cp "$ff" "$loc/${stubname[i]}_cpm_${type}_all_genes.txt"
				sf=0
				printf "${stubname[i]}\t$sf\n" >> "$outdir/_summary_stats/cpm_size_factors_${type}.txt"
			fi			
		done
	done	
	
	# merge together all counts into single matrix again
	rm -f "$outdir/count_matrices/file_lists/all_cpm_in_genes_w_dup.txt" "$outdir/count_matrices/file_lists/all_cpm_in_genes_dedup.txt" $outdir/count_matrices/all_cpm_in_genes_w_dup.txt $outdir/count_matrices/all_cpm_in_genes_dedup.txt
	for type in "w_dup" "dedup"; do
		for ((i=0;i<${#stubname[@]};++i)); do
			loc="$outdir/${project[i]}_map/${stubname[i]}_map/cpm"
			echo "$loc/${stubname[i]}_cpm_${type}_all_genes.txt" >> "$outdir/count_matrices/file_lists/all_cpm_in_genes_${type}.txt"
		done
	done
		
	cmd="$path_to_scripts/merge_many_files.sh $outdir/count_matrices/file_lists/all_cpm_in_genes_w_dup.txt $outdir/count_matrices/all_cpm_in_genes_w_dup.txt locus_name"
	bsub -o "$lsf/merge_many_files_all_cpm_in_genes_w_dup.txt" -K "$cmd" & pids1+=( $! )
	cmd="$path_to_scripts/merge_many_files.sh $outdir/count_matrices/file_lists/all_cpm_in_genes_dedup.txt $outdir/count_matrices/all_cpm_in_genes_dedup.txt locus_name"
	bsub -o "$lsf/merge_many_files_all_cpm_in_genes_dedup.txt" -K "$cmd" & pids2+=( $! )
	
	wait "$pids1" || { err_msg "error with merge_many_files, see $lsf/merge_many_files_all_cpm_in_genes_w_dup.txt" "$log"; }
	wait "$pids2" || { err_msg "error with merge_many_files, see $lsf/merge_many_files_all_cpm_in_genes_dedup.txt" "$log"; }
					
	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	

fi


# ----------------------
# Step 3: identify nuclei passing basic QC, filter out neg/pos ctrls and examine spike-ins
# ----------------------
if "$skip3"; then printf "\nSkipping step 3 by user request...\n" | tee -a "$log"
else printf "\nStep 3: identifying high-quality nuclei libraries and examining spike-ins\n" | tee -a "$log"; ts=$(date +%s)

	# store all nuclei that pass QC here; also remember how many passed, failed etc.
	passQC=()
	n_3N_pass=0; n_3N_fail=0; n_6N_pass=0; n_6N_fail=0; n_4N_pass=0; n_4N_fail=0
	n_negctrl_pass=0; n_negctrl_fail=0; n_2x3N_pass=0; n_2x3N_fail=0
	
	# nuclei pass QC if they detect at least $mingenescov genes and at least $mingenescov5 genes with 5 or more reads
	# (both are given as parameters above, mingenescov=1500 and mingenescov5=1000)
	echo "stubname	cross	dap	num_nuc	peak	stage	batch	n_genes_cov	n_genes_cov_min5	totreadsuniq	log10totreadsuniq" > "$summ/qcfilt_summary.txt"
	for ((i=0;i<${#stubname[@]};++i)); do
		[ "$i" = "0" ] && echo "Processing nuclei in batch ${sbatch[i]}..."
		[ "$i" != "0" ] && { [ "${sbatch[$i]}" != "${sbatch[$i-1]}" ] && echo "Processing nuclei in batch ${sbatch[i]}..."; }
		loc="$outdir/${project[i]}_map/${stubname[i]}_map"
		[ -f "$loc/htseq_count/${stubname[i]}_counts_w_dup_all_genes_cleaned.txt" ] && cuniq=$( awk -F$'\t' '{OFS=FS} $2 > 0 && $1 !~ /_/ {ss+=1} END {if (length(ss) == 0) {print 0} else {print ss}}' "$loc/htseq_count/${stubname[i]}_counts_w_dup_all_genes_cleaned.txt" ) || cuniq=0
		[ -f "$loc/htseq_count/${stubname[i]}_counts_w_dup_all_genes_cleaned.txt" ] && cuniq_5=$( awk -F$'\t' '{OFS=FS} $2 >= 5 && $1 !~ /_/ {ss+=1} END {if (length(ss) == 0) {print 0} else {print ss}}' "$loc/htseq_count/${stubname[i]}_counts_w_dup_all_genes_cleaned.txt" ) || cuniq_5=0
		
		totreadsuniq=$( awk '/mapped uniquely. Of these/ {print $1}' "$loc/${stubname[i]}_log.txt" )
		logtotreadsuniq=$( echo $totreadsuniq | awk '{printf "%11.4f\n",log($1)/log(10)}' )
						
		if [ "$cuniq" -gt $mingenescov ] && [ "$cuniq_5" -ge $mingenescov5 ]; then
			passQC[i]="Y"
			[ "${dap[i]}" = "nocell" ] && n_negctrl_pass=$(( n_negctrl_pass + 1 ))
			[ "${dap[i]}" = "blank" ] && n_negctrl_pass=$(( n_negctrl_pass + 1 ))
			[ "${peak[i]}" = "3N" ] && { [ "${num_nuc[i]}" = "1" ] && n_3N_pass=$(( n_3N_pass + 1 )); }
			[ "${peak[i]}" = "6N" ] && { [ "${num_nuc[i]}" = "1" ] && n_6N_pass=$(( n_6N_pass + 1 )); }
			[ "${peak[i]}" = "4N" ] && { [ "${num_nuc[i]}" = "1" ] && n_4N_pass=$(( n_4N_pass + 1 )); }
			[ "${peak[i]}" = "3N" ] && { [ "${num_nuc[i]}" = "2" ] && n_2x3N_pass=$(( n_2x3N_pass + 1 )); }
		else
			passQC[i]="N"
			[ "${dap[i]}" = "nocell" ] && n_negctrl_fail=$(( n_negctrl_fail + 1 ))
			[ "${dap[i]}" = "blank" ] && n_negctrl_fail=$(( n_negctrl_fail + 1 ))
			[ "${peak[i]}" = "3N" ] && { [ "${num_nuc[i]}" = "1" ] && n_3N_fail=$(( n_3N_fail + 1 )); }
			[ "${peak[i]}" = "6N" ] && { [ "${num_nuc[i]}" = "1" ] && n_6N_fail=$(( n_6N_fail + 1 )); }
			[ "${peak[i]}" = "4N" ] && { [ "${num_nuc[i]}" = "1" ] && n_4N_fail=$(( n_4N_fail + 1 )); }
			[ "${peak[i]}" = "3N" ] && { [ "${num_nuc[i]}" = "2" ] && n_2x3N_fail=$(( n_2x3N_fail + 1 )); }
		fi
		
		echo "${stubname[i]}	${cross[i]}	${dap[i]}	${num_nuc[i]}	${peak[i]}	${stage[i]}	${sbatch[i]}	$cuniq	$cuniq_5	$totreadsuniq	$logtotreadsuniq" >> "$summ/qcfilt_summary.txt"		
	done
		
	# summarize results
	rm -f "$summ/qcfilt_summary_counts.txt"
	echo "Summary of libraries passing QC:" | tee -a "$log"	"$summ/qcfilt_summary_counts.txt"
	echo " -> $n_3N_pass of $(( $n_3N_pass + $n_3N_fail )) single 3N nuclei passed basic QC" | tee -a "$log" "$summ/qcfilt_summary_counts.txt"
	echo " -> $n_6N_pass of $(( $n_6N_pass + $n_6N_fail )) single 6N nuclei passed basic QC" | tee -a "$log" "$summ/qcfilt_summary_counts.txt"
	echo " -> $n_4N_pass of $(( $n_4N_pass + $n_4N_fail )) single 4N nuclei (from 2 DAP) passed basic QC" | tee -a "$log" "$summ/qcfilt_summary_counts.txt"
	echo " -> $n_2x3N_pass of $(( $n_2x3N_pass + $n_2x3N_fail )) 2 x 3N pos. ctrls passed basic QC" | tee -a "$log" "$summ/qcfilt_summary_counts.txt"
	echo " -> $n_negctrl_pass of $(( $n_negctrl_pass + $n_negctrl_fail )) negative controls passed basic QC" | tee -a "$log" "$summ/qcfilt_summary_counts.txt"

	# output summary file of all libraries passing QC that are not from 2 nuclei
	echo "stubname	project	seqtype	cross	dap	peak	stage	batch	platform	ercc" > "$summ/singlenuc_passQC.txt"
	for ((i=0;i<${#stubname[@]};++i)); do
		[ "${passQC[i]}" = "Y" ] && [ "${num_nuc[i]}" = "1" ] && echo "${stubname[i]}	${project[i]}	${seqtype[i]}	${cross[i]}	${dap[i]}	${peak[i]}	${stage[i]}	${sbatch[i]}	${platform[i]}	${ercc[i]}" >> "$summ/singlenuc_passQC.txt"
	done
	
	# subset the large count and cpm matrices to keep only nuclei that passed QC
	mkdir -p "$outdir/count_matrices/all_nuclei"
	mv "$outdir/count_matrices/"*.txt "$outdir/count_matrices/all_nuclei"
	
	cut -f1 "$summ/singlenuc_passQC.txt" | tail -n+2 > "$summ/passQC_list.txt"
	
	mkdir -p "$outdir/count_matrices/logs"
	for allele in "all" "Col" "Cvi" "Ler"; do	
		for type in "w_dup" "dedup"; do
			$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_nuclei/${allele}_counts_in_genes_${type}.txt" "$outdir/count_matrices/${allele}_counts_in_genes_${type}_passQC.txt" --cols "$summ/passQC_list.txt" > "$outdir/count_matrices/logs/log_subset_${allele}_${type}_passQC.txt"
			[ $? != 0 ] && err_msg "error during subset_large_file, see $outdir/count_matrices/logs/log_subset_${allele}_${type}_passQC.txt" "$log"

			$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_nuclei/${allele}_counts_in_TEs_${type}.txt" "$outdir/count_matrices/${allele}_counts_in_TEs_${type}_passQC.txt" --cols "$summ/passQC_list.txt" > "$outdir/count_matrices/logs/log_subset_${allele}_${type}_passQC.txt"
			[ $? != 0 ] && err_msg "error during subset_large_file, see $outdir/count_matrices/logs/log_subset_${allele}_${type}_passQC.txt" "$log"

			$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_nuclei/all_cpm_in_genes_${type}.txt" "$outdir/count_matrices/all_nuclei/all_cpm_in_genes_${type}_passQC.txt" --cols "$summ/passQC_list.txt" > "$outdir/count_matrices/logs/log_subset_all_cpm_${type}_passQC.txt"
			[ $? != 0 ] && err_msg "error during subset_large_file, see $outdir/count_matrices/logs/log_subset_${allele}_${type}_passQC.txt" "$log"

			# in TE and intergenic matrices, drop all rows where no nuclei had counts	
			awk -F$'\t' '{if(NR==1) {print $0} else { for(i=1; i<=NF;i++) j+=$i; {if(j > 0) {print $0}}; j=0 }}' "$outdir/count_matrices/${allele}_counts_in_TEs_${type}_passQC.txt" > "$outdir/count_matrices/${allele}_counts_in_TEs_${type}_passQC_nozero.txt"

			cat "$outdir/count_matrices/${allele}_counts_in_genes_${type}_passQC.txt" "$outdir/count_matrices/${allele}_counts_in_TEs_${type}_passQC_nozero.txt" "$outdir/count_matrices/${allele}_counts_in_intergenic_${type}_passQC_nozero.txt" > "$outdir/count_matrices/${allele}_counts_in_all_${type}_passQC.txt"
		done
	done	
fi


# read-in stubname, etc. of libraries that passed QC (and are from just 1 nucleus); proceed with these from now on;
i=0
stubname=()
project=()
seqtype=()
cross=()
dap=()
peak=()
stage=()
sbatch=()
strain2=()

while IFS=$'\t' read -r -a aa; do
	stubname[i]="${aa[0]}"
	project[i]="${aa[1]}"
	seqtype[i]="${aa[2]}"
	cross[i]="${aa[3]}"
	dap[i]="${aa[4]}"
	peak[i]="${aa[5]}"
	stage[i]="${aa[6]}"
	sbatch[i]="${aa[7]}"
	[[ ${cross[i]} == *"L"* ]] && strain2[i]="Ler"
	[[ ${cross[i]} == *"V"* ]] && strain2[i]="Cvi"
	i=$(( $i + 1 ))
done < <( tail -n+2 "$summ/singlenuc_passQC.txt" )


# ----------------------
# Step 4: clustering by tSNE and SC3, assign tissue of origin based on clustering
# ----------------------
if "$skip4"; then printf "\nSkipping step 4 by user request...\n" | tee -a "$log"
else printf "\nStep 4: initial clustering\n" | tee -a "$log"; ts=$(date +%s)
			
	# PCA and preliminary t-SNE analysis
	# -------------------
	mkdir -p "$outdir/clustering/PCA/passQC_genes"
	mkdir -p "$outdir/clustering/PCA/passQC_TEs"
	mkdir -p "$outdir/clustering/PCA/passQC_all"
	
	# calculate total % maternal for each nucleus (using SNP-assigned reads) and add to summary stats file
	cut -f1,2,19,20 "$summ/rna_seq_map_summary.txt" | awk -F$'\t' '{OFS=FS} {if($3+$4==0){fCol=0} else {fCol = $3/($3+$4)}; if ($2 ~ /Cx/) {fMat=fCol} else {fMat=(1-fCol)}; if (NR==1){print $1,"fMat","sortorder"} else {print $1,fMat,NR}}' > "$outdir/clustering/overall_pMat.txt"
	$path_to_scripts/merge_by_column.R "$summ/singlenuc_passQC.txt" "$outdir/clustering/overall_pMat.txt" stubname "-" --tokeep allx | sort -k12n,12n -t$'\t' | cut -f1-11 > "$summ/singlenuc_passQC_wPmat.txt"
			
	# run tSNE
	$path_to_scripts/single_cell_cluster_PCA_tSNE.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$summ/singlenuc_passQC_wPmat.txt" "$outdir/clustering/PCA/passQC_genes/passQC_genes_k3" --k 3 --perplexity 50 --tSNE_PCs 5 > "$outdir/clustering/PCA/passQC_genes_log.txt"
	$path_to_scripts/single_cell_cluster_PCA_tSNE.R "$outdir/count_matrices/all_counts_in_TEs_dedup_passQC_nozero.txt" "$summ/singlenuc_passQC_wPmat.txt" "$outdir/clustering/PCA/passQC_TEs/passQC_TEs_k3" --k 3 --perplexity 10 --minexpr 100 --tSNE_PCs 3 > "$outdir/clustering/PCA/passQC_TEs_log.txt"	
	
	# based on tSNE result, identify seed coat genes (cluster 2 with default seed)
	tissue=(); i=0	
	while IFS=$'' read -r line; do
		tissue[i]=$( echo "$line" | awk -F$'\t' '{OFS=FS} {if($12 == 2 || $11 > 0.85) {print "seedcoat"} else if($11<0.6) {print "embryo"} else {print "endo"}}' )
		i=$(( $i + 1 ))
	done < <( tail -n+2 "$outdir/clustering/PCA/passQC_genes/passQC_genes_k3_tSNE_perp50_5PCs_clusters.txt" )
	
	echo "stubname	tissue	1" > "$outdir/_summary_stats/tissue_based_on_clus.txt"
	for ((i=0;i<${#tissue[@]};++i)); do
		echo "${stubname[i]}	${tissue[i]}	$(( $i + 1 ))" >> "$outdir/_summary_stats/tissue_based_on_clus.txt"
	done
		
	# add preliminary tissue assignments to summary stats file
	$path_to_scripts/merge_by_column.R "$summ/singlenuc_passQC_wPmat.txt" "$outdir/_summary_stats/tissue_based_on_clus.txt" stubname "-" --tokeep allx | sort -k13n,13n -t$'\t' | cut -f1-12 > "$summ/singlenuc_passQC_wTissue.txt"


	# SC3 clustering
	# -------------------
	# run SC3 over all nuclei that passed QC
	mkdir -p "$outdir/clustering/SC3_all"
	
	$path_to_scripts/single_cell_cluster_SC3.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$summ/singlenuc_passQC_wTissue.txt" "$outdir/clustering/SC3_all/SC3_all" > "$outdir/clustering/SC3_all/SC3_all_log.txt"

	# add info on a small variant in snRNA-seq prep protocol (wash vs. no wash, see methods) to summary stats file
	awk -F$'\t' '{OFS=FS} {if ($2 == "project") {print $0,"washed"} else if ($2 == "180928Geh") {print $0,"yes"} else {print $0,"no"}}' "$outdir/_summary_stats/singlenuc_passQC_wTissue.txt" > "$outdir/_summary_stats/singlenuc_passQC_wWash.txt"

	# remake SC3 plots using custom script and the main factors of interest (tissue, genotype, ploidy peak, DAP, % mat, sequencing type, wash)
	cut -f1,4,5,6,12,13 "$outdir/_summary_stats/singlenuc_passQC_wWash.txt" > "$outdir/clustering/SC3_all/newfactors.txt"	
	ff="$outdir/clustering/SC3_all/SC3_all_23_consensus_matrix.txt"
	factorsfile="$outdir/clustering/SC3_all/newfactors.txt"
	bb=${ff%"_consensus_matrix.txt"}
	bbd=$( dirname "$bb" )		
	consensus="${bb}_consensus_matrix.txt"
	clusfile="${bb}_cluster_info.txt"	
	$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig" --metadata "$factorsfile"
	$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig_nocluscolors" --metadata "$factorsfile" --nocluscolors

	# try to subcluster all clusters (some will fail b/c of too few nuclei, etc.)
	pids=()
	for ((k=1;k<=23;++k)); do
		echo "Trying SC3 on cluster $k"
		loc="$outdir/clustering/SC3_all_subclus/SC3_all_subclus${k}"; mkdir -p "$loc"
		awk -F$'\t' -v clusID="$k" '$2 == clusID' "$outdir/clustering/SC3_all/SC3_all_23_cluster_info.txt" | cut -f1 > "$loc/orig_cluster_nuclei.txt"

		$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$loc/orig_cluster_counts.txt" --cols "$loc/orig_cluster_nuclei.txt"		
		$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_wTissue.txt" "$loc/orig_cluster_nucdesc.txt" --rows "$loc/orig_cluster_nuclei.txt"

		cmd="$path_to_scripts/single_cell_cluster_SC3.R $loc/orig_cluster_counts.txt $loc/orig_cluster_nucdesc.txt $loc/SC3_all_subclus${k} --maxK 10 --mincells 0 --minreads 5 > $loc/SC3_all_subclus${k}_log.txt 2>&1"
		bsub -o "$loc/LSF_log.txt" -K "$cmd" & pids[k]=$!
	done
		
	for ((k=1;k<=23;++k)); do
		wait ${pids[k]} || echo "SC3 failed when subclustering primary cluster ${k}, see $loc/LSF_log.txt; script will continue anyway"
	done

	# the only subcluster we'll use further is cluster 21, which is really not fully clustered and contains some smaller tight clusters
	# - 21 has 5 subclusters - these will be referred to as 21_1,...,21_5
	awk -F$'\t' '{OFS=FS} NR==FNR{a[$1] = $2; next} {if ($1 in a) {print $1,$2"_"a[$1]} else {print $0}}' "$outdir/clustering/SC3_all_subclus/SC3_all_subclus21/SC3_all_subclus21_5_cluster_info.txt" "$outdir/clustering/SC3_all/SC3_all_23_cluster_info.txt" > "$outdir/clustering/SC3_all_final_clusters.txt"

	# remake plot for subclustered cluster 21 too with custom script and factors of interest
	k=21
	loc="$outdir/clustering/SC3_all_subclus/SC3_all_subclus${k}"
	ff="$loc/SC3_all_subclus21_5_consensus_matrix.txt"
	$path_to_scripts/subset_large_file.R "$outdir/clustering/SC3_all/newfactors.txt" "$loc/newfactors.txt" --rows "$loc/orig_cluster_nuclei.txt"
	factorsfile="$loc/newfactors.txt"
	bb=${ff%"_consensus_matrix.txt"}
	bbd=$( dirname "$bb" )
		
	consensus="${bb}_consensus_matrix.txt"
	clusfile="${bb}_cluster_info.txt"
	
	$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig" --metadata "$factorsfile"
	$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig_nocluscolors" --metadata "$factorsfile" --nocluscolors


	# use SC3 results to refine tissue assignments
	# -------------------
	
	# two nuclei cluster with embryo but are listed as seed coat (based on tSNE) with slightly more than 60% maternal reads;
	# re-assign to embryo
	awk -F$'\t' '$25 == 5 && $11 < 0.7 {print $1}' "$outdir/clustering/SC3_all_subclus/SC3_all_subclus21/SC3_all_subclus21_all_cluster_info.txt" > "$outdir/clustering/embryos.txt"
	embryos=$( cat $outdir/clustering/embryos.txt | tr '\n' ' ' )

	# for nearly homogeneous seed coat clusters (5,6,7,9,16,18,20,21_1,21_3) - assign all members of these clusters to seed coat
	awk -F$'\t' '$2 == 5 || $2 == 6 || $2 == 7 || $2 == 9 || $2 == 16 || $2 == 18 || $2 == 20 || $2 == 21_1 || $2 == 21_3' "$outdir/clustering/SC3_all_final_clusters.txt" > "$outdir/clustering/seedcoat.txt"
	seedcoat=$( cat $outdir/clustering/seedcoat.txt | tr '\n' ' ' )

	# for nearly homogeneous endosperm clusters (1,2,3,4,8,10,11,12,13,14,15,17,19,21_2,21_4,22,23) - assign all members of these clusters to endosperm
	awk -F$'\t' '$2 == 1 || $2 == 2 || $2 == 3 || $2 == 4 || $2 == 8 || $2 == 10 || $2 == 11 || $2 == 12 || $2 == 13 || $2 == 14 || $2 == 15 || $2 == 17 || $2 == 19 || $2 == 21_2 || $2 == 21_4 || $2 == 22 || $2 == 23' "$outdir/clustering/SC3_all_final_clusters.txt" > "$outdir/clustering/endosperm.txt"
	endosperm=$( cat $outdir/clustering/endosperm.txt | tr '\n' ' ' )

	echo "stubname	project	seqtype	cross	dap	peak	stage	batch	fmat	tissue	washed" > "$outdir/_summary_stats/singlenuc_passQC_final.txt"
	
	washed=()
	for ((i=0;i<${#stubname[@]};++i)); do
		for ss in $embryos; do
			if [ "${stubname[i]}" = "$ss" ]; then
				tissue[i]="embryo"
			fi
		done
		for ss in $endosperm; do
			if [ "${stubname[i]}" = "$ss" ]; then
				tissue[i]="endo"
			fi
		done
		for ss in $seedcoat; do
			if [ "${stubname[i]}" = "$ss" ]; then
				tissue[i]="seedcoat"
			fi
		done
		[ "${project[i]}" = "180928Geh" ] && washed[i]="yes" || washed[i]="no"
		echo "${stubname[i]}	${project[i]}	${seqtype[i]}	${cross[i]}	${dap[i]}	${peak[i]}	${stage[i]}	${sbatch}	${fmat[i]}	${tissue[i]}	${washed[i]}" >> "$outdir/_summary_stats/singlenuc_passQC_final.txt"
	done
	
	
	i=0
	stubname=()
	project=()
	seqtype=()
	cross=()
	dap=()
	peak=()
	stage=()
	sbatch=()
	strain2=()
	fmat=()
	tissue=()
	wash=()

	while IFS=$'\t' read -r -a aa; do
		stubname[i]="${aa[0]}"
		project[i]="${aa[1]}"
		seqtype[i]="${aa[2]}"
		cross[i]="${aa[3]}"
		dap[i]="${aa[4]}"
		peak[i]="${aa[5]}"
		stage[i]="${aa[6]}"
		sbatch[i]="${aa[7]}"
		fmat[i]="${aa[8]}"
		tissue[i]="${aa[9]}"
		wash[i]="${aa[10]}"
		[[ ${cross[i]} == *"L"* ]] && strain2[i]="Ler"
		[[ ${cross[i]} == *"V"* ]] && strain2[i]="Cvi"
		i=$(( $i + 1 ))
	done < <( tail -n+2 "$summ/singlenuc_passQC_final.txt" )

	
	# regenerate SC3_all plots using custom script and updated tissue info
	cut -f1,4,5,6,10,11 "$outdir/_summary_stats/singlenuc_passQC_final.txt" > "$outdir/clustering/SC3_all/newfactors.txt"	
	ff="$outdir/clustering/SC3_all/SC3_all_23_consensus_matrix.txt"
	factorsfile="$outdir/clustering/SC3_all/newfactors.txt"
	bb=${ff%"_consensus_matrix.txt"}
	bbd=$( dirname "$bb" )	
	consensus="${bb}_consensus_matrix.txt"
	clusfile="${bb}_cluster_info.txt"	
	$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus" --metadata "$factorsfile"
	$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_nocluscolors" --metadata "$factorsfile" --nocluscolors
		
	# also make plot with single factor per plot
	numfactors=$( awk -F$'\t' '{print NF-1}' "$factorsfile" | head -1 )
	for ((i=1;i<=${numfactors};++i)); do
		fac=$(( $i + 1 ))
		cut -f1,$fac "$factorsfile" > "$bbd/tmp.txt"
		factor=$( head -1 "$bbd/tmp.txt" | cut -f2 )
		$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_${factor}" --metadata "$bbd/tmp.txt" --nocluscolors
	done
	
	# repeat w/ cluster 21 subset only
	awk -F$'\t' '$2 == 21 {print $1}' "$outdir/clustering/SC3_all/SC3_all_23_cluster_info.txt" > "$outdir/clustering/SC3_all/SC3_all_list_of_clus21_nuclei.txt"
	awk -F$'\t' 'NR==1 || $2 == 21' "$outdir/clustering/SC3_all/SC3_all_23_cluster_info.txt" > "$outdir/clustering/SC3_all/SC3_all_list_of_clus21_nuclei_wclus.txt"
	$path_to_scripts/subset_large_file.R "$outdir/clustering/SC3_all/SC3_all_23_consensus_matrix.txt" "$outdir/clustering/SC3_all/SC3_all_23_consensus_matrix_clus21only.txt" --cols "$outdir/clustering/SC3_all/SC3_all_list_of_clus21_nuclei.txt" --rows "$outdir/clustering/SC3_all/SC3_all_list_of_clus21_nuclei.txt"
	$path_to_scripts/subset_large_file.R "$outdir/clustering/SC3_all/newfactors.txt" "$outdir/clustering/SC3_all/newfactors_clus21only.txt" --rows "$outdir/clustering/SC3_all/SC3_all_list_of_clus21_nuclei.txt"	
	$path_to_scripts/plot_SC3_heatmap.py "$outdir/clustering/SC3_all/SC3_all_23_consensus_matrix_clus21only.txt" "$outdir/clustering/SC3_all/SC3_all_list_of_clus21_nuclei_wclus.txt" "$outdir/clustering/SC3_all/SC3_all_clus21only" --metadata "$outdir/clustering/SC3_all/newfactors_clus21only.txt" --nocluscolors
	
	
	# repeat SC3 analysis on all nuclei subsets of interest (e.g. CxV endo, VxC endo, CxV 4DAP endo, VxC 4DAP endo, CxL endo, LxC endo, and all of those SC.)
	# in paper we focus primarily on CxV 4 DAP endo, VxC 4 DAP endo, CxV 4 DAP seedcoat and VxC 4 DAP seedcoat
	# -------------------
	rm -f $outdir/clustering/SC3_subsets/*.txt
	mkdir -p "$outdir/clustering/SC3_subsets"
	for ((i=0;i<${#stubname[@]};++i)); do
		if [ "${tissue[i]}" = "endo" ]; then
			[ "${cross[i]}" = "CxV" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/CxV_endo.txt"
			[ "${cross[i]}" = "VxC" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/VxC_endo.txt"
			[ "${cross[i]}" = "CxL" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/CxL_endo.txt"
			[ "${cross[i]}" = "LxC" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/LxC_endo.txt"

			[ "${cross[i]}" = "CxV" ] && [ "${dap[i]}" = "4" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/CxV_endo_4DAP.txt"
			[ "${cross[i]}" = "VxC" ] && [ "${dap[i]}" = "4" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/VxC_endo_4DAP.txt"
		fi
		if [ "${tissue[i]}" = "seedcoat" ]; then
			[ "${cross[i]}" = "CxV" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/CxV_seedcoat.txt"
			[ "${cross[i]}" = "VxC" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/VxC_seedcoat.txt"
			[ "${cross[i]}" = "CxL" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/CxL_seedcoat.txt"
			[ "${cross[i]}" = "LxC" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/LxC_seedcoat.txt"

			[ "${cross[i]}" = "CxV" ] && [ "${dap[i]}" = "4" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP.txt"
			[ "${cross[i]}" = "VxC" ] && [ "${dap[i]}" = "4" ] && echo "${stubname[i]}" >> "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP.txt"
		fi
	done

	# also run CxV + VxC 4 DAP endo and seedcoat
	cat "$outdir/clustering/SC3_subsets/CxV_endo_4DAP.txt" "$outdir/clustering/SC3_subsets/VxC_endo_4DAP.txt" > "$outdir/clustering/SC3_subsets/CandV_endo_4DAP.txt"
	cat "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP.txt" "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP.txt" > "$outdir/clustering/SC3_subsets/CandV_seedcoat_4DAP.txt"
		
	# for each of these subsets, get subset from main counts file and run SC3
	pids=()
	mv "$outdir/clustering/SC3_subsets/overall_pMat.txt" "$summ"
	for ff in $outdir/clustering/SC3_subsets/*.txt; do
		echo "Processing $ff"
		bb=$(basename $ff .txt)
		mkdir -p "$outdir/clustering/SC3_subsets/$bb"
		
		$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$outdir/clustering/SC3_subsets/$bb/${bb}_counts.txt" --cols "$ff"		
		$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_final.txt" "$outdir/clustering/SC3_subsets/$bb/${bb}_descfile.txt" --rows "$ff"
		
		cmd="$path_to_scripts/single_cell_cluster_SC3.R $outdir/clustering/SC3_subsets/$bb/${bb}_counts.txt $outdir/clustering/SC3_subsets/$bb/${bb}_descfile.txt $outdir/clustering/SC3_subsets/$bb/${bb} --maxK 15 > $outdir/clustering/SC3_subsets/$bb/${bb}_log.txt 2>&1"
		bsub -o "$outdir/clustering/SC3_subsets/$bb/${bb}_LSF_log.txt" -K "$cmd" & pids+=( $! )
	done

	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "SC3 failed, see $outdir/clustering logs" "$log"
	done
		
	# for CxV 4DAP endo, re-cluster clusters 10:
	awk -F$'\t' '{OFS=FS} NR!=1 && $2 == 10 {print $1}' "$outdir/clustering/SC3_subsets/CxV_endo_4DAP/CxV_endo_4DAP_10_cluster_info.txt" > "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo.txt"
	mkdir -p "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo"
	
	$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_counts.txt" --cols "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo.txt"		
	$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_final.txt" "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_descfile.txt" --rows "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo.txt"
	$path_to_scripts/single_cell_cluster_SC3.R $outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_counts.txt $outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_descfile.txt $outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo --maxK 15 > $outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_log.txt 2>&1

	# for VxC 4DAP endo, re-cluster clusters 6-8:
	awk -F$'\t' '{OFS=FS} NR!=1 && $2 >= 6 {print $1}' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP/VxC_endo_4DAP_8_cluster_info.txt" > "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo.txt"
	mkdir -p "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo"
	
	$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo/VxC_endo_4DAP_6to8_redo_counts.txt" --cols "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo.txt"		
	$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_final.txt" "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo/VxC_endo_4DAP_6to8_redo_descfile.txt" --rows "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo.txt"
	$path_to_scripts/single_cell_cluster_SC3.R $outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo/VxC_endo_4DAP_6to8_redo_counts.txt $outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo/VxC_endo_4DAP_6to8_redo_descfile.txt $outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo/VxC_endo_4DAP_6to8_redo --maxK 15 > $outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo/VxC_endo_4DAP_6to8_redo_log.txt 2>&1

	# for VxC 4DAP seedcoat, re-cluster cluster 7:
	awk -F$'\t' '{OFS=FS} NR!=1 && $2 == 7 {print $1}' "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP/VxC_seedcoat_4DAP_7_cluster_info.txt" > "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo.txt"
	mkdir -p "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo"
	
	$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo/VxC_seedcoat_4DAP_7_redo_counts.txt" --cols "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo.txt"		
	$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_final.txt" "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo/VxC_seedcoat_4DAP_7_redo_descfile.txt" --rows "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo.txt"
	$path_to_scripts/single_cell_cluster_SC3.R $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo/VxC_seedcoat_4DAP_7_redo_counts.txt $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo/VxC_seedcoat_4DAP_7_redo_descfile.txt $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo/VxC_seedcoat_4DAP_7_redo --maxK 15 > $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo/VxC_seedcoat_4DAP_7_redo_log.txt 2>&1

	# note - tried seeing if different normalization approaches (e.g. DEnorm approach from DEsingle) affected
	# clustering; effects were very minor, so proceeded w/ CPM normalization as it's simpler

	# re-generate plots for all consensus plots in clustering/SC3_subsets using custom script
	for ff in $outdir/clustering/SC3_subsets/*/*_consensus_matrix.txt; do
		bb=${ff%"_consensus_matrix.txt"}
		bbd=$( dirname "$bb" )
		cut -f1,4,5,6,10,11 $( ls "$bbd"/*_descfile.txt ) > $bbd/factors_for_plot.txt
		todrop=""

		# filter out factors with no variation in this group
		for ((i=2;i<=6;++i)); do
			numvals=$( cut -f${i} "$bbd/factors_for_plot.txt" | sort -u | wc -l )
			if [ ${numvals} -le 2 ]; then
				todrop="${todrop},${i}"
			fi
		done
		todrop2=$( echo "$todrop" | cut -c 2- )
		cut -f${todrop2} --complement "$bbd/factors_for_plot.txt" > "$bbd/tmp.txt"
		mv "$bbd/tmp.txt" "$bbd/factors_for_plot.txt"
		factorsfile="$bbd/factors_for_plot.txt"
				
		consensus="${bb}_consensus_matrix.txt"
		clusfile="${bb}_cluster_info.txt"
	
		$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus" --metadata "$factorsfile"
		$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_nocluscolors" --metadata "$factorsfile" --nocluscolors
		
		# also make plot with single factor per plot
		numfactors=$( awk -F$'\t' '{print NF-1}' "$factorsfile" | head -1 )
		for ((i=1;i<=${numfactors};++i)); do
			fac=$(( $i + 1 ))
			cut -f1,$fac "$factorsfile" > "$bbd/tmp.txt"
			factor=$( head -1 "$bbd/tmp.txt" | cut -f2 )
			$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_${factor}" --metadata "$bbd/tmp.txt" --nocluscolors
		done
	done	
	
	# what do the SC3_all (Fig. S4) clusters look like in the per-genotype clusters (Fig. 1B, Fig. S5A)?	
	echo "stubname	SC3_all" > "$outdir/clustering/SC3_all_final_clusters_fxheader.txt"
	tail -n+2 "$outdir/clustering/SC3_all_final_clusters.txt" >> "$outdir/clustering/SC3_all_final_clusters_fxheader.txt"
	
	for type in CxV_endo_4DAP CxV_endo_4DAP_10_redo VxC_endo_4DAP VxC_endo_4DAP_6to8_redo CxV_seedcoat_4DAP VxC_seedcoat_4DAP VxC_seedcoat_4DAP_7_redo; do
		$path_to_scripts/merge_by_column.R "$outdir/clustering/SC3_subsets/$type/factors_for_plot.txt" "$outdir/clustering/SC3_all_final_clusters_fxheader.txt" stubname "-" --tokeep allx | cut -f1,4 > "$outdir/clustering/SC3_subsets/$type/factors_for_plot2.txt"
		
		consensus=$( ls $outdir/clustering/SC3_subsets/$type/${type}*_consensus_matrix.txt )
		clusfile=$( ls $outdir/clustering/SC3_subsets/$type/${type}*_cluster_info.txt | head -1 )
		$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "$outdir/clustering/SC3_subsets/$type/${type}_consensus_SC3all" --metadata "$outdir/clustering/SC3_subsets/$type/factors_for_plot2.txt" --nocluscolors
	done
	
	# based on visual inspection, further refine clustering
	# CxV 4DAP endo cluster 12 has two clear subclusters but is too small to re-cluster with SC3, so separate
	# those manually; saved to "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/manual_clus_4_list.txt"
	echo "Original:"
	cat "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_4_cluster_info.txt"
	
	echo "New:"
	awk -F$'\t' '{OFS=FS} NR==FNR {a[$1]=1; next} {if ($1 in a) {print $1,"4"} else if ($2 == 4) {print $1,"5"} else {print $1,$2}}' "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/manual_clus_4_list.txt" "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_4_cluster_info.txt" > "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_final_clusters.txt"
	cat "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_final_clusters.txt"
	
	# remake plots 
	factorsfile="$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/factors_for_plot.txt"
	consensus="$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_4_consensus_matrix.txt"
	clusfile="$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_final_clusters.txt"
	numfactors=$( awk -F$'\t' '{print NF-1}' "$factorsfile" | head -1 )
	for ((i=1;i<=${numfactors};++i)); do
		fac=$(( $i + 1 ))
		cut -f1,$fac "$factorsfile" > "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/tmp.txt"
		factor=$( head -1 "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/tmp.txt" | cut -f2 )
		$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_final_consensus_${factor}" --metadata "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/tmp.txt" --nocluscolors
	done

	type="CxV_endo_4DAP_10_redo"
	$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "$outdir/clustering/SC3_subsets/$type/${type}_consensus_SC3all_final" --metadata "$outdir/clustering/SC3_subsets/$type/factors_for_plot2.txt" --nocluscolors
	
	# VxC endo 4DAP cluster 3 also has additional subsets inside, based on SC3_all clustering; separate those now
	awk -F$'\t' '{OFS=FS} NR!=1 && $2 == 3 {print $1}' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP/VxC_endo_4DAP_8_cluster_info.txt" > "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo.txt"
	mkdir -p "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo"
	
	$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_counts.txt" --cols "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo.txt"		
	$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_final.txt" "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_descfile.txt" --rows "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo.txt"
	$path_to_scripts/single_cell_cluster_SC3.R $outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_counts.txt $outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_descfile.txt $outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo --maxK 6 > $outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_log.txt 2>&1
		
	# remake plots
	cut -f1,6,11 "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_descfile.txt" > $outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/factors_for_plot.txt 
	factorsfile="$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/factors_for_plot.txt"
	consensus="$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_2_consensus_matrix.txt"
	clusfile="$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_2_cluster_info.txt"
	numfactors=$( awk -F$'\t' '{print NF-1}' "$factorsfile" | head -1 )
	for ((i=1;i<=${numfactors};++i)); do
		fac=$(( $i + 1 ))
		cut -f1,$fac "$factorsfile" > "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/tmp.txt"
		factor=$( head -1 "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/tmp.txt" | cut -f2 )
		$path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_consensus_${factor}" --metadata "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/tmp.txt" --nocluscolors
	done
	

	# make final cluster lists for CxV, VxC 4DAP SC and endo
	# -------------------
	
	# CxV 4DAP endo
	awk '$2 != 10' "$outdir/clustering/SC3_subsets/CxV_endo_4DAP/CxV_endo_4DAP_10_cluster_info.txt" > "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt"
	awk -F$'\t' '{OFS=FS} NR != 1 {print $1,$2+9}' "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_10_redo/CxV_endo_4DAP_10_redo_final_clusters.txt" >> "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt"
	
	# CxV 4DAP SC
	cat "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP/CxV_seedcoat_4DAP_6_cluster_info.txt" > "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt"

	# VxC 4DAP endo	
	awk '$2 <= 2 || NR == 1' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP/VxC_endo_4DAP_8_cluster_info.txt" > "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt"
	awk -F$'\t' '{OFS=FS} NR != 1 {print $1,$2+2}' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_3_redo/VxC_endo_4DAP_3_redo_2_cluster_info.txt" >> "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt"
	awk -F$'\t' '{OFS=FS} ($2 > 3 && $2 <= 5) {print $1,$2+1}' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP/VxC_endo_4DAP_8_cluster_info.txt" >> "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt"
	awk -F$'\t' '{OFS=FS} NR != 1 {print $1,$2+6}' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_6to8_redo/VxC_endo_4DAP_6to8_redo_5_cluster_info.txt" >> "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt"
	
	# VxC 4DAP SC	
	awk '$2 <= 6 || NR == 1' "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP/VxC_seedcoat_4DAP_7_cluster_info.txt" > "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt"
	awk -F$'\t' '{OFS=FS} NR != 1 {print $1,$2+6}' "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_7_redo/VxC_seedcoat_4DAP_7_redo_2_cluster_info.txt" >> "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt"
	
	# make versions of each file w/ no header, for later
	tail -n+2 $outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt > $outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt
	tail -n+2 $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt > $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt
	tail -n+2 $outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt > $outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt
	tail -n+2 $outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt > $outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt
		
	# for nuclei in each cluster, make BAM files of all reads in that cluster
	rm -rf "$outdir/per_cluster_BAMs"
	for type in "CxV_endo_4DAP" "CxV_seedcoat_4DAP" "VxC_endo_4DAP" "VxC_seedcoat_4DAP"; do
		mkdir -p "$outdir/per_cluster_BAMs/$type"
	done
	
	for ((i=0;i<${#stubname[@]};++i)); do
		if grep -q "${stubname[i]}" $outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt; then
			clus=$( grep "${stubname[i]}" $outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt | cut -f2 )
			echo "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_all_dedup.bam" >> "$outdir/per_cluster_BAMs/CxV_endo_4DAP/CxV_endo_4DAP_cluster_${clus}_filelist.txt"
		elif grep -q "${stubname[i]}" $outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt; then
			clus=$( grep "${stubname[i]}" $outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt | cut -f2 )
			echo "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_all_dedup.bam" >> "$outdir/per_cluster_BAMs/CxV_seedcoat_4DAP/CxV_seedcoat_4DAP_cluster_${clus}_filelist.txt"
		elif grep -q "${stubname[i]}" $outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt; then 
			clus=$( grep "${stubname[i]}" $outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt | cut -f2 )
			echo "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_all_dedup.bam" >> "$outdir/per_cluster_BAMs/VxC_endo_4DAP/VxC_endo_4DAP_cluster_${clus}_filelist.txt"
		elif grep -q "${stubname[i]}" $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt; then
			clus=$( grep "${stubname[i]}" $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt | cut -f2 )
			echo "$outdir/${project[i]}_map/${stubname[i]}_map/assign_to_allele/${stubname[i]}_all_dedup.bam" >> "$outdir/per_cluster_BAMs/VxC_seedcoat_4DAP/VxC_seedcoat_4DAP_cluster_${clus}_filelist.txt"
		fi
	done
	
	# merge BAM files into one file per cluster (just for browsing)
	pids=()
	mkdir -p "$lsf/per_cluster_BAMs"
	for ff in $outdir/per_cluster_BAMs/*/*_filelist.txt; do
		bb=${ff%"_filelist.txt"}
		bbd=$( basename "$bb" )
		bbdd=${bbd%"_filelist.txt"}
		tomerge=$( cat "$ff" | tr '\n' ' ' )
		
		cmd="samtools merge -f ${bb}.bam $tomerge"
		bsub -o "$lsf/per_cluster_BAMs/${bbdd}.txt" -K "$cmd" & pids+=( $! )
	done

	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "samtools merge failed, see log files in $lsf/per_cluster_BAMs" "$log"
	done
	
	# index all BAM files
	for ff in $outdir/per_cluster_BAMs/*/*.bam; do
		samtools index "$ff"
	done

fi

i=0
stubname=()
project=()
seqtype=()
cross=()
dap=()
peak=()
stage=()
sbatch=()
strain2=()
fmat=()
tissue=()
wash=()

while IFS=$'\t' read -r -a aa; do
	stubname[i]="${aa[0]}"
	project[i]="${aa[1]}"
	seqtype[i]="${aa[2]}"
	cross[i]="${aa[3]}"
	dap[i]="${aa[4]}"
	peak[i]="${aa[5]}"
	stage[i]="${aa[6]}"
	sbatch[i]="${aa[7]}"
	fmat[i]="${aa[8]}"
	tissue[i]="${aa[9]}"
	wash[i]="${aa[10]}"
	[[ ${cross[i]} == *"L"* ]] && strain2[i]="Ler"
	[[ ${cross[i]} == *"V"* ]] && strain2[i]="Cvi"
	i=$(( $i + 1 ))
done < <( tail -n+2 "$summ/singlenuc_passQC_final.txt" )


# ----------------------
# Step 5: cell cycle analysis
# ----------------------
if "$skip5"; then printf "\nSkipping step 5 by user request...\n" | tee -a "$log"
else printf "\nStep 5: cell cycle analysis\n" | tee -a "$log"; ts=$(date +%s)

	# run trajectory analysis to assign nuclei to cell cycle phases & identify genes varying significantly w/ cell cycle (Ext. data Fig. 7)
	# -------------------

	# Ben showed that t-SNE over cell cycle markers recovered a nice circle-ish pattern
	# repeating his analysis here with same params to confirm replication and adding some
	# custom trajectory analysis stuff
	mkdir -p "$outdir/cell_cycle"
	matrixfile="$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt"
	cellinfo="$outdir/_summary_stats/singlenuc_passQC_final.txt"
	outprefix="$outdir/cell_cycle/traj_analysis"
	genelist="/lab/solexa_gehring/colette/single_cell_seq/cell_cycle_markers_Ben.txt"

	# do trajectory analysis with final params
	$path_to_scripts/single_cell_trajectory_analysis.R "$matrixfile" "$outprefix" --genelist "$genelist" --seed 1234567 --CPM --log2 --minreads 10 --mincells 5 --minexpr 1000 --perp 100 --PCs 2 --k 6 --smoother 'dijkstra' --fracdraw 0.5 --numiter 200 --neighbors 10 --sigthreshold 0.01 --kgenes 12 --cell_info "$cellinfo"

	# assign each k-means nuclei cluster to a cell cycle phase (note with fixed seed this
	# clustering is fully reproducible, but will change if you change the random seed --seed parameter in
	# command above)
	awk -F$'\t' '{OFS=FS} { if (NR == 1) {
		print "locus_name","cc_phase","cc_trajectory"
	} else {
		if ($4 == 1) { print $1,"G2",$5 }
		else if ($4 == 2) { print $1,"G1toS",$5 }
		else if ($4 == 3) { print $1,"S",$5 }
		else if ($4 == 4) { print $1,"M",$5 }
		else if ($4 == 5) { print $1,"G1",$5 }
		else if ($4 == 6) { print $1,"G0",$5 }
	}}' "$outdir/cell_cycle/traj_analysis_trajectory.txt" > "$outdir/cell_cycle/traj_analysis_cc_assignments.txt"

	# merge back to full nuclei stats table
	$path_to_scripts/merge_by_column.R <( sed '1s/stubname/locus_name/g' "$outdir/_summary_stats/singlenuc_passQC_final.txt" ) "$outdir/cell_cycle/traj_analysis_cc_assignments.txt" locus_name "-" | awk -F$'\t' '{OFS=FS} $12 == "" {$12 = "unk"; $13 = "NA"} {print $0}' | sed '1s/locus_name/libname/g' > "$outdir/_summary_stats/singlenuc_passQC_final_wCC.txt"


	# repeat SC3 clustering while omitting genes that vary significantly with CC (Fig. S7)
	# -------------------
	mkdir -p "$outdir/clustering/repeat_noCC"
	
	# drop the sig. genes from the count matrix
	awk -F$'\t' '{OFS=FS}ARGIND == 1 { del[$2]++ } ARGIND == 2 && !del[$1]' "$outdir/cell_cycle/traj_analysis_significant_genes.txt" "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" > "$outdir/clustering/repeat_noCC/all_counts_in_genes_dedup_passQC_noCC.txt"
	
	# all nuclei
	mkdir -p "$outdir/clustering/repeat_noCC/SC3_all"
	cmd="$path_to_scripts/single_cell_cluster_SC3.R $outdir/clustering/repeat_noCC/all_counts_in_genes_dedup_passQC_noCC.txt $outdir/_summary_stats/singlenuc_passQC_final.txt $outdir/clustering/repeat_noCC/SC3_all/SC3_all"
	bsub -o "$outdir/clustering/repeat_noCC/SC3_all/SC3_all_log.txt" -K "$cmd" & pid1=$!
	
	# CxV endo 4 DAP	
	mkdir -p "$outdir/clustering/repeat_noCC/CxV_endo_4DAP"	
	$path_to_scripts/subset_large_file.R "$outdir/clustering/repeat_noCC/all_counts_in_genes_dedup_passQC_noCC.txt" "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_counts_in_genes_dedup_passQC_noCC.txt" --cols "$outdir/clustering/SC3_subsets/CxV_endo_4DAP.txt"		
	$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_final.txt" "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_descfile.txt" --rows "$outdir/clustering/SC3_subsets/CxV_endo_4DAP.txt"
	
	cmd="$path_to_scripts/single_cell_cluster_SC3.R $outdir/clustering/repeat_noCC/CxV_endo_4DAP_counts_in_genes_dedup_passQC_noCC.txt $outdir/clustering/repeat_noCC/CxV_endo_4DAP_descfile.txt $outdir/clustering/repeat_noCC/CxV_endo_4DAP/CxV_endo_4DAP --maxK 15"
	bsub -o "$outdir/clustering/repeat_noCC/CxV_endo_4DAP/CxV_endo_4DAP_log.txt" -K "$cmd" & pid2=$!
	
	# VxC endo 4 DAP
	mkdir -p "$outdir/clustering/repeat_noCC/VxC_endo_4DAP"	
	$path_to_scripts/subset_large_file.R "$outdir/clustering/repeat_noCC/all_counts_in_genes_dedup_passQC_noCC.txt" "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_counts_in_genes_dedup_passQC_noCC.txt" --cols "$outdir/clustering/SC3_subsets/VxC_endo_4DAP.txt"		
	$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_final.txt" "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_descfile.txt" --rows "$outdir/clustering/SC3_subsets/VxC_endo_4DAP.txt"
	
	cmd="$path_to_scripts/single_cell_cluster_SC3.R $outdir/clustering/repeat_noCC/VxC_endo_4DAP_counts_in_genes_dedup_passQC_noCC.txt $outdir/clustering/repeat_noCC/VxC_endo_4DAP_descfile.txt $outdir/clustering/repeat_noCC/VxC_endo_4DAP/VxC_endo_4DAP --maxK 15"
	bsub -o "$outdir/clustering/repeat_noCC/VxC_endo_4DAP/VxC_endo_4DAP_log.txt" -K "$cmd" & pid3=$!

	# wait for all jobs to finish
	wait ${pid1} || err_msg "SC3 noCC run on all nuclei failed, see $outdir/clustering/repeat_noCC/SC3_all/SC3_all_log.txt" "$log"
	wait ${pid2} || err_msg "SC3 noCC run on CxV 4 DAP nuclei failed, see $outdir/clustering/repeat_noCC/CxV_endo_4DAP/CxV_endo_4DAP_log.txt" "$log"
	wait ${pid3} || err_msg "SC3 noCC run on VxC 4 DAP nuclei failed, see $outdir/clustering/repeat_noCC/VxC_endo_4DAP/VxC_endo_4DAP_log.txt" "$log"

	# make consensus plots prettier	
	ff="$outdir/clustering/repeat_noCC/SC3_all/SC3_all_21_consensus_matrix.txt"
	factorsfile="$outdir/clustering/SC3_all/newfactors.txt"
	bb=${ff%"_consensus_matrix.txt"}
	bbd=$( dirname "$bb" )		
	consensus="${bb}_consensus_matrix.txt"
	clusfile="${bb}_cluster_info.txt"
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig" --metadata "$factorsfile"
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig_nocluscolors" --metadata "$factorsfile" --nocluscolors

	ff="$outdir/clustering/repeat_noCC/CxV_endo_4DAP/CxV_endo_4DAP_9_consensus_matrix.txt"
	cut -f1,6,11 "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_descfile.txt" > "$outdir/clustering/repeat_noCC/CxV_endo_4DAP/newfactors.txt"
	factorsfile="$outdir/clustering/repeat_noCC/CxV_endo_4DAP/newfactors.txt"
	bb=${ff%"_consensus_matrix.txt"}
	bbd=$( dirname "$bb" )
		
	consensus="${bb}_consensus_matrix.txt"
	clusfile="${bb}_cluster_info.txt"
	
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig" --metadata "$factorsfile"
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig_nocluscolors" --metadata "$factorsfile" --nocluscolors
	
	ff="$outdir/clustering/repeat_noCC/VxC_endo_4DAP/VxC_endo_4DAP_7_consensus_matrix.txt"
	cut -f1,6,11 "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_descfile.txt" > "$outdir/clustering/repeat_noCC/VxC_endo_4DAP/newfactors.txt"
	factorsfile="$outdir/clustering/repeat_noCC/VxC_endo_4DAP/newfactors.txt"
	bb=${ff%"_consensus_matrix.txt"}
	bbd=$( dirname "$bb" )
		
	consensus="${bb}_consensus_matrix.txt"
	clusfile="${bb}_cluster_info.txt"
	
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig" --metadata "$factorsfile"
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig_nocluscolors" --metadata "$factorsfile" --nocluscolors
		
	# sub-cluster a few clusters
	# CxV endo 4DAP - redo cluster 9
	awk -F$'\t' '{OFS=FS} NR!=1 && $2 == 9 {print $1}' "$outdir/clustering/repeat_noCC/CxV_endo_4DAP/CxV_endo_4DAP_9_cluster_info.txt" > "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo.txt"
	mkdir -p "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo"
	
	$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo_counts.txt" --cols "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo.txt"		
	$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_final.txt" "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo_descfile.txt" --rows "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo.txt"

	$path_to_scripts/single_cell_cluster_SC3.R "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo_counts.txt" "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo_descfile.txt" "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo/CxV_endo_4DAP_9_redo" --maxK 8 #> "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo/CxV_endo_4DAP_9_redo_log.txt" 2>&1
	
	ff="$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo/CxV_endo_4DAP_9_redo_4_consensus_matrix.txt"
	cut -f1,6,11 "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo_descfile.txt" > "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo/newfactors.txt"
	factorsfile="$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo/newfactors.txt"
	bb=${ff%"_consensus_matrix.txt"}
	bbd=$( dirname "$bb" )
		
	consensus="${bb}_consensus_matrix.txt"
	clusfile="${bb}_cluster_info.txt"
	
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig" --metadata "$factorsfile"
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig_nocluscolors" --metadata "$factorsfile" --nocluscolors
	
	# VxC endo 4DAP - redo clusters 1-3
	awk -F$'\t' '{OFS=FS} NR!=1 && ($2 == 1 || $2 == 2 || $2 == 3) {print $1}' "$outdir/clustering/repeat_noCC/VxC_endo_4DAP/VxC_endo_4DAP_7_cluster_info.txt" > "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo.txt"
	mkdir -p "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo"
	
	$path_to_scripts/subset_large_file.R "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo_counts.txt" --cols "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo.txt"		
	$path_to_scripts/subset_large_file.R "$summ/singlenuc_passQC_final.txt" "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo_descfile.txt" --rows "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo.txt"

	$path_to_scripts/single_cell_cluster_SC3.R "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo_counts.txt" "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo_descfile.txt" "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo/VxC_endo_4DAP_1to3_redo" --maxK 8 #> "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo/VxC_endo_4DAP_1to3_redo_log.txt" 2>&1
		
	ff="$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo/VxC_endo_4DAP_1to3_redo_5_consensus_matrix.txt"
	cut -f1,6,11 "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo_descfile.txt" > "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo/newfactors.txt"
	factorsfile="$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo/newfactors.txt"
	bb=${ff%"_consensus_matrix.txt"}
	bbd=$( dirname "$bb" )
		
	consensus="${bb}_consensus_matrix.txt"
	clusfile="${bb}_cluster_info.txt"
	
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig" --metadata "$factorsfile"
	python2 $path_to_scripts/plot_SC3_heatmap.py "$consensus" "$clusfile" "${bb}_consensus_orig_nocluscolors" --metadata "$factorsfile" --nocluscolors
	
	# Get final lists of clusters
	# CxV 4DAP endo
	awk '$2 != 9' "$outdir/clustering/repeat_noCC/CxV_endo_4DAP/CxV_endo_4DAP_9_cluster_info.txt" > "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_final_clusters.txt"
	awk -F$'\t' '{OFS=FS} NR != 1 {print $1,$2+8}' "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_9_redo/CxV_endo_4DAP_9_redo_4_cluster_info.txt" >> "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_final_clusters.txt"
	
	cat "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_1to3_redo/VxC_endo_4DAP_1to3_redo_5_cluster_info.txt" > "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_final_clusters.txt"
	awk -F$'\t' '{OFS=FS} NR != 1 && $2 > 3 {print $1,$2+2}' "$outdir/clustering/repeat_noCC/VxC_endo_4DAP/VxC_endo_4DAP_7_cluster_info.txt" >> "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_final_clusters.txt"
	
	# compare clusterings obtained with and without cell cycle-correlated genes (Fig. S7)
	# -------------------
	# Test overlap between these new clusters and the old ones (that included the CC-correlated genes)
	
	# CxV 4DAP endo
	mkdir -p "$outdir/clustering/repeat_noCC/gene_overlaps/CxV_endo_4DAP"
	ff_wCC=""; nn_wCC=""
	for ((i=1;i<=14;++i)); do
		awk -v a="$i" '$2 == a {print $1}' "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt" > "$outdir/clustering/repeat_noCC/gene_overlaps/CxV_endo_4DAP/cluster${i}_wCC.txt"
		ff_wCC="${ff_wCC},${outdir}/clustering/repeat_noCC/gene_overlaps/CxV_endo_4DAP/cluster${i}_wCC.txt"; nn_wCC="${nn_wCC},${i}_wCC"
	done
	ff_wCC="${ff_wCC:1}"; nn_wCC="${nn_wCC:1}"	

	ff_noCC=""; nn_noCC=""
	for ((i=1;i<=12;++i)); do
		awk -v a="$i" '$2 == a {print $1}' "$outdir/clustering/repeat_noCC/CxV_endo_4DAP_final_clusters.txt" > "$outdir/clustering/repeat_noCC/gene_overlaps/CxV_endo_4DAP/cluster${i}_noCC.txt"
		ff_noCC="${ff_noCC},${outdir}/clustering/repeat_noCC/gene_overlaps/CxV_endo_4DAP/cluster${i}_noCC.txt"; nn_noCC="${nn_noCC},${i}_noCC"
	done
	ff_noCC="${ff_noCC:1}"; nn_noCC="${nn_noCC:1}"	
	
	$path_to_scripts/gene_overlaps_dotplot.R "$ff_wCC" "$ff_noCC" "$outdir/clustering/repeat_noCC/gene_overlaps/CxV_endo_4DAP_dotplot" --namesA $nn_wCC --namesB $nn_noCC --popsize 399 --scoretype A --width 20 --height 15 --sizeupper 50
		
	# VxC 4DAP endo
	mkdir -p "$outdir/clustering/repeat_noCC/gene_overlaps/VxC_endo_4DAP"
	ff_wCC=""; nn_wCC=""
	for ((i=1;i<=11;++i)); do
		awk -v a="$i" '$2 == a {print $1}' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt" > "$outdir/clustering/repeat_noCC/gene_overlaps/VxC_endo_4DAP/cluster${i}_wCC.txt"
		ff_wCC="${ff_wCC},${outdir}/clustering/repeat_noCC/gene_overlaps/VxC_endo_4DAP/cluster${i}_wCC.txt"; nn_wCC="${nn_wCC},${i}_wCC"
	done
	ff_wCC="${ff_wCC:1}"; nn_wCC="${nn_wCC:1}"	

	ff_noCC=""; nn_noCC=""
	for ((i=1;i<=9;++i)); do
		awk -v a="$i" '$2 == a {print $1}' "$outdir/clustering/repeat_noCC/VxC_endo_4DAP_final_clusters.txt" > "$outdir/clustering/repeat_noCC/gene_overlaps/VxC_endo_4DAP/cluster${i}_noCC.txt"
		ff_noCC="${ff_noCC},${outdir}/clustering/repeat_noCC/gene_overlaps/VxC_endo_4DAP/cluster${i}_noCC.txt"; nn_noCC="${nn_noCC},${i}_noCC"
	done
	ff_noCC="${ff_noCC:1}"; nn_noCC="${nn_noCC:1}"	
	
	$path_to_scripts/gene_overlaps_dotplot.R "$ff_wCC" "$ff_noCC" "$outdir/clustering/repeat_noCC/gene_overlaps/VxC_endo_4DAP_dotplot" --namesA $nn_wCC --namesB $nn_noCC --popsize 399 --scoretype A --width 15 --height 10 --sizeupper 50

	# all nuclei
	mkdir -p "$outdir/clustering/repeat_noCC/gene_overlaps/SC3_all"
	ff_wCC=""; nn_wCC=""
	for ((i=1;i<=23;++i)); do
		awk -v a="$i" '$2 == a {print $1}' "$outdir/clustering/SC3_all/SC3_all_23_cluster_info.txt" > "$outdir/clustering/repeat_noCC/gene_overlaps/SC3_all/cluster${i}_wCC.txt"
		ff_wCC="${ff_wCC},${outdir}/clustering/repeat_noCC/gene_overlaps/SC3_all/cluster${i}_wCC.txt"; nn_wCC="${nn_wCC},${i}_wCC"
	done
	ff_wCC="${ff_wCC:1}"; nn_wCC="${nn_wCC:1}"	

	ff_noCC=""; nn_noCC=""
	for ((i=1;i<=21;++i)); do
		awk -v a="$i" '$2 == a {print $1}' "$outdir/clustering/repeat_noCC/SC3_all/SC3_all_21_cluster_info.txt" > "$outdir/clustering/repeat_noCC/gene_overlaps/SC3_all/cluster${i}_noCC.txt"
		ff_noCC="${ff_noCC},${outdir}/clustering/repeat_noCC/gene_overlaps/SC3_all/cluster${i}_noCC.txt"; nn_noCC="${nn_noCC},${i}_noCC"
	done
	ff_noCC="${ff_noCC:1}"; nn_noCC="${nn_noCC:1}"	
	
	$path_to_scripts/gene_overlaps_dotplot.R "$ff_wCC" "$ff_noCC" "$outdir/clustering/repeat_noCC/gene_overlaps/SC3_all_dotplot" --namesA $nn_wCC --namesB $nn_noCC --popsize 399 --scoretype A --width 20 --height 15 --sizeupper 50

fi


# ----------------------
# Step 6: characterizing the CxV, VxC 4 DAP clusters by identifying DE genes
# ----------------------
if "$skip6"; then printf "\nSkipping step 6 by user request...\n" | tee -a "$log"
else printf "\nStep 6: characterizing the CxV, VxC 4 DAP clusters by identifying DE genes\n" | tee -a "$log"; ts=$(date +%s)
	
	# identify DE genes using DEsingle (https://bioconductor.org/packages/release/bioc/html/DEsingle.html)
	# -------------------

	# differential expression analysis between clusters (all pairwise); do several sets of comparisons:
	#  - all CxV endo clusters against themselves
	#  - all CxV endo batches/plates against themselves (control for technical variability)
	#  - all VxC endo clusters against themselves
	#  - all VxC endo batches/plates against themselves (control for technical variability)
	#  - all CxV endo clusters against VxC endo clusters
	#  - all CxV seedcoat clusters against themselves
	#  - all CxV seedcoat batches/plates against themselves (control for technical variability)
	#  - all VxC seedcoat clusters against themselves
	#  - all VxC seedcoat batches/plates against themselves (control for technical variability)
	#  - all CxV seedcoat clusters against VxC seedcoat clusters
	#  - all CxV endo clusters against all CxV seedcoat clusters
	#  - all VxC endo clusters against all VxC seedcoat clusters

	mkdir -p "$outdir/DE_analysis"
	echo "stubname	cluster" > "$outdir/DE_analysis/all_endo_clusters.txt"
	cat <( tail -n+2 "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt" ) <( tail -n+2 "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt" ) >> "$outdir/DE_analysis/all_endo_clusters.txt"	
	$path_to_scripts/merge_by_column.R "$outdir/DE_analysis/all_endo_clusters.txt" "$outdir/_summary_stats/singlenuc_passQC_final.txt" stubname "-" --tokeep allx | awk -F$'\t' '{OFS=FS} {print $1,$5"_"$2,$9}' > "$outdir/DE_analysis/all_endo_factors.txt"
		
	mkdir -p "$lsf/DE_analysis"
	mkdir -p "$outdir/DE_analysis/cluslists"
	
	# CxV clusters
	pids=()
	echo "DE analyses of CxV endo clusters"
	mkdir -p "$outdir/DE_analysis/CxV_endo_4DAP_within"
	for ((i=1;i<14;++i)); do
		for ((j=$(( $i + 1 ));j<=14;++j)); do

			# get lists of nuclei in these clusters
			[ ! -f "$outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${i}.txt" ] && { awk -v ii="${i}" '$2 == ii {print $1}' "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt" > "$outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${i}.txt"; }
			[ ! -f "$outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${j}.txt" ] && { awk -v ii="${j}" '$2 == ii {print $1}' "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt" > "$outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${j}.txt"; }
		
			# compare these two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/clustering/SC3_subsets/CxV_endo_4DAP/CxV_endo_4DAP_counts.txt $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${i}.txt $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${j}.txt $outdir/DE_analysis/CxV_endo_4DAP_within/CxV_C${i}_vs_CxV_C${j}"
			bsub -o "$lsf/DE_analysis/CxV_endo_C${i}_vs_C${j}.txt" -K "$cmd" & pids+=( $! )
		done
	done
	
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in $lsf/DE_analysis/CxV_endo_C*" "$log"
	done
	
	# batch effects that may confound DE analysis, so identify DE genes pairwise between all batches too
	pids=()
	echo "DE analyses of CxV endo batch effects"
	mkdir -p "$outdir/DE_analysis/CxV_endo_4DAP_within_batch"
	
	echo "geneID	batch" > "$outdir/DE_analysis/batches.txt"
	cut -f1,8 "$outdir/_summary_stats/singlenuc_passQC_final.txt" | tail -n+2 >> "$outdir/DE_analysis/batches.txt"
	
	$path_to_scripts/merge_by_column.R "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt" "$outdir/DE_analysis/batches.txt" geneID - --tokeep allx | cut -f1,3 > "$outdir/DE_analysis/CxV_endo_4DAP_within_batch/batchlist.txt"
	batchlist=( $( tail -n+2 "$outdir/DE_analysis/CxV_endo_4DAP_within_batch/batchlist.txt" | cut -f2 | sort -k 1n,1 | uniq ) )
	
	for ((i=0;i<${#batchlist[@]}-1;++i)); do
		for ((j=$(( $i + 1 ));j<${#batchlist[@]};++j)); do
			# get lists of nuclei in these batches
			[ ! -f "$outdir/DE_analysis/cluslists/CxV_4DAP_endo_B${batchlist[i]}.txt" ] && { awk -v ii="${batchlist[i]}" '$2 == ii {print $1}' "$outdir/DE_analysis/CxV_endo_4DAP_within_batch/batchlist.txt" > "$outdir/DE_analysis/cluslists/CxV_4DAP_endo_B${batchlist[i]}.txt"; }
			[ ! -f "$outdir/DE_analysis/cluslists/CxV_4DAP_endo_B${batchlist[j]}.txt" ] && { awk -v ii="${batchlist[j]}" '$2 == ii {print $1}' "$outdir/DE_analysis/CxV_endo_4DAP_within_batch/batchlist.txt" > "$outdir/DE_analysis/cluslists/CxV_4DAP_endo_B${batchlist[j]}.txt"; }
		
			# compare the two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/clustering/SC3_subsets/CxV_endo_4DAP/CxV_endo_4DAP_counts.txt $outdir/DE_analysis/cluslists/CxV_4DAP_endo_B${batchlist[i]}.txt $outdir/DE_analysis/cluslists/CxV_4DAP_endo_B${batchlist[j]}.txt $outdir/DE_analysis/CxV_endo_4DAP_within_batch/CxV_B${batchlist[i]}_vs_CxV_B${batchlist[j]}"
			bsub -o "$lsf/DE_analysis/CxV_endo_B${batchlist[i]}_vs_CxV_B${batchlist[j]}.txt" -K "$cmd" & pids+=( $! )
		done
	done
	
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in $lsf/DE_analysis/CxV_endo_B*" "$log"
	done

	# VxC clusters
	pids=()
	echo "DE analyses of VxC endo clusters"
	mkdir -p "$outdir/DE_analysis/VxC_endo_4DAP_within"
	for ((i=1;i<11;++i)); do
		for ((j=$(( $i + 1 ));j<=11;++j)); do

			# get lists of nuclei in these clusters
			[ ! -f "$outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${i}.txt" ] && { awk -v ii="${i}" '$2 == ii {print $1}' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt" > "$outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${i}.txt"; }
			[ ! -f "$outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${j}.txt" ] && { awk -v ii="${j}" '$2 == ii {print $1}' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt" > "$outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${j}.txt"; }
		
			# compare these two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/clustering/SC3_subsets/VxC_endo_4DAP/VxC_endo_4DAP_counts.txt $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${i}.txt $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${j}.txt $outdir/DE_analysis/VxC_endo_4DAP_within/VxC_C${i}_vs_VxC_C${j}"
			bsub -o "$lsf/DE_analysis/VxC_endo_C${i}_vs_C${j}.txt" -K "$cmd" & pids+=( $! )
		done
	done
	
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in $lsf/DE_analysis/VxC_endo_C*" "$log"
	done
	
	# batch effects that may confound DE analysis (VxC endo)
	pids=()
	echo "DE analyses of VxC endo batch effects"
	mkdir -p "$outdir/DE_analysis/VxC_endo_4DAP_within_batch"
		
	$path_to_scripts/merge_by_column.R "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt" "$outdir/DE_analysis/batches.txt" geneID - --tokeep allx | cut -f1,3 > "$outdir/DE_analysis/VxC_endo_4DAP_within_batch/batchlist.txt"
	batchlist=( $( tail -n+2 "$outdir/DE_analysis/VxC_endo_4DAP_within_batch/batchlist.txt" | cut -f2 | sort -k 1n,1 | uniq ) )
	
	for ((i=0;i<${#batchlist[@]}-1;++i)); do
		for ((j=$(( $i + 1 ));j<${#batchlist[@]};++j)); do
			# get lists of nuclei in these batches
			[ ! -f "$outdir/DE_analysis/cluslists/VxC_4DAP_endo_B${batchlist[i]}.txt" ] && { awk -v ii="${batchlist[i]}" '$2 == ii {print $1}' "$outdir/DE_analysis/VxC_endo_4DAP_within_batch/batchlist.txt" > "$outdir/DE_analysis/cluslists/VxC_4DAP_endo_B${batchlist[i]}.txt"; }
			[ ! -f "$outdir/DE_analysis/cluslists/VxC_4DAP_endo_B${batchlist[j]}.txt" ] && { awk -v ii="${batchlist[j]}" '$2 == ii {print $1}' "$outdir/DE_analysis/VxC_endo_4DAP_within_batch/batchlist.txt" > "$outdir/DE_analysis/cluslists/VxC_4DAP_endo_B${batchlist[j]}.txt"; }
		
			# compare the two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/clustering/SC3_subsets/VxC_endo_4DAP/VxC_endo_4DAP_counts.txt $outdir/DE_analysis/cluslists/VxC_4DAP_endo_B${batchlist[i]}.txt $outdir/DE_analysis/cluslists/VxC_4DAP_endo_B${batchlist[j]}.txt $outdir/DE_analysis/VxC_endo_4DAP_within_batch/VxC_B${batchlist[i]}_vs_VxC_B${batchlist[j]}"
			bsub -o "$lsf/DE_analysis/VxC_endo_B${batchlist[i]}_vs_VxC_B${batchlist[j]}.txt" -K "$cmd" & pids+=( $! )
		done
	done
	
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in $lsf/DE_analysis/VxC_endo_B*" "$log"
	done

	# all CxV clusters against all VxC clusters
	echo "DE analyses of between cross endo clusters"
	mkdir -p "$outdir/DE_analysis/CxV_and_VxC_endo_4DAP_between"
	pids=()
	for ((i=1;i<=14;++i)); do
		for ((j=1;j<=11;++j)); do		
			# compare these two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${i}.txt $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${j}.txt $outdir/DE_analysis/CxV_and_VxC_endo_4DAP_between/CxV_C${i}_vs_VxC_C${j}"
			bsub -o "$lsf/DE_analysis/btw_clus_endo_C${i}_vs_C${j}.txt" -K "$cmd" & pids+=( $! )
		done
	done

	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in $lsf/DE_analysis/btw_clus_endo*" "$log"
	done
	
	# CxV seed coat clusters
	pids=()
	echo "DE analyses of CxV seedcoat clusters"
	mkdir -p "$outdir/DE_analysis/CxV_seedcoat_4DAP_within"
	for ((i=1;i<6;++i)); do
		for ((j=$(( $i + 1 ));j<=6;++j)); do

			# get lists of nuclei in these clusters
			[ ! -f "$outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_C${i}.txt" ] && { awk -v ii="${i}" '$2 == ii {print $1}' "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt" > "$outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_C${i}.txt"; }
			[ ! -f "$outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_C${j}.txt" ] && { awk -v ii="${j}" '$2 == ii {print $1}' "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt" > "$outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_C${j}.txt"; }
		
			# compare these two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP/CxV_seedcoat_4DAP_counts.txt $outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_C${i}.txt $outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_C${j}.txt $outdir/DE_analysis/CxV_seedcoat_4DAP_within/CxV_C${i}_vs_CxV_C${j}"
			bsub -o "$lsf/DE_analysis/CxV_seedcoat_C${i}_vs_C${j}.txt" -K "$cmd" & pids+=( $! )
		done
	done
	
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in$lsf/DE_analysis/CxV_C*" "$log"
	done
	
	# batch effects that may confound DE analysis (CxV seed coat)
	pids=()
	echo "DE analyses of CxV seedcoat batch effects"
	mkdir -p "$outdir/DE_analysis/CxV_seedcoat_4DAP_within_batch"
		
	$path_to_scripts/merge_by_column.R "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt" "$outdir/DE_analysis/batches.txt" geneID - --tokeep allx | cut -f1,3 > "$outdir/DE_analysis/CxV_seedcoat_4DAP_within_batch/batchlist.txt"
	batchlist=( $( tail -n+2 "$outdir/DE_analysis/CxV_seedcoat_4DAP_within_batch/batchlist.txt" | cut -f2 | sort -k 1n,1 | uniq ) )
	
	for ((i=0;i<${#batchlist[@]}-1;++i)); do
		for ((j=$(( $i + 1 ));j<${#batchlist[@]};++j)); do
			# get lists of nuclei in these batches
			[ ! -f "$outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_B${batchlist[i]}.txt" ] && { awk -v ii="${batchlist[i]}" '$2 == ii {print $1}' "$outdir/DE_analysis/CxV_seedcoat_4DAP_within_batch/batchlist.txt" > "$outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_B${batchlist[i]}.txt"; }
			[ ! -f "$outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_B${batchlist[j]}.txt" ] && { awk -v ii="${batchlist[j]}" '$2 == ii {print $1}' "$outdir/DE_analysis/CxV_seedcoat_4DAP_within_batch/batchlist.txt" > "$outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_B${batchlist[j]}.txt"; }
		
			# compare the two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP/CxV_seedcoat_4DAP_counts.txt $outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_B${batchlist[i]}.txt $outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_B${batchlist[j]}.txt $outdir/DE_analysis/CxV_seedcoat_4DAP_within_batch/CxV_B${batchlist[i]}_vs_CxV_B${batchlist[j]}"
			bsub -o "$lsf/DE_analysis/CxV_B${batchlist[i]}_vs_CxV_B${batchlist[j]}.txt" -K "$cmd" & pids+=( $! )
		done
	done
	
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in$lsf/DE_analysis/CxV_B*" "$log"
	done
	
	# VxC seedcoat clusters
	pids=()
	echo "DE analyses of VxC seedcoat clusters"
	mkdir -p "$outdir/DE_analysis/VxC_seedcoat_4DAP_within"
	for ((i=1;i<8;++i)); do
		for ((j=$(( $i + 1 ));j<=8;++j)); do

			# get lists of nuclei in these clusters
			[ ! -f "$outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_C${i}.txt" ] && { awk -v ii="${i}" '$2 == ii {print $1}' "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt" > "$outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_C${i}.txt"; }
			[ ! -f "$outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_C${j}.txt" ] && { awk -v ii="${j}" '$2 == ii {print $1}' "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt" > "$outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_C${j}.txt"; }
		
			# compare these two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP/VxC_seedcoat_4DAP_counts.txt $outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_C${i}.txt $outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_C${j}.txt $outdir/DE_analysis/VxC_seedcoat_4DAP_within/VxC_C${i}_vs_VxC_C${j}"
			bsub -o "$lsf/DE_analysis/VxC_seedcoat_C${i}_vs_C${j}.txt" -K "$cmd" & pids+=( $! )
		done
	done
	
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in$lsf/DE_analysis/VxC_C*" "$log"
	done
	
	# batch effects that may confound DE analysis (VxC seedcoat)
	pids=()
	echo "DE analyses of VxC seedcoat batch effects"
	mkdir -p "$outdir/DE_analysis/VxC_seedcoat_4DAP_within_batch"
	
	$path_to_scripts/merge_by_column.R "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt" "$outdir/DE_analysis/batches.txt" geneID - --tokeep allx | cut -f1,3 > "$outdir/DE_analysis/VxC_seedcoat_4DAP_within_batch/batchlist.txt"
	batchlist=( $( tail -n+2 "$outdir/DE_analysis/VxC_seedcoat_4DAP_within_batch/batchlist.txt" | cut -f2 | sort -k 1n,1 | uniq ) )
	
	for ((i=0;i<${#batchlist[@]}-1;++i)); do
		for ((j=$(( $i + 1 ));j<${#batchlist[@]};++j)); do
			# get lists of nuclei in these batches
			[ ! -f "$outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_B${batchlist[i]}.txt" ] && { awk -v ii="${batchlist[i]}" '$2 == ii {print $1}' "$outdir/DE_analysis/VxC_seedcoat_4DAP_within_batch/batchlist.txt" > "$outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_B${batchlist[i]}.txt"; }
			[ ! -f "$outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_B${batchlist[j]}.txt" ] && { awk -v ii="${batchlist[j]}" '$2 == ii {print $1}' "$outdir/DE_analysis/VxC_seedcoat_4DAP_within_batch/batchlist.txt" > "$outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_B${batchlist[j]}.txt"; }
		
			# compare the two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP/VxC_seedcoat_4DAP_counts.txt $outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_B${batchlist[i]}.txt $outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_B${batchlist[j]}.txt $outdir/DE_analysis/VxC_seedcoat_4DAP_within_batch/VxC_B${batchlist[i]}_vs_VxC_B${batchlist[j]}"
			bsub -o "$lsf/DE_analysis/VxC_B${batchlist[i]}_vs_VxC_B${batchlist[j]}.txt" -K "$cmd" & pids+=( $! )
		done
	done
	
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in$lsf/DE_analysis/VxC_B*" "$log"
	done
	
	# all CxV clusters against all VxC clusters
	echo "DE analyses of between cross seedcoat clusters"
	mkdir -p "$outdir/DE_analysis/CxV_and_VxC_seedcoat_4DAP_between"
	pids=()
	for ((i=1;i<=6;++i)); do
		for ((j=1;j<=8;++j)); do		
			# compare these two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt $outdir/DE_analysis/cluslists/CxV_4DAP_seedcoat_C${i}.txt $outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_C${j}.txt $outdir/DE_analysis/CxV_and_VxC_seedcoat_4DAP_between/CxV_C${i}_vs_VxC_C${j}"
			bsub -o "$lsf/DE_analysis/btw_clus_seedcoat_C${i}_vs_C${j}.txt" -K "$cmd" & pids+=( $! )
		done
	done

	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in $lsf/DE_analysis/btw_clus_seedcoat*" "$log"
	done

	# all CxV endo clusters against all CxV seedcoat
	echo "DE analyses of between tissue CxV clusters"
	mkdir -p "$outdir/DE_analysis/CxV_endo_vs_seedcoat_between"
	pids=()
	for ((i=1;i<=14;++i)); do
		for ((j=1;j<=6;++j)); do		
			# compare these two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${i}.txt $outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_C${j}.txt $outdir/DE_analysis/CxV_endo_vs_seedcoat_between/CxV_endo_C${i}_vs_seedcoat_C${j}"
			bsub -o "$lsf/DE_analysis/btw_clus_endo_C${i}_vs_seedcoat_C${j}.txt" -K "$cmd" & pids+=( $! )
		done
	done

	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in $lsf/DE_analysis/btw_clus_endo_C*" "$log"
	done
	
	# all VxC endo clusters against all VxC seedcoat
	echo "DE analyses of between tissue VxC clusters"
	mkdir -p "$outdir/DE_analysis/VxC_endo_vs_seedcoat_between"
	pids=()
	for ((i=1;i<=11;++i)); do
		for ((j=7;j<=8;++j)); do		
			# compare these two sets of nuclei using DEsingle
			cmd="$path_to_scripts/run_DEsingle.R $outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${i}.txt $outdir/DE_analysis/cluslists/VxC_4DAP_seedcoat_C${j}.txt $outdir/DE_analysis/VxC_endo_vs_seedcoat_between/VxC_endo_C${i}_vs_seedcoat_C${j}"
			bsub -o "$lsf/DE_analysis/btw_clus_endo_C${i}_vs_seedcoat_C${j}.txt" -K "$cmd" & pids+=( $! )
		done
	done

	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "run_DEsingle failed, see log files in $lsf/DE_analysis/btw_clus_endo_C*" "$log"
	done


	# Pull out DE genes across all these comparisons
	# -------------------
	# Across all comparisons, get smallest p-value for each gene and largest abs(log2FC)
	# pool together all cluster lists from CxV, VxC endo and seedcoat
	cat <( awk -F$'\t' '{OFS=FS} NR==1 {print $0} NR!=1 {print $1,"CxV_endo_"$2}' "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt" ) <( awk -F$'\t' '{OFS=FS} NR!=1 {print $1,"VxC_endo_"$2}' "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt" ) <( awk -F$'\t' '{OFS=FS} NR!=1 {print $1,"CxV_seedcoat_"$2}' "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt" ) <( awk -F$'\t' '{OFS=FS} NR!=1 {print $1,"VxC_seedcoat_"$2}' "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt" ) > "$outdir/clustering/SC3_subsets/all_4DAP_nuclei_clusters.txt"	
	
	# first, get all comparisons that pass selected cutoffs
	# a single comparison must pass both the pval and log2fc cutoffs for that gene to pass overall
	# note - cutoff settings set at top of script	
	mkdir -p $outdir/DE_analysis/summaries

	rm -f "$outdir/DE_analysis/summaries/all_batch_effects.txt"; rm -f "$outdir/DE_analysis/summaries/all_cluster_effects.txt"
	for ff in $outdir/DE_analysis/*/*_all.txt; do
			printf "appending $ff"
		if [[ $ff  == *"batch"* ]]; then
				printf " -> batch\n"
			cut -f1,15,22 "$ff" | awk -F$'\t' -v p="$pvalcutoff" -v l="$log2fccutoff" '(log($2)/log(2) > l || log($2)/log(2) < -1*l || $2 == 0 || $2 == "Inf") && $3 < p' >> "$outdir/DE_analysis/summaries/all_batch_effects.txt"
		else
				printf " -> cluster\n"
			cut -f1,15,22 "$ff" | awk -F$'\t' -v p="$pvalcutoff" -v l="$log2fccutoff" '(log($2)/log(2) > l || log($2)/log(2) < -1*l || $2 == 0 || $2 == "Inf") && $3 < p' >> "$outdir/DE_analysis/summaries/all_cluster_effects.txt"
		fi
	done
	
	# get all genes that pass these criteria among the cluster comparisons and among the batch/plate comparisons
	cut -f1 "$outdir/DE_analysis/summaries/all_batch_effects.txt" | sort -u > "$outdir/DE_analysis/summaries/genes_pass_batch.txt"
	cut -f1 "$outdir/DE_analysis/summaries/all_cluster_effects.txt" | sort -u > "$outdir/DE_analysis/summaries/genes_pass_cluster.txt"
	
	# get all cluster-DE genes not in the batch DE list
	diff "$outdir/DE_analysis/summaries/genes_pass_cluster.txt" "$outdir/DE_analysis/summaries/genes_pass_batch.txt" --new-line-format="" --unchanged-line-format="" > "$outdir/DE_analysis/summaries/genes_pass_cluster_nobatch.txt"
	
	numC=$( cat "$outdir/DE_analysis/summaries/genes_pass_cluster.txt" | wc -l )
	numB=$( cat "$outdir/DE_analysis/summaries/genes_pass_batch.txt" | wc -l )
	numF=$( cat "$outdir/DE_analysis/summaries/genes_pass_cluster_nobatch.txt" | wc -l )
	
	echo "For p-value cutoff ${pvalcutoffs[i]} and log2FC cutoff ${log2fccutoffs[i]}:"
	echo " - Cluster-specific DE genes: $numC"
	echo " - Batch-specific DE genes: $numB"
	echo " - Cluster(-batch)-specific DE genes: $numF"

	# Result (03/11/2020):
	# For p-value cutoff  and log2FC cutoff :
	#  - Cluster-specific DE genes: 4540
	#  - Batch-specific DE genes: 97
	#  - Cluster(-batch)-specific DE genes: 4500


	# Separately, identify genes significantly DE between endosperm and seed coat
	# -------------------
	# do analysis in CxV 4 DAP and VxC 4 DAP separately
	mkdir -p $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP
	awk -F$'\t' '{OFS=FS} $5 == 4 && $4 == "CxV" && $10 == "endo" {print $1,$4}' $outdir/_summary_stats/singlenuc_passQC_final.txt > $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_4DAP_endo_nuclei.txt
	awk -F$'\t' '{OFS=FS} $5 == 4 && $4 == "CxV" && $10 == "seedcoat" {print $1,$4}' $outdir/_summary_stats/singlenuc_passQC_final.txt > $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_4DAP_seedcoat_nuclei.txt
	awk -F$'\t' '{OFS=FS} $5 == 4 && $4 == "VxC" && $10 == "endo" {print $1,$4}' $outdir/_summary_stats/singlenuc_passQC_final.txt > $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_4DAP_endo_nuclei.txt
	awk -F$'\t' '{OFS=FS} $5 == 4 && $4 == "VxC" && $10 == "seedcoat" {print $1,$4}' $outdir/_summary_stats/singlenuc_passQC_final.txt > $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_4DAP_seedcoat_nuclei.txt
	
	cmd1="$path_to_scripts/run_DEsingle.R $outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_4DAP_endo_nuclei.txt $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_4DAP_seedcoat_nuclei.txt $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_4DAP_endo_vs_seedcoat"
	cmd2="$path_to_scripts/run_DEsingle.R $outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_4DAP_endo_nuclei.txt $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_4DAP_seedcoat_nuclei.txt $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_4DAP_endo_vs_seedcoat"
	bsub -o "$lsf/DE_analysis/CxV_4DAP_endo_vs_seedcoat.txt" -K "$cmd1" & pid1=$!
	bsub -o "$lsf/DE_analysis/VxC_4DAP_endo_vs_seedcoat.txt" -K "$cmd2" & pid2=$!
	
	wait ${pid1} || err_msg "run_DEsingle failed, see $lsf/DE_analysis/CxV_4DAP_endo_vs_seedcoat.txt" "$log"
	wait ${pid2} || err_msg "run_DEsingle failed, see $lsf/DE_analysis/VxC_4DAP_endo_vs_seedcoat.txt" "$log"
	
	# use less stringent cutoff compared to above
	cut -f1,15,22 $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_4DAP_endo_vs_seedcoat_all.txt | awk -F$'\t' -v p="$pvalcutoff_sc" -v l="$log2fccutoff_sc" 'log($2)/log(2) > l && $3 < p {print $1}' > "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_genes_endo_gr_seedcoat.txt"
	cut -f1,15,22 $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_4DAP_endo_vs_seedcoat_all.txt | awk -F$'\t' -v p="$pvalcutoff_sc" -v l="$log2fccutoff_sc" 'log($2)/log(2) < -1*l && $3 < p {print $1}' > "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_genes_endo_le_seedcoat.txt"
	numA=$( cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_genes_endo_gr_seedcoat.txt" | wc -l )
	numB=$( cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_genes_endo_le_seedcoat.txt" | wc -l )
	cut -f1,15,22 $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_4DAP_endo_vs_seedcoat_all.txt | awk -F$'\t' -v p="$pvalcutoff_sc" -v l="$log2fccutoff_sc" 'log($2)/log(2) > l && $3 < p {print $1}' > "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_genes_endo_gr_seedcoat.txt"
	cut -f1,15,22 $outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_4DAP_endo_vs_seedcoat_all.txt | awk -F$'\t' -v p="$pvalcutoff_sc" -v l="$log2fccutoff_sc" 'log($2)/log(2) < -1*l && $3 < p {print $1}' > "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_genes_endo_le_seedcoat.txt"
	numC=$( cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_genes_endo_gr_seedcoat.txt" | wc -l )
	numD=$( cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_genes_endo_le_seedcoat.txt" | wc -l )
	
	# get intersection (all genes endo > seedcoat, etc. in both crosses)
	awk 'NR==FNR { lines[$0]=1; next } $0 in lines' "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_genes_endo_gr_seedcoat.txt" "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_genes_endo_gr_seedcoat.txt" > "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/intersect_genes_endo_gr_seedcoat.txt"
	awk 'NR==FNR { lines[$0]=1; next } $0 in lines' "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_genes_endo_le_seedcoat.txt" "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_genes_endo_le_seedcoat.txt" > "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/intersect_genes_endo_le_seedcoat.txt"
	numE=$( cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/intersect_genes_endo_gr_seedcoat.txt" | wc -l )
	numF=$( cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/intersect_genes_endo_le_seedcoat.txt" | wc -l )	
	
	# get union (all genes endo > seedcoat, etc. in either cross)
	cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_genes_endo_gr_seedcoat.txt" "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_genes_endo_gr_seedcoat.txt" | sort -u > "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/union_genes_endo_gr_seedcoat.txt"
	cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/CxV_genes_endo_le_seedcoat.txt" "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/VxC_genes_endo_le_seedcoat.txt" | sort -u > "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/union_genes_endo_le_seedcoat.txt"
	numG=$( cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/union_genes_endo_gr_seedcoat.txt" | wc -l )
	numH=$( cat "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/union_genes_endo_le_seedcoat.txt" | wc -l )	
	
	echo ""
	echo "For comparing endo and seedcoat, with p-value cutoff $pvalcutoff_sc and log2FC cutoff $log2fccutoff_sc:"
	echo "Genes endo > seedcoat:"
	echo " - CxV: $numA"
	echo " - VxC: $numC"
	echo " - intersect: $numE"
	echo " - union: $numG"
	echo ""
	echo "Genes endo < seedcoat:"
	echo " - CxV: $numB"
	echo " - VxC: $numD"
	echo " - intersect: $numF"
	echo " - union: $numH"
	
	# Result (03/12/2020):
	# For comparing endo and seedcoat, with p-value cutoff  and log2FC cutoff :
	# Genes endo > seedcoat:
	#  - CxV: 1482
	#  - VxC: 2500
	#  - intersect: 802
	#  - union: 3180
	# 
	# Genes endo < seedcoat:
	#  - CxV: 859
	#  - VxC: 975
	#  - intersect: 442
	#  - union: 1392
		
fi


# ----------------------
# Step 7: identify imprinted genes from endosperm nuclei
# ----------------------
if "$skip7"; then printf "\nSkipping step 7 (identify imprinted genes from endosperm nuclei) by user request...\n" | tee -a "$log"
else printf "\nStep 7: identify imprinted genes from endosperm nuclei\n" | tee -a "$log"; ts=$(date +%s)

	# Identify imprinted genes from mat and pat count matrices
	# -------------------
	# get CxV and VxC endosperm count matrices across all clusters and across each of the 14 CxV and 11 VxC clusters
	mkdir -p "$outdir/imprinting_analysis/matrices"
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/CxV_all_counts.txt $outdir/imprinting_analysis/matrices/CxV_all_counts_C1.txt --cols $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C1.txt 
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/CxV_Col_counts.txt $outdir/imprinting_analysis/matrices/CxV_Col_counts_C1.txt --cols $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C1.txt 
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/CxV_Cvi_counts.txt $outdir/imprinting_analysis/matrices/CxV_Cvi_counts_C1.txt --cols $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C1.txt

	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/VxC_all_counts.txt $outdir/imprinting_analysis/matrices/VxC_all_counts_C1.txt --cols $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C1.txt 
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/VxC_Col_counts.txt $outdir/imprinting_analysis/matrices/VxC_Col_counts_C1.txt --cols $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C1.txt 
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/VxC_Cvi_counts.txt $outdir/imprinting_analysis/matrices/VxC_Cvi_counts_C1.txt --cols $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C1.txt 

	for cc in 1 2 3 4 5 6 7 8 9 10 11 12 13 14; do
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/CxV_all_counts.txt $outdir/imprinting_analysis/matrices/CxV_all_counts_C${cc}.txt --cols $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${cc}.txt 
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/CxV_Col_counts.txt $outdir/imprinting_analysis/matrices/CxV_Col_counts_C${cc}.txt --cols $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${cc}.txt 
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/CxV_Cvi_counts.txt $outdir/imprinting_analysis/matrices/CxV_Cvi_counts_C${cc}.txt --cols $outdir/DE_analysis/cluslists/CxV_4DAP_endo_C${cc}.txt
	done

	for cc in 1 2 3 4 5 6 7 8 9 10 11; do
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/VxC_all_counts.txt $outdir/imprinting_analysis/matrices/VxC_all_counts_C${cc}.txt --cols $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${cc}.txt 
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/VxC_Col_counts.txt $outdir/imprinting_analysis/matrices/VxC_Col_counts_C${cc}.txt --cols $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${cc}.txt 
	$path_to_scripts/subset_large_file.R $outdir/imprinting_analysis/matrices/VxC_Cvi_counts.txt $outdir/imprinting_analysis/matrices/VxC_Cvi_counts_C${cc}.txt --cols $outdir/DE_analysis/cluslists/VxC_4DAP_endo_C${cc}.txt 
	done
	
	# run the ASE calling script
	rm -rf "$lsf/ASE"; mkdir -p "$lsf/ASE"
	rm -rf $outdir/imprinting_analysis/ASE_CxV; rm -rf $outdir/imprinting_analysis/ASE_VxC
	mkdir $outdir/imprinting_analysis/ASE_CxV; mkdir $outdir/imprinting_analysis/ASE_VxC
	cmd1="$path_to_scripts/single_cell_ASE_analysis.R $outdir/imprinting_analysis/matrices/CxV_all_counts.txt $outdir/imprinting_analysis/matrices/CxV_Col_counts.txt $outdir/imprinting_analysis/matrices/CxV_Cvi_counts.txt $outdir/imprinting_analysis/ASE_CxV/ASE_CxV --nameA Col --nameB Cvi"
	cmd2="$path_to_scripts/single_cell_ASE_analysis.R $outdir/imprinting_analysis/matrices/VxC_all_counts.txt $outdir/imprinting_analysis/matrices/VxC_Cvi_counts.txt $outdir/imprinting_analysis/matrices/VxC_Col_counts.txt $outdir/imprinting_analysis/ASE_VxC/ASE_VxC --nameA Cvi --nameB Col"
	bsub -o "$lsf/ASE/ASE_CxV_log.txt" -K "$cmd1" & pid1=$!
	bsub -o "$lsf/ASE/ASE_VxC_log.txt" -K "$cmd2" & pid2=$!

	wait ${pid1} || err_msg "single_cell_ASE_analysis failed for full CxV dataset, see $lsf/ASE/ASE_CxV_log.txt" "$log"
	wait ${pid2} || err_msg "single_cell_ASE_analysis failed for full VxC dataset, see $lsf/ASE/ASE_VxC_log.txt" "$log"

	# classify each gene into overall category based on CxV and VxC results
	
	# merge together the CxV and VxC classification results
	sed '1s/.*/locus_name\tCxV_status\tCxV_btype/' $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_classified.txt > $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_classified_fxheader.txt
	sed '1s/.*/locus_name\tVxC_status\tVxC_btype/' $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_classified.txt > $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_classified_fxheader.txt
	
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_classified_fxheader.txt $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_classified_fxheader.txt locus_name "-" --tokeep all | awk -F$'\t' '{OFS=FS} $2 == "" {$2 = "no data"; $3 = "N/A"} $4 == "" {$4 = "no data"; $5 = "N/A"} {print $0}' > $outdir/imprinting_analysis/final_status.txt

	# add encoding of statuses (1-8 for all 8 possible states in each cross):
	# 1 = no data, 2 = no bias, 3 = weak mat bias, 4 = mat bias, 5 = strong mat bias, 6 = weak pat bias, 7 = pat bias, 8 = strong pat bias
	awk -F$'\t' '{OFS=FS} {
		if (NR == 1) {print $0, "CxV_status_code"}
		else if ($2 == "no data") {print $0,"1"}
		else if ($2 == "no bias" || $2 == "pot. mat bias" || $2 == "pot. pat bias") {print $0,"2"}
		else if ($2 == "weak mat bias") {print $0,"3"}
		else if ($2 == "mat bias") {print $0,"4"}
		else if ($2 == "strong mat bias") {print $0,"5"}
		else if ($2 == "weak pat bias") {print $0,"6"}
		else if ($2 == "pat bias") {print $0,"7"}
		else if ($2 == "strong pat bias") {print $0,"8"}		
		else {print $0,"ERROR"}
		}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/final_status_tmp.txt
	
	awk -F$'\t' '{OFS=FS} {
		if (NR == 1) {print $0, "VxC_status_code"}
		else if ($4 == "no data") {print $0,"1"}
		else if ($4 == "no bias" || $4 == "pot. mat bias" || $4 == "pot. pat bias") {print $0,"2"}
		else if ($4 == "weak mat bias") {print $0,"3"}
		else if ($4 == "mat bias") {print $0,"4"}
		else if ($4 == "strong mat bias") {print $0,"5"}
		else if ($4 == "weak pat bias") {print $0,"6"}
		else if ($4 == "pat bias") {print $0,"7"}
		else if ($4 == "strong pat bias") {print $0,"8"}		
		else {print $0,"ERROR"}
		}' $outdir/imprinting_analysis/final_status_tmp.txt > $outdir/imprinting_analysis/final_status.txt
	
	# get final overall status (accounting for both CxV and VxC results) for each gene
	rm $outdir/imprinting_analysis/final_status_tmp.txt
	
	# any gene with 1 in either cross -> '1' (can't determine imprinting status w/o that info)
	awk -F$'\t' '{OFS=FS} {if (NR==1) {print $0,"final_status"} else if ($6 == "1" || $7 == "1") {print $0,"no data"} else {print $0,""}}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt

	# any gene with 2 in at least one direction and 2 or potential or weak bias in the other -> '2'
	awk -F$'\t' '{OFS=FS} $8 == "" && (($6 == "2" && ($7 == "2" || $7 == "3" || $7 == "6")) || ($7 == "2" && ($6 == "3" || $6 == "6"))) {$8 = "no bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt

	# Imprinted genes:
	# any gene with weak mat. bias in one cross direction and either weak/med/strong in the other -> 'weak MEG'
	awk -F$'\t' '{OFS=FS} $8 == "" && (($6 == "3" && ($7 == "3" || $7 == "4" || $7 == "5")) || ($7 == "3" && ($6 == "4" || $6 == "5"))) {$8 = "weak MEG"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt

	# any gene with mat. bias in one cross direction and either med/strong in the other -> 'MEG'
	awk -F$'\t' '{OFS=FS} $8 == "" && (($6 == "4" && ($7 == "4" || $7 == "5")) || ($7 == "4" && $6 == "5")) {$8 = "MEG"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# any gene with strong mat bias in both crosses -> 'strong MEG'
	awk -F$'\t' '{OFS=FS} $8 == "" && ($6 == "5" && $7 == "5") {$8 = "strong MEG"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt

	# any gene with weak pat. bias in one cross direction and either weak/med/strong in the other -> 'weak PEG'
	awk -F$'\t' '{OFS=FS} $8 == "" && (($6 == "6" && ($7 == "6" || $7 == "7" || $7 == "8")) || ($7 == "6" && ($6 == "7" || $6 == "8"))) {$8 = "weak PEG"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt

	# any gene with pat. bias in one cross direction and either med/strong in the other -> 'PEG'
	awk -F$'\t' '{OFS=FS} $8 == "" && (($6 == "7" && ($7 == "7" || $7 == "8")) || ($7 == "7" && $6 == "8")) {$8 = "PEG"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# any gene with strong pat bias in both crosses -> 'strong PEG'
	awk -F$'\t' '{OFS=FS} $8 == "" && ($6 == "8" && $7 == "8") {$8 = "strong PEG"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# Strain-biased genes:
	# any gene with weak Col bias in one cross direction and either weak/med/strong in the other -> 'weak Col bias'
	awk -F$'\t' '{OFS=FS} $8 == "" && (($6 == "3" && ($7 == "6" || $7 == "7" || $7 == "8")) || ($7 == "6" && ($6 == "3" || $6 == "4" || $6 == "5"))) {$8 = "weak Col bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt

	# any gene with Col bias in one cross direction and either med/strong in the other -> 'Col bias'
	awk -F$'\t' '{OFS=FS} $8 == "" && (($6 == "4" && ($7 == "7" || $7 == "8")) || ($6 == "5" && $7 == "7")) {$8 = "Col bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# any gene with strong Col bias in both crosses -> 'strong Col bias'
	awk -F$'\t' '{OFS=FS} $8 == "" && ($6 == "5" && $7 == "8") {$8 = "strong Col bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt

	# any gene with weak Cvi bias in one cross direction and either weak/med/strong in the other -> 'weak Cvi bias'
	awk -F$'\t' '{OFS=FS} $8 == "" && (($7 == "3" && ($6 == "6" || $6 == "7" || $6 == "8")) || ($6 == "6" && ($7 == "3" || $7 == "4" || $7 == "5"))) {$8 = "weak Cvi bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt

	# any gene with Cvi bias in one cross direction and either med/strong in the other -> 'Cvi bias'
	awk -F$'\t' '{OFS=FS} $8 == "" && (($7 == "4" && ($6 == "7" || $6 == "8")) || ($7 == "5" && $6 == "7")) {$8 = "Cvi bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# any gene with strong Cvi bias in both crosses -> 'strong Cvi bias'
	awk -F$'\t' '{OFS=FS} $8 == "" && ($7 == "5" && $6 == "8") {$8 = "strong Cvi bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# Genes that may be exhibiting strain-specific imprinting
	# any gene with 4 or 5 in CxV but 2 in VxC -> CxV Col bias
	awk -F$'\t' '{OFS=FS} $8 == "" && (($6 == "4" || $6 == "5") && $7 == "2") {$8 = "CxV Col bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# any gene with 7 or 8 in CxV but 2 in VxC -> CxV Cvi bias
	awk -F$'\t' '{OFS=FS} $8 == "" && (($6 == "7" || $6 == "8") && $7 == "2") {$8 = "CxV Cvi bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# any gene with 4 or 5 in VxC but 2 in VxC -> VxC Col bias
	awk -F$'\t' '{OFS=FS} $8 == "" && (($7 == "4" || $7 == "5") && $6 == "2") {$8 = "VxC Cvi bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# any gene with 7 or 8 in VxC but 2 in VxC -> VxC Cvi bias
	awk -F$'\t' '{OFS=FS} $8 == "" && (($7 == "7" || $7 == "8") && $6 == "2") {$8 = "VxC Col bias"} {print $0}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/tmp.txt; mv $outdir/imprinting_analysis/tmp.txt $outdir/imprinting_analysis/final_status.txt
	
	# print out the entire crosstab matrix	
	echo "Summary of imprinting results in CxV (rows) and VxC (columns):"
	echo "	no data	no bias	weak mat bias	mat bias	strong mat bias	weak pat bias	pat bias	strong pat bias"
	statuses=( "" "no data" "no bias" "weak mat bias" "mat bias" "strong mat bias" "weak pat bias" "pat bias" "strong pat bias" )
	for i in 1 2 3 4 5 6 7 8; do
		printf "${statuses[i]}\t"
		cut -f6,7 $outdir/imprinting_analysis/final_status.txt | awk -F$'\t' -v i="$i" '$1 == i {print $2}' | sort | uniq -c | sed 's/^ *//g' | cut -f1 -d' ' | tr '\n' '\t'
		echo ""
	done
	
	echo ""
	echo "Final classifications:"
	cut -f 8 $outdir/imprinting_analysis/final_status.txt | sort | uniq -c
	cut -f1,8  $outdir/imprinting_analysis/final_status.txt >  $outdir/imprinting_analysis/final_status_only.txt


	# Make a few different example plots of imprinted, nonimprinted genes
	# -------------------
	
	mkdir -p "$outdir/imprinting_analysis/plots"
	
	# merge the normalized m/pcount files from CxV and VxC
	sed '1s/^/locus_name/' $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_mcounts_norm.txt > $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_mcounts_norm_fxheader.txt
	sed '1s/^/locus_name/' $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_pcounts_norm.txt > $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_pcounts_norm_fxheader.txt
	sed '1s/^/locus_name/' $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_mcounts_norm.txt > $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_mcounts_norm_fxheader.txt
	sed '1s/^/locus_name/' $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_pcounts_norm.txt > $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_pcounts_norm_fxheader.txt
		
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_mcounts_norm_fxheader.txt $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_mcounts_norm_fxheader.txt locus_name $outdir/imprinting_analysis/all_CV_mcounts_norm.txt --tokeep all
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_pcounts_norm_fxheader.txt $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_pcounts_norm_fxheader.txt locus_name $outdir/imprinting_analysis/all_CV_pcounts_norm.txt --tokeep all
	
	# get list of all CxV and all VxC endosperm samples, make into a single samplefile
	cat <( tail -n+2 $outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt | sed 's/\t.*$/\tCxV/g' ) <( tail -n+2 $outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt | sed 's/\t.*$/\tVxC/g' ) > $outdir/imprinting_analysis/all_CV_nuclist.txt

	# (1) 'nuc' plot of a few example genes, showing allelic bias in each nucleus with detected allelic reads (3 MEGs, 3 PEGs, 3 non-imprinted genes)
	$path_to_scripts/single_cell_RNAseq_plots.R nuc $outdir/imprinting_analysis/plots/nuc_examples1 --mcounts $outdir/imprinting_analysis/all_CV_mcounts_norm.txt --pcounts $outdir/imprinting_analysis/all_CV_pcounts_norm.txt --sampfile $outdir/imprinting_analysis/all_CV_nuclist.txt --genes AT2G35670,AT5G35490,AT1G09540,AT5G39550,AT5G26210,AT4G13460,AT5G63370,AT5G67380,AT5G66320

	# (2) make dotplot of log2FC in CxV vs. VxC for figure
	cut -f1,8 $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_results.txt | sed '1s/.*/locus_name\tCxV_log2fc/' > $outdir/imprinting_analysis/plots/CxV_log2fc.txt
	cut -f1,8 $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_results.txt | sed '1s/.*/locus_name\tVxC_log2fc/' > $outdir/imprinting_analysis/plots/VxC_log2fc.txt

	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/plots/CxV_log2fc.txt $outdir/imprinting_analysis/plots/VxC_log2fc.txt locus_name $outdir/imprinting_analysis/plots/scatterplotdata.txt --tokeep all
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/plots/scatterplotdata.txt $outdir/imprinting_analysis/final_status.txt locus_name $outdir/imprinting_analysis/plots/scatterplotdata.txt --tokeep all

	# add color and pointsize info for plot
	awk -F$'\t' '{OFS=FS} {
		if ($10 == "no bias") {print $0,"gray60","0.15","0"}
		else if ($10 == "weak MEG") {print $0,"\"#F4B5A6\"","0.5","11"}
		else if ($10 == "MEG") {print $0,"\"#DF796A\"","0.5","12"}
		else if ($10 == "strong MEG") {print $0,"\"#C00000\"","0.5","13"}
		else if ($10 == "weak PEG") {print $0,"\"#BCD0F5\"","0.5","14"}
		else if ($10 == "PEG") {print $0,"\"#6C87BE\"","0.5","15"}
		else if ($10 == "strong PEG") {print $0,"\"#2F5597\"","0.5","16"}
		else if ($10 == "weak Col bias") {print $0,"\"#FFE699\"","0.5","5"}
		else if ($10 == "Col bias") {print $0,"\"#F6CE64\"","0.5","6"}
		else if ($10 == "strong Col bias") {print $0,"\"#B89230\"","0.5","7"}
		else if ($10 == "weak Cvi bias") {print $0,"\"#C2E3B8\"","0.5","8"}
		else if ($10 == "Cvi bias") {print $0,"\"#92C964\"","0.5","9"}
		else if ($10 == "strong Cvi bias") {print $0,"\"#5E813F\"","0.5","10"}

		else if ($10 == "CxV Col bias") {print $0,"gray60","0.15","1"}
		else if ($10 == "CxV Cvi bias") {print $0,"gray60","0.15","2"}

		else if ($10 == "VxC Col bias") {print $0,"gray60","0.15","3"}
		else if ($10 == "VxC Cvi bias") {print $0,"gray60","0.15","4"}
		}' $outdir/imprinting_analysis/plots/scatterplotdata.txt | sort -k13n,13 > $outdir/imprinting_analysis/plots/scatterplotdata_tmp.txt
		
	mv $outdir/imprinting_analysis/plots/scatterplotdata_tmp.txt $outdir/imprinting_analysis/plots/scatterplotdata.txt 
		
	$path_to_scripts/scatterplot.R $outdir/imprinting_analysis/plots/scatterplotdata.txt $outdir/imprinting_analysis/plots/scatterplot_log2fc.pdf --xval 2 --yval 3 --color 11 --size 12 --xline 1.44552118384104,1.94552118384104,3.44552118384104,6.44552118384104,0.94552118384104,-0.55447881615896,-3.55447881615896 --yline 1.31535544096388,1.81535544096388,3.31535544096388,6.31535544096388,0.81535544096388,-0.68464455903612,-3.68464455903612 --PDF

	# (3) pie charts comparing Pignatta et al. results, other studies, to our results (Fig. S12)
	
	# downloaded Daniela's data from 2014 eLife paper (Col/Cvi crosses only):
	# elife-03198-fig1-data2-v4.xlsx -> saved to /lab/solexa_gehring/colette/single_cell_seq/pignatta_2014_CV_data/Col_Cvi_rep1.txt, rep2.txt and rep3.txt
	# for each rep, get 'final status'
	mkdir -p "$outdir/imprinting_analysis/pignatta_status" 
	
	# encodes 'status' of gene in Pignatta data:
	# 0 = no data
	# 1 = failed pval
	# 2 = failed IF (matbias)
	# 3 = failed pmat (matbias)
	# 4 = MEG
	# -2 = failed IF (patbias)
	# -3 = failed pmat (patbias)
	# -4 = PEG
	process_pignatta () {
		awk -F$'\t' '{OFS=FS} {
			if (NR == 1) { status = "status" }
			else if ($7 < 0.01) {
				if ($8 >= 2) {
					pmat1 = $2 / ($2+$3)
					pmat2 = $5 / ($4+$5)
					if ($6 == "mother") {
						if (pmat1 >= 0.85 && pmat2 >= 0.85) {
							status = "4"
						} else {
							status = "3"
						}
					} else if ($6 == "father") {
						if (pmat1 <= 0.5 && pmat2 <= 0.5) {
							status = "-4"
						} else {
							status = "-3"
						}
					}
				} else {
					if ($6 == "mother") { status = "2" } else { status = "-2" }
				}
			} else {
				if ($2 + $3 + $4 + $5 < 15) { status = "0" }
				else { status = "1" }
			}} {print $0,status}' $1 > $2
	}
	
	process_pignatta /lab/solexa_gehring/colette/single_cell_seq/pignatta_2014_CV_data/Col_Cvi_rep1.txt $outdir/imprinting_analysis/pignatta_status/Col_Cvi_rep1.txt
	process_pignatta /lab/solexa_gehring/colette/single_cell_seq/pignatta_2014_CV_data/Col_Cvi_rep2.txt $outdir/imprinting_analysis/pignatta_status/Col_Cvi_rep2.txt
	process_pignatta /lab/solexa_gehring/colette/single_cell_seq/pignatta_2014_CV_data/Col_Cvi_rep3.txt $outdir/imprinting_analysis/pignatta_status/Col_Cvi_rep3.txt

	# also process the Belmonte filter data
	awk -F$'\t' '{OFS=FS} {if (NR == 1) {print $1,"belmonte_status"} else if ($5 == "no_data") {print $1, "belmonte_passed"} else if ($5 > 1) {print $1,"belmonte_failed"} else {print $1,"belmonte_passed"}}' /lab/solexa_gehring/colette/single_cell_seq/pignatta_2014_CV_data/Belmonte_filter_data.txt > $outdir/imprinting_analysis/pignatta_status/Belmonte_filter.txt

	# combine all three reps + belmonte
	$path_to_scripts/merge_by_column.R <( cut -f1,12 $outdir/imprinting_analysis/pignatta_status/Col_Cvi_rep1.txt | sed '1s/.*/locus_name\trep1_status/' ) <( cut -f1,12 $outdir/imprinting_analysis/pignatta_status/Col_Cvi_rep2.txt | sed '1s/.*/locus_name\trep2_status/' ) locus_name $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt --tokeep all
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt <( cut -f1,12 $outdir/imprinting_analysis/pignatta_status/Col_Cvi_rep3.txt | sed '1s/.*/locus_name\trep3_status/' ) locus_name $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt --tokeep all
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt $outdir/imprinting_analysis/pignatta_status/Belmonte_filter.txt locus_name $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt --tokeep all
	
	# status of the gene is the middle sorted value of the three replicates
	# add status == 6 for seedcoat filtered genes
	awk -F$'\t' '{OFS=FS} {
		if (NR == 1) {print $0, "final_status"}
		else {
			a[1] = $2; a[2] = $3; a[3] = $4
			asort(a)
			median = a[2]
			if ($5 == "belmonte_failed" && median >= 2) {print $0,"6"} else {print $0,median}
		}
	}' $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt > $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses_tmp.txt
	
	mv $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses_tmp.txt $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt

	# OK, now for all weak/med./strong MEGs/PEGs, get status in Pignatta et al.
	for ss in "weak MEG" "MEG" "strong MEG" "weak PEG" "PEG" "strong PEG"; do		
		tt=$( echo "$ss" | tr ' ' '_' )
		awk -F$'\t' -v s="$ss" 'NR == 1 {print "locus_name"} $8 == s {print $1}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/final_list_${tt}s.txt
		tot=$( cat $outdir/imprinting_analysis/final_list_${tt}s.txt | wc -l ); tot=$(( $tot - 1 ))
		echo "Total number of ${ss}s: $tot"
		
		$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/final_list_${tt}s.txt $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt locus_name "-" --tokeep allx | awk -F$'\t' '{OFS=FS} $6 == "" {$6 = "0"} {print $0}' > $outdir/imprinting_analysis/plots/pignatta_status_${tt}s.txt

		cut -f6 $outdir/imprinting_analysis/plots/pignatta_status_${tt}s.txt | tail -n+2 | sort | uniq -c | sed 's/^ *//g' | tr ' ' '\t' | sort -k2n,2 > $outdir/imprinting_analysis/plots/pignatta_status_${tt}s_summary.txt
		$path_to_scripts/piechart.R $outdir/imprinting_analysis/plots/pignatta_status_${tt}s_summary.txt $outdir/imprinting_analysis/plots/pignatta_status_${tt}s.pdf --PDF --val 1 --labels 2	
	done
	
	# repeat for the endo-specific MEGs/PEGs (actually obtained in step 8 but hey)
	for ss in "weak MEG" "MEG" "strong MEG"; do
		tt=$( echo "$ss" | tr ' ' '_' )
		awk -F$'\t' -v s="$ss" 'NR == 1 {print "locus_name"} $2 == s {print $1}' $outdir/expression_variability/genelists/endo_MEGs_v2.txt > $outdir/imprinting_analysis/final_list_endo_specific_${tt}s.txt
		tot=$( cat $outdir/imprinting_analysis/final_list_endo_specific_${tt}s.txt | wc -l ); tot=$(( $tot - 1 ))
		echo "Total number of endo specific ${ss}s: $tot"
		
		$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/final_list_endo_specific_${tt}s.txt $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt locus_name "-" --tokeep allx | awk -F$'\t' '{OFS=FS} $6 == "" {$6 = "0"} {print $0}' > $outdir/imprinting_analysis/plots/pignatta_status_endo_specific_${tt}s.txt

		cut -f6 $outdir/imprinting_analysis/plots/pignatta_status_endo_specific_${tt}s.txt | tail -n+2 | sort | uniq -c | sed 's/^ *//g' | tr ' ' '\t' | sort -k2n,2 > $outdir/imprinting_analysis/plots/pignatta_status_endo_specific_${tt}s_summary.txt
		$path_to_scripts/piechart.R $outdir/imprinting_analysis/plots/pignatta_status_endo_specific_${tt}s_summary.txt $outdir/imprinting_analysis/plots/pignatta_status_endo_specific_${tt}s.pdf --PDF --val 1 --labels 2	
	done	
	
	
	# other direction: for the Pignatta Col-Cvi MEGs/PEGs, get status in my analysis
	awk -F$'\t' '{OFS=FS} NR == 1 || $6 == 4 {print $1}' $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt > $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_MEGs.txt
	awk -F$'\t' '{OFS=FS} NR == 1 || $6 == -4 {print $1}' $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_statuses.txt > $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_PEGs.txt
	
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_MEGs.txt $outdir/imprinting_analysis/final_status.txt locus_name "-" --tokeep allx | awk -F$'\t' '{OFS=FS} $8 == "" {$8 = "no data"} {print $0}' > $outdir/imprinting_analysis/plots/pignatta_MEGs_mystatus.txt
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/pignatta_status/all_Col_Cvi_PEGs.txt $outdir/imprinting_analysis/final_status.txt locus_name "-" --tokeep allx | awk -F$'\t' '{OFS=FS} $8 == "" {$8 = "no data"} {print $0}' > $outdir/imprinting_analysis/plots/pignatta_PEGs_mystatus.txt
	
	# encode the values
	encode_final_status () {
		awk -F$'\t' '{OFS=FS} {
			if (NR == 1) { status = "status" }
			else if ($8 == "no data") { status = 0 }
			else if ($8 == "no bias") { status = 1 }
			else if ($8 == "weak Col bias") { status = 2 }
			else if ($8 == "Col bias") { status = 3 }
			else if ($8 == "strong Col bias") { status = 4 }
			else if ($8 == "weak Cvi bias") { status = 5 }
			else if ($8 == "Cvi bias") { status = 6 }
			else if ($8 == "strong Cvi bias") { status = 7 }
			else if ($8 == "weak MEG") { status = 8 }
			else if ($8 == "MEG") { status = 9 }
			else if ($8 == "strong MEG") { status = 10 }
			else if ($8 == "weak PEG") { status = 11 }
			else if ($8 == "PEG") { status = 12 }
			else if ($8 == "strong PEG") { status = 13 }
			else if ($8 == "CxV Col bias") { status = 14 }
			else if ($8 == "CxV Cvi bias") { status = 15 }
			else if ($8 == "VxC Col bias") { status = 16 }
			else if ($8 == "VxC Cvi bias") { status = 17 }
			else { status = "ERROR" }
			print $0,status
		}' $1 > $2
	}
	
	encode_final_status $outdir/imprinting_analysis/plots/pignatta_MEGs_mystatus.txt $outdir/imprinting_analysis/plots/pignatta_MEGs_mystatus_tmp.txt
	mv $outdir/imprinting_analysis/plots/pignatta_MEGs_mystatus_tmp.txt $outdir/imprinting_analysis/plots/pignatta_MEGs_mystatus.txt
	cut -f9 $outdir/imprinting_analysis/plots/pignatta_MEGs_mystatus.txt | tail -n+2 | sort | uniq -c | sed 's/^ *//g' | tr ' ' '\t' | sort -k2n,2 > $outdir/imprinting_analysis/plots/pignatta_MEGs_mystatus_summary.txt
	$path_to_scripts/piechart.R $outdir/imprinting_analysis/plots/pignatta_MEGs_mystatus_summary.txt $outdir/imprinting_analysis/plots/pignatta_MEGs_mystatus.pdf --PDF --val 1 --labels 2	
	
	encode_final_status $outdir/imprinting_analysis/plots/pignatta_PEGs_mystatus.txt $outdir/imprinting_analysis/plots/pignatta_PEGs_mystatus_tmp.txt
	mv $outdir/imprinting_analysis/plots/pignatta_PEGs_mystatus_tmp.txt $outdir/imprinting_analysis/plots/pignatta_PEGs_mystatus.txt
	cut -f9 $outdir/imprinting_analysis/plots/pignatta_PEGs_mystatus.txt | tail -n+2 | sort | uniq -c | sed 's/^ *//g' | tr ' ' '\t' | sort -k2n,2 > $outdir/imprinting_analysis/plots/pignatta_PEGs_mystatus_summary.txt
	$path_to_scripts/piechart.R $outdir/imprinting_analysis/plots/pignatta_PEGs_mystatus_summary.txt $outdir/imprinting_analysis/plots/pignatta_PEGs_mystatus.pdf --PDF --val 1 --labels 2	
	
	# (4) get summary table of degree of bias of strong/med/weak MEGs and PEGs in CxV and VxC
	cut -f1,8 "$outdir/imprinting_analysis/ASE_CxV/ASE_CxV_results.txt" | sed '1s/.*/locus_name\tCxV_log2fc/' > "$outdir/imprinting_analysis/ASE_CxV/ASE_CxV_log2fc.txt"
	cut -f1,8 "$outdir/imprinting_analysis/ASE_VxC/ASE_VxC_results.txt" | sed '1s/.*/locus_name\tCxV_log2fc/' > "$outdir/imprinting_analysis/ASE_VxC/ASE_VxC_log2fc.txt"
	
	$path_to_scripts/merge_by_column.R "snRNAseq_analysis/imprinting_analysis/final_status.txt" "$outdir/imprinting_analysis/ASE_CxV/ASE_CxV_log2fc.txt" locus_name "snRNAseq_analysis/imprinting_analysis/_summary_log2fc.txt" --tokeep all
	$path_to_scripts/merge_by_column.R "snRNAseq_analysis/imprinting_analysis/_summary_log2fc.txt" "$outdir/imprinting_analysis/ASE_VxC/ASE_VxC_log2fc.txt" locus_name "snRNAseq_analysis/imprinting_analysis/_summary_log2fc.txt" --tokeep all
	
	awk -F$'\t' '{OFS=FS} {
		if (NR != 1) {
			CxV[$8] += $9
			VxC[$8] += $10
			tot[$8] += 1
			allCxV += $9
			allVxC += $10
			alltot += 1
		}} END {
			print "status","tot_genes","CxV_log2fc_avg", "VxC_log2fc_avg"
			for (a in CxV) {
				print a,tot[a],CxV[a]/tot[a],VxC[a]/tot[a]
			}
			print "all",alltot,
		}' "snRNAseq_analysis/imprinting_analysis/_summary_log2fc.txt"
			

	# Run simulations to get a sense of how the model works under various conditions (few cells, low expression etc.), 
	# using values estimated from the full CxV and VxC runs done in testing - see Ext. data Fig. 8
	# -------------------

	# set up values to use for simulations
	mkdir -p "$outdir/imprinting_analysis/simulations"
	echo "Getting params to run"
	log2_ratios=( -4 -1 0 0.7 1.5 2.3 3 4 7 )		# log2(m/p) ratios to simulate, with 1.5 ~ not imprinted (mimicking slight maternal bias dataset-wide, see Fig. S9)
	log2_ratios_str=( neg4 neg1 0 pos0pt7 pos1pt4 pos2pt3 pos2pt8 pos4 pos7 )	# string version of log2 ratios for filenames
	ncells=( 400 200 100 50 25 10 5 )				# number of cells (nuclei) to simulate
	tot_expr=( 0.01 0.05 0.1 0.15 0.25 0.5 0.75 1.0 1.5 2 3.5 5 15 50 )		# target total expression levels
	tot_expr_str=( 0pt01 0pt05 0pt1 0pt15 0pt25 0pt5 0pt75 1pt0 1pt5 2 3pt5 5 15 50 )	# string version of total expression levels for filenames
	sigma="0.05,0.5,0.95"		
	
	# get coefs from the CxV data (see Fig. S9)
	coef_nu=0.85
	coef_mu_ZINB=0.82
	
	# reset if re-running this section
	echo "Removing previous log files..."
	rm -rf $lsf/ASE_simulations; mkdir -p $lsf/ASE_simulations; pids=(); n=0
	echo "Done"

	# simulate reads with each combination of these params 1000x each
	for ((i=0;i<${#log2_ratios[@]};++i)); do
		for ((j=0;j<${#tot_expr[@]};++j)); do
			for ((k=0;k<${#ncells[@]};++k)); do
				echo "Submitting job $n"
				cmd="$path_to_scripts/single_cell_simulate_reads.R $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_fits.txt $outdir/imprinting_analysis/simulations/r${log2_ratios_str[i]}_m${tot_expr_str[j]}_n${ncells[k]}_sims.txt --n_reps 200 --ncells ${ncells[k]} --log2fc ${log2_ratios[i]} --mean_expr ${tot_expr[j]} --sigma_pctl $sigma --coef_nu $coef_nu --coef_mu_ZINB $coef_mu_ZINB"
				bsub -o "$lsf/ASE_simulations/r${log2_ratios_str[i]}_m${tot_expr_str[j]}_n${ncells[k]}_sims_log.txt" -K "$cmd" & pids[n]=$!
				n=$(( $n + 1 ))
			done
		done
	done
	
	# wait for all jobs to finish
	n=0
	for ((i=0;i<${#log2_ratios[@]};++i)); do
		for ((j=0;j<${#tot_expr[@]};++j)); do
			for ((k=0;k<${#ncells[@]};++k)); do
				wait ${pids[n]} || err_msg "single_cell_simulate_reads failed, see log file $lsf/ASE_simulations/r${log2_ratios_str[i]}_m${tot_expr_str[j]}_n${ncells[k]}_sims_log.txt" "$log"
				n=$(( $n + 1 ))
			done
		done
	done
		
	# combine the various datasets (one dataset per cell size)
	rm -f "$outdir/imprinting_analysis/simulations/allsims_"*
	for ((i=0;i<${#log2_ratios[@]};++i)); do
		echo "Combining log2 ratios ${log2_ratios[i]}"
		for ((j=0;j<${#tot_expr[@]};++j)); do
			echo "Combining tot_expr ${tot_expr[j]}"
			for ((k=0;k<${#ncells[@]};++k)); do
				if [ ! -f "$outdir/imprinting_analysis/simulations/allsims_n${ncells[k]}.txt" ]; then
					cp "$outdir/imprinting_analysis/simulations/r${log2_ratios_str[i]}_m${tot_expr_str[j]}_n${ncells[k]}_sims.txt" "$outdir/imprinting_analysis/simulations/allsims_n${ncells[k]}.txt"
				else
					cat "$outdir/imprinting_analysis/simulations/allsims_n${ncells[k]}.txt" <( tail -n+2 "$outdir/imprinting_analysis/simulations/r${log2_ratios_str[i]}_m${tot_expr_str[j]}_n${ncells[k]}_sims.txt" ) >| "$outdir/imprinting_analysis/simulations/tmp.txt"
					mv -f "$outdir/imprinting_analysis/simulations/tmp.txt" "$outdir/imprinting_analysis/simulations/allsims_n${ncells[k]}.txt"				
				fi
			done
		done
	done
	# plots in Ext data Fig. 8 made from these datasets in Stata, see stata_code.txt
fi


# ----------------------
# Step 8: characterize clusters by examining patterns of overall and allele-specific expression over cell cycle & clusters
# ----------------------
if "$skip8"; then printf "\nSkipping step 8 (examining patterns of overall & allelic expression across clusters, cell cycle) by user request...\n" | tee -a "$log"
else printf "\nStep 8: examining patterns of overall & allelic expression across clusters, cell cycle\n" | tee -a "$log"; ts=$(date +%s)

	# (0) expression of LCM data marker genes in different clusters - see Fig. 1c
	# ----------------
	# look at expression of preglobular, globular, and heart-stage marker genes from Schon and Nodine 2017
	# (based on Belmonte et al. microarray data); note g = globular, h = heart, and g_to_h combines markers
	# from both timepoints since our data are mostly from globular-heart stage
	mkdir -p "$outdir/clustering/Belmonte_marker_genes"
		
	stage="g_to_h"
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/clustering/Belmonte_marker_genes/${stage}_marker_genes/SC3_all_fxaxes" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_all_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Belmonte_marker_genes/Belmonte_all_markers_${stage}.txt" --yorder PEN,MCE,CZE,EP,SUS,CZSC,GSC --fillupper 6 --sizeupper 0.4 --dotsize 15 --allowmissinggenes --xorder 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21_1,21_2,21_3,21_4,21_5,22,23
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/clustering/Belmonte_marker_genes/${stage}_marker_genes/CxV_endo_4DAP_fxaxes" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Belmonte_marker_genes/Belmonte_all_markers_${stage}.txt" --yorder PEN,MCE,CZE,EP,SUS,CZSC,GSC --fillupper 6 --sizeupper 0.4 --dotsize 15 --allowmissinggenes
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/clustering/Belmonte_marker_genes/${stage}_marker_genes/CxV_seedcoat_4DAP_fxaxes" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Belmonte_marker_genes/Belmonte_all_markers_${stage}.txt" --yorder PEN,MCE,CZE,EP,SUS,CZSC,GSC --fillupper 6 --sizeupper 0.4 --dotsize 15 --allowmissinggenes
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/clustering/Belmonte_marker_genes/${stage}_marker_genes/VxC_endo_4DAP_fxaxes" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Belmonte_marker_genes/Belmonte_all_markers_${stage}.txt" --yorder PEN,MCE,CZE,EP,SUS,CZSC,GSC --fillupper 6 --sizeupper 0.4 --dotsize 15 --allowmissinggenes
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/clustering/Belmonte_marker_genes/${stage}_marker_genes/VxC_seedcoat_4DAP_fxaxes" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Belmonte_marker_genes/Belmonte_all_markers_${stage}.txt" --yorder PEN,MCE,CZE,EP,SUS,CZSC,GSC --fillupper 6 --sizeupper 0.4 --dotsize 15 --allowmissinggenes
	
	# get list of all genes detected in >= 1 nucleus (as null set of detected genes, for GO-enrichment analyses later):
	awk -F$'\t' '{if(NR==1) {print $0} else { for(i=1; i<=NF;i++) j+=$i; {if(j > 0) {print $0}}; j=0 }}' "$outdir/count_matrices/all_counts_in_genes_dedup_passQC.txt" | cut -f1 > "$outdir/all_detected_genes.txt"


	# (1) where are various groups of genes expressed in our data? - see Fig. S4
	# ----------------
	mkdir -p "$outdir/expression_variability"

	# (1a) get lists of genes, nuclei clusters, factors, etc. to use in analysis
	
	# get lists of factors (wash/no wash, genotype, tissue) to control for in analysis
	mkdir -p "$outdir/expression_variability/factorlists"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 != "embryo" {print $1,$4,$10}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_tissue.txt"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 != "embryo" {print $1,$4,$10,$11}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_tissue_wash.txt"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 != "embryo" {print $1,$4,$6,$10,$11}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_tissue_wash_ploidy.txt"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 != "embryo" {print $1,$11}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/all_4DAP_nuclei_wash.txt"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 != "embryo" {print $1,$10}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/all_4DAP_nuclei_tissue.txt"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 != "embryo" {print $1,$4}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross.txt"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 != "embryo" {print $1,$4,$11}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_wash.txt"

	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 == "endo" {print $1,$4}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross.txt"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 == "endo" {print $1,$4,$11}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross_wash.txt"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 == "endo" {print $1,$4,$6,$11}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross_wash_ploidy.txt"
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 == "endo" {print $1,$11}' $outdir/_summary_stats/singlenuc_passQC_final.txt > "$outdir/expression_variability/factorlists/endo_4DAP_nuclei_wash.txt"

	# get nuclei clusterings to look for variation over:
	mkdir -p "$outdir/expression_variability/nuclei_lists"	

	# CxV and VxC endosperm / seedcoat clusters
	tail -n+2 "$outdir/clustering/SC3_subsets/all_4DAP_endo_nuclei_clusters.txt" > "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_clusters.txt"

	# endo clusters only
	tail -n+2 "$outdir/clustering/SC3_subsets/all_4DAP_nuclei_clusters.txt" > "$outdir/expression_variability/nuclei_lists/all_4DAP_nuclei_clusters.txt"
	
	# cell cycle phase (all tissues)
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 != "embryo" {print $1,$12}' $path_to_files/trajectory_nuclei_all_stats.txt | tail -n+2 > "$outdir/expression_variability/nuclei_lists/all_4DAP_nuclei_cell_cycle.txt"

	# cell cycle phase (endo only)
	awk -F$'\t' '{OFS=FS} $5 == 4 && ($4 == "CxV" || $4 == "VxC") && $10 == "endo" {print $1,$12}' $path_to_files/trajectory_nuclei_all_stats.txt | tail -n+2 > "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt"

	# get various groups of genes to look at:
	listname=(); listlist=()

	# strong/med/weak MEGs/PEGs
	mkdir -p "$outdir/expression_variability/genelists"
	for ff in "weak MEG" "MEG" "strong MEG" "weak PEG" "PEG" "strong PEG"; do	
		touse=$( echo "$ff" | sed 's/ /_/g' )
		awk -F$'\t' -v a="$ff" '{OFS=FS} $2 == a {print $1}' "$outdir/imprinting_analysis/final_status_only.txt" > "$outdir/expression_variability/genelists/${touse}_list.txt"
		listname+=( "$ff" ); listlist+=( "$outdir/expression_variability/genelists/${touse}_list.txt" )
	done

	# all MEGs
	cat <( awk -F$'\t' '{OFS=FS} {print $1,"strong MEG"}' "$outdir/expression_variability/genelists/strong_MEG_list.txt" ) <( awk -F$'\t' '{OFS=FS} {print $1,"MEG"}' "$outdir/expression_variability/genelists/MEG_list.txt" ) <( awk -F$'\t' '{OFS=FS} {print $1,"weak MEG"}' "$outdir/expression_variability/genelists/weak_MEG_list.txt" ) > "$outdir/expression_variability/genelists/all_MEG_list.txt"
	listname+=( "all MEG" ); listlist+=( "$outdir/expression_variability/genelists/all_MEG_list.txt" )

	# all PEGs
	cat <( awk -F$'\t' '{OFS=FS} {print $1,"strong PEG"}' "$outdir/expression_variability/genelists/strong_PEG_list.txt" ) <( awk -F$'\t' '{OFS=FS} {print $1,"PEG"}' "$outdir/expression_variability/genelists/PEG_list.txt" ) <( awk -F$'\t' '{OFS=FS} {print $1,"weak PEG"}' "$outdir/expression_variability/genelists/weak_PEG_list.txt" ) > "$outdir/expression_variability/genelists/all_PEG_list.txt"
	listname+=( "all PEG" ); listlist+=( "$outdir/expression_variability/genelists/all_PEG_list.txt" )

	# all DE genes across CxV, VxC clusters
	cat "$outdir/DE_analysis/summaries/genes_pass_cluster_nobatch.txt" >  "$outdir/expression_variability/genelists/genes_pass_cluster_nobatch.txt"
	listname+=( "cluster DEG" ); listlist+=( "$outdir/expression_variability/genelists/genes_pass_cluster_nobatch.txt" )

	# all genes that varied significantly by cell cycle
	cut -f2 "$outdir/cell_cycle/traj_analysis_significant_genes.txt" | tail -n+2 > "$outdir/expression_variability/genelists/CC_sig_genes.txt"
	listname+=( "CC DEG" ); listlist+=( "$outdir/expression_variability/genelists/CC_sig_genes.txt" )

	# all 'no bias' genes that were significantly mat/pat biased in at least one cluster
	rm -f "$outdir/expression_variability/genelists/sig_in_one_clus.txt"
	for cc in 1 2 3 4 5 6 7 8 9 10 11 12 13 14; do
		awk -F$'\t' '$2 == "weak mat bias" || $2 == "mat bias" || $2 == "strong mat bias" || $2 == "weak pat bias" || $2 == "pat bias" || $2 == "strong pat bias" {print $1}' $outdir/imprinting_analysis/ASE_CxV_C${cc}/ASE_CxV_C${cc}_classified.txt >> "$outdir/expression_variability/genelists/sig_in_one_clus.txt"
	done
	for cc in 1 2 3 4 5 6 7 8 9 10 11; do
		awk -F$'\t' '$2 == "weak mat bias" || $2 == "mat bias" || $2 == "strong mat bias" || $2 == "weak pat bias" || $2 == "pat bias" || $2 == "strong pat bias" {print $1}' $outdir/imprinting_analysis/ASE_VxC_C${cc}/ASE_VxC_C${cc}_classified.txt >> "$outdir/expression_variability/genelists/sig_in_one_clus.txt"
	done
	echo "locus_name" > "$outdir/expression_variability/genelists/tmp.txt"
	sort "$outdir/expression_variability/genelists/sig_in_one_clus.txt" | uniq >> "$outdir/expression_variability/genelists/tmp.txt"
	mv "$outdir/expression_variability/genelists/tmp.txt" "$outdir/expression_variability/genelists/sig_in_one_clus.txt"
	awk -F$'\t' '{OFS=FS} NR==1 {print $1,"flagged"} NR!=1 {print $1,"Y"}' "$outdir/expression_variability/genelists/sig_in_one_clus.txt" > "$outdir/expression_variability/genelists/tmp.txt"
	$path_to_scripts/merge_by_column.R "$outdir/expression_variability/genelists/tmp.txt" "$outdir/imprinting_analysis/final_status_only.txt" locus_name "-" --tokeep all | awk -F$'\t' '{OFS=FS} $2 == "Y" || $3 ~ /EG/ {if ($3 ~ /EG/ || $3 == "no bias") {print $1,$3} else {print $1,"other"}}' > "$outdir/expression_variability/genelists/sig_in_one_clus.txt"
	awk -F$'\t' '{OFS=FS} $2 == "no bias" {print $1}' "$outdir/expression_variability/genelists/sig_in_one_clus.txt" > "$outdir/expression_variability/genelists/sig_in_one_clus_no_bias.txt"
	listname+=( "sig no bias" ); listlist+=( "$outdir/expression_variability/genelists/sig_in_one_clus_no_bias.txt" )

	# random subset of 500 'no bias' genes that were not significant in at least one cluster
	$path_to_scripts/merge_by_column.R "$outdir/expression_variability/genelists/tmp.txt" "$outdir/imprinting_analysis/final_status_only.txt" locus_name "-" --tokeep all | awk -F$'\t' '{OFS=FS} $2 != "Y" && $3 == "no bias" {print $1}' > "$outdir/expression_variability/genelists/no_bias_all.txt"
	sort -R "$outdir/expression_variability/genelists/no_bias_all.txt" | head -500 > "$outdir/expression_variability/genelists/no_bias_subset.txt"
	listname+=( "no bias subset" ); listlist+=( "$outdir/expression_variability/genelists/no_bias_subset.txt" )
	
	# endo-specific MEGs (based on an initial run of the analysis on all MEGs)
	factype="all_nuc_allfac"; clustype="expr_var_over_factors"
	factors="$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_tissue_wash_ploidy.txt"
	genelist="$outdir/expression_variability/genelists/all_MEG_list.txt"
	touse="all_MEG"
	mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
	$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --genefile "$genelist" --factors "$factors" --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 60 --colorsZ "#3D42B6,#4F74B0,#80ABCD,#B4D8E7,#E4F2F7,white,#FFFEC6,#F9E19B,#F2B26E,#E3754F,#D1382C" > "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt"
	
	echo "locus_name" > $outdir/expression_variability/genelists/tmp.txt
	awk -F$'\t' '{OFS=FS} NR!=1 && $2 > 0 {print $1}' $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier_factor3_zscores.txt >> $outdir/expression_variability/genelists/tmp.txt
	sed -e '1i\locus_name\ttype' $outdir/expression_variability/genelists/all_MEG_list.txt > $outdir/expression_variability/genelists/tmp2.txt
	$path_to_scripts/merge_by_column.R $outdir/expression_variability/genelists/tmp.txt $outdir/expression_variability/genelists/tmp2.txt locus_name - --tokeep allx | tail -n+2 | sort -k2,2 > "$outdir/expression_variability/genelists/endo_MEGs.txt"
	rm $outdir/expression_variability/genelists/tmp.txt $outdir/expression_variability/genelists/tmp2.txt
	listname+=( "endo MEG" ); listlist+=( "$outdir/expression_variability/genelists/endo_MEGs.txt" )
	
	factype="all_nuc_threefac"; clustype="expr_var_over_factors"
	factors="$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_tissue_wash.txt"
	genelist="$outdir/expression_variability/genelists/all_MEG_list.txt"
	touse="all_MEG"
	mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
	$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --genefile "$genelist" --factors "$factors" --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 60 --colorsZ "#3D42B6,#4F74B0,#80ABCD,#B4D8E7,#E4F2F7,white,#FFFEC6,#F9E19B,#F2B26E,#E3754F,#D1382C" > "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt"
	
	echo "locus_name" > $outdir/expression_variability/genelists/tmp.txt
	awk -F$'\t' '{OFS=FS} NR!=1 && $2 > 0 {print $1}' $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier_factor2_zscores.txt >> $outdir/expression_variability/genelists/tmp.txt
	sed -e '1i\locus_name\ttype' $outdir/expression_variability/genelists/all_MEG_list.txt > $outdir/expression_variability/genelists/tmp2.txt
	$path_to_scripts/merge_by_column.R $outdir/expression_variability/genelists/tmp.txt $outdir/expression_variability/genelists/tmp2.txt locus_name - --tokeep allx | tail -n+2 | sort -k2,2 > "$outdir/expression_variability/genelists/endo_MEGs_v2.txt"
	rm $outdir/expression_variability/genelists/tmp.txt $outdir/expression_variability/genelists/tmp2.txt
	listname+=( "endo MEG v2" ); listlist+=( "$outdir/expression_variability/genelists/endo_MEGs_v2.txt" )
	
	# genes up or downreg. in endosperm relative to seed coat
	listname+=( "endo gr SC" ); listlist+=( "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/intersect_genes_endo_gr_seedcoat.txt" )
	listname+=( "endo le SC" ); listlist+=( "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/intersect_genes_endo_le_seedcoat.txt" )

	# a list of chromatin genes from Mary
	# get the subset of them that are in our expression matrix
	echo "locus_name" > "$outdir/tmp.txt"
	cat "$path_to_files/chromdb_loci.txt" >> "$outdir/tmp.txt"
	cut -f1 $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt > "$outdir/tmp2.txt"
	
	$path_to_scripts/merge_by_column.R "$outdir/tmp.txt" "$outdir/tmp2.txt" locus_name - --tokeep merged | tail -n+2 > "$outdir/chromdb_loci_touse.txt"
	listname+=( "chromatin related" ); listlist+=( "$outdir/chromdb_loci_touse.txt" )

	# (1b) run analysis on all gene groups, nuclei groups, controlling for different sets of factors
	# ----------------
	pids=(); logfiles=()
	for ((i=0;i<${#listname[@]};++i)); do
		echo "Analyzing expression pattern of genes in list: ${listname[i]}s"
		touse=$( echo "${listname[i]}" | sed 's/ /_/g' )

		# expression patterns in all nuclei, over nuclei clusters, controlling for cross + tissue
		factype="all_nuc_control_cross_tissue"; clustype="expr_var_over_clusters"
		cluslist="$outdir/expression_variability/nuclei_lists/all_4DAP_nuclei_clusters.txt"
		factors="$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_tissue.txt"
		mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_kmeans --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --kmeans --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" )
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans_log.txt" )
		if [ "${listname[i]}" = "CC DEG" ] || [ "${listname[i]}" = "cluster DEG" ]; then
			cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_kmeans60 --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 8 --k 60 --kmeans --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --height 16"
			bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans60_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans60_log.txt" )
		fi
		
		# expression patterns in endo nuclei, over nuclei clusters, controlling for cross only
		factype="endo_nuc_control_cross"; clustype="expr_var_over_clusters"
		cluslist="$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_clusters.txt"
		factors="$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross.txt"
		mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_kmeans --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --kmeans --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" )
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans_log.txt" )
		if [ "${listname[i]}" = "CC DEG" ] || [ "${listname[i]}" = "cluster DEG" ]; then
			cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_kmeans60 --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 8 --k 60 --kmeans --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --height 16"
			bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans60_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans60_log.txt" )
		fi		

		# expression patterns in all nuclei, by all factors
		factype="all_nuc"; clustype="expr_var_over_factors"
		factors="$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_tissue_wash_ploidy.txt"
		mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 20 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" )
				
		# expression patterns in all nuclei, by three main factors
		factype="all_nuc_threefac"; clustype="expr_var_over_factors"
		factors="$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_tissue_wash.txt"
		mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 20 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" )
				
		# expression patterns in endo nuclei, by all factors
		factype="endo_nuc"; clustype="expr_var_over_factors"
		factors="$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross_wash_ploidy.txt"
		mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 20 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" )
			
		# expression patterns in all nuclei, over cell cycle, controlling for cross + tissue
		factype="all_nuc_control_cross_tissue"; clustype="expr_var_over_CC"
		cluslist="$outdir/expression_variability/nuclei_lists/all_4DAP_nuclei_cell_cycle.txt"
		factors="$outdir/expression_variability/factorlists/all_4DAP_nuclei_cross_tissue.txt"
		mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
		mkdir -p "$outdir/expression_variability/$clustype/$factype/${touse}_max3"
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_kmeans --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --kmeans --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" )
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans_log.txt" )
		if [ "${listname[i]}" = "CC DEG" ] || [ "${listname[i]}" = "cluster DEG" ]; then
			cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_kmeans60 --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --kmeans --clustersep 8 --k 60 --colorder G0,G1,G1toS,S,G2,M --kmeans --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --height 16"
			bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans60_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans60_log.txt" )
		fi
		# repeat with max Z = 3
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/${touse}_max3/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/${touse}_max3/${touse}_kmeans --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --kmeans --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_hier_log.txt" )
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_kmeans_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_kmeans_log.txt" )
		if [ "${listname[i]}" = "CC DEG" ] || [ "${listname[i]}" = "cluster DEG" ]; then
			cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/${touse}_max3/${touse}_kmeans60 --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --kmeans --clustersep 8 --k 60 --colorder G0,G1,G1toS,S,G2,M --kmeans --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --height 16"
			bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_kmeans60_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_kmeans60_log.txt" )
		fi
		
		# expression patterns in endo nuclei, over cell cycle, controlling for cross only
		factype="endo_nuc_control_cross"; clustype="expr_var_over_CC"
		cluslist="$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt"
		factors="$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross.txt"
		mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
		mkdir -p "$outdir/expression_variability/$clustype/$factype/${touse}_max3"
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_kmeans --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --kmeans --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" )
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans_log.txt" )
		if [ "${listname[i]}" = "CC DEG" ] || [ "${listname[i]}" = "cluster DEG" ]; then
			cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_kmeans60 --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --kmeans --clustersep 8 --k 60 --colorder G0,G1,G1toS,S,G2,M --kmeans --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --height 16"
			bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans60_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans60_log.txt" )
		fi
		# repeat with max Z = 3
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/${touse}_max3/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/${touse}_max3/${touse}_kmeans --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --kmeans --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_hier_log.txt" )
		bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_kmeans_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_kmeans_log.txt" )
		if [ "${listname[i]}" = "CC DEG" ] || [ "${listname[i]}" = "cluster DEG" ]; then
			cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/${touse}_max3/${touse}_kmeans60 --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --kmeans --clustersep 8 --k 60 --colorder G0,G1,G1toS,S,G2,M --kmeans --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --height 16"
			bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_kmeans60_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_kmeans60_log.txt" )
		fi
	done
		
	# wait for all jobs to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "cluster_gene_expression for assessing expression variability failed, see logfile ${logfiles[i]}" "$log"
	done

	# repeat the clustering over cell cycle while omitting VxC E9, which has depleted MEG expression already,
	# to make sure the S-phase effect in cell cycle isn't just because of something weird with that one cluster
	printf "nucID\tccphase\n" >| tmp1.txt
	cat "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" >> tmp1.txt

	printf "nucID\tcluster\n" >| tmp2.txt
	cat "$outdir/expression_variability/nuclei_lists/all_4DAP_nuclei_clusters.txt" >> tmp2.txt

	$path_to_scripts/merge_by_column.R tmp1.txt tmp2.txt nucID - --tokeep merged | awk '$3 != "VxC_endo_9"' | cut -f1,2 | tail -n+2 > $outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle_noVxCE9.txt

	touse="all_MEG"
	factype="endo_nuc_control_cross"; clustype="expr_var_over_CC"
	cluslist="$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle_noVxCE9.txt"
	factors="$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross.txt"
	genefile="$outdir/expression_variability/genelists/all_MEG_list.txt"
	mkdir -p "$outdir/expression_variability/$clustype/$factype/${touse}_max3_noVxCE9_v2"
	$path_to_wipscripts/cluster_gene_expression.R $outdir/expression_variability/$clustype/$factype/${touse}_max3_noVxCE9_v2/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile $genefile --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C


	# (1c) for the clustering with k = 60 of the cluster DEGs, more detailed analyses
	# ----------------
	mkdir -p "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/belmonte_dotplot"
		
	clusfiles=""; clusnames=""
	for ((i=1;i<=60;++i)); do
		awk -v a="$i" '$2 == a {print $1}' "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_heatmap_k60_zscores_roworder.txt" > "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/belmonte_dotplot/DE_k60_cluster${i}.txt"
		clusfiles="${clusfiles},$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/belmonte_dotplot/DE_k60_cluster${i}.txt"
		clusnames="${clusnames},${i}"
	done
	
	clusfiles="${clusfiles:1}"
	clusnames="${clusnames:1}"	
	
	# for each of the K=60 clusters, look at overlap w/ the Belmonte markers (Fig. S4b)
	# using only the Belmonte genes that are also in our background set of 29427 detected genes
	mkdir -p "$outdir/DE_analysis/summaries/belmonte_dotplot"
	printf "locus_name\ttissue\n" > $outdir/tmp.txt
	cat /lab/solexa_gehring/colette/single_cell_seq/Belmonte_marker_genes/Belmonte_all_markers_g_to_h.txt >> $outdir/tmp.txt
	
	$path_to_scripts/merge_by_column.R $outdir/tmp.txt $outdir/all_detected_genes.txt locus_name - --tokeep merged | tail -n+2 >| $outdir/Belmonte_markers_touse_g_to_h.txt
	
	for rr in PEN MCE CZE EP CZSC GSC; do
		awk -F$'\t' -v a="$rr" '{OFS=FS} $2 == a {print $1}' $outdir/Belmonte_markers_touse_g_to_h.txt > $outdir/DE_analysis/summaries/belmonte_dotplot/Belmonte_${rr}_markers_g_to_h_inBG.txt
		belmontefiles="${belmontefiles},$outdir/DE_analysis/summaries/belmonte_dotplot/Belmonte_${rr}_markers_g_to_h_inBG.txt"
		belmontenames="${belmontenames},${rr}"
	done
	
	belmontefiles="${belmontefiles:1}"
	belmontenames="${belmontenames:1}"
			
	$path_to_scripts/gene_overlaps_dotplot.R "$belmontefiles" "$clusfiles" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/belmonte_dotplot/DE_k60_dotplot_fx" --namesA $belmontenames --namesB $clusnames --popsize 29427 --scoretype hgeo --height 20 --dotsize 80
			
	# also get, for every nuclei cluster, a list of all sig. up and down genes (Table S2)
	mkdir -p $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists
	rm $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up.txt
	rm $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn.txt
	
	sed '1s/^/locus_name\t/' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_sig.txt > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_sig_fxheader.txt
	sed '1s/^/locus_name\t/' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_zscores.txt > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_zscores_fxheader.txt
	CxV_endo_up=(); CxV_endo_dwn=()
	VxC_endo_up=(); VxC_endo_dwn=()
	CxV_SC_up=(); CxV_SC_dwn=()
	VxC_SC_up=(); VxC_SC_dwn=()
	
	for ((i=2;i<=40;++i)); do
		sname=$( head -1 $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_sig_fxheader.txt | cut -f${i} )		
		cut -f1,${i} $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_sig_fxheader.txt | awk -F$'\t' -v s="$sname" -v l="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists" '{OFS=FS} $2 == 1 {print $1 > l"/"s"_up.txt"} $2 == -1 {print $1 > l"/"s"_dwn.txt"}' 
		nup=$( cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_up.txt | wc -l )
		ndwn=$( cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_dwn.txt | wc -l )
		echo "$sname	$nup	$ndwn"
		sed "s/^/${sname}\t/g" $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_up.txt >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up.txt
		sed "s/^/${sname}\t/g" $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_dwn.txt >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn.txt
		
		# test no duplicates
		nupdd=$( cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_up.txt | sort -u | wc -l )
		ndwndd=$( cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_dwn.txt | sort -u | wc -l )
		[ "$nup" != "$nupdd" ] && echo "For ${sname} upreg. genes, w/ dup has $nup and dedup has $nupdd"
		[ "$ndwn" != "$ndwndd" ] && echo "For ${sname} dwnreg. genes, w/ ddwn has $ndwn and dedup has $ndwndd"
		
		# save to arrays	
		snum=$( echo "$sname" | cut -f3 -d'_' )
		[[ "$sname" = *"CxV_endo"* ]] && { CxV_endo_up["$snum"]="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_up.txt"; CxV_endo_dwn["$snum"]="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_dwn.txt"; }
		[[ "$sname" = *"VxC_endo"* ]] && { VxC_endo_up["$snum"]="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_up.txt"; VxC_endo_dwn["$snum"]="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_dwn.txt"; }
		[[ "$sname" = *"CxV_seedcoat"* ]] && { CxV_SC_up["$snum"]="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_up.txt"; CxV_SC_dwn["$snum"]="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_dwn.txt"; }
		[[ "$sname" = *"VxC_seedcoat"* ]] && { VxC_SC_up["$snum"]="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_up.txt"; VxC_SC_dwn["$snum"]="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/${sname}_dwn.txt"; }
	done

	# also get lists of unique markers for each cluster (e.g. not present for any other cluster)
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up.txt | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >| $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_up.txt
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn.txt | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >| $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_dwn.txt

	# repeat within genotype & tissue
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up.txt | awk '$1 ~ /CxV_endo/' | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >| $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_up_bygrp.txt
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up.txt | awk '$1 ~ /VxC_endo/' | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_up_bygrp.txt
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up.txt | awk '$1 ~ /CxV_seedcoat/' | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_up_bygrp.txt
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up.txt | awk '$1 ~ /VxC_seedcoat/' | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_up_bygrp.txt

	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn.txt | awk '$1 ~ /CxV_endo/' | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >| $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_dwn_bygrp.txt
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn.txt | awk '$1 ~ /VxC_endo/' | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_dwn_bygrp.txt
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn.txt | awk '$1 ~ /CxV_seedcoat/' | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_dwn_bygrp.txt
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn.txt | awk '$1 ~ /VxC_seedcoat/' | awk -F$'\t' '{OFS=FS} {if ($2 in uniq) {rep[$2]} else {uniq[$2]=$1}} END {for (a in uniq) {if (! (a in rep)) {print a,uniq[a]}}}' | sort -k2,2 -k1,1 >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_uniq_dwn_bygrp.txt
	
	CxV_endo_up_ll=""; CxV_endo_dwn_ll=""; CxV_endo_nn=""
	for ((i=1;i<=${#CxV_endo_up[@]};++i)); do
		CxV_endo_up_ll="${CxV_endo_up_ll},${CxV_endo_up[i]}"
		CxV_endo_dwn_ll="${CxV_endo_dwn_ll},${CxV_endo_dwn[i]}"
		CxV_endo_nn="${CxV_endo_nn},${i}"
	done
	CxV_endo_up_ll="${CxV_endo_up_ll:1}"; CxV_endo_dwn_ll="${CxV_endo_dwn_ll:1}"; CxV_endo_nn="${CxV_endo_nn:1}"
	
	VxC_endo_up_ll=""; VxC_endo_dwn_ll=""; VxC_endo_nn=""
	for ((i=1;i<=${#VxC_endo_up[@]};++i)); do
		VxC_endo_up_ll="${VxC_endo_up_ll},${VxC_endo_up[i]}"
		VxC_endo_dwn_ll="${VxC_endo_dwn_ll},${VxC_endo_dwn[i]}"
		VxC_endo_nn="${VxC_endo_nn},${i}"
	done
	VxC_endo_up_ll="${VxC_endo_up_ll:1}"; VxC_endo_dwn_ll="${VxC_endo_dwn_ll:1}"; VxC_endo_nn="${VxC_endo_nn:1}"
	
	CxV_SC_up_ll=""; CxV_SC_dwn_ll=""; CxV_SC_nn=""
	for ((i=1;i<=${#CxV_SC_up[@]};++i)); do
		CxV_SC_up_ll="${CxV_SC_up_ll},${CxV_SC_up[i]}"
		CxV_SC_dwn_ll="${CxV_SC_dwn_ll},${CxV_SC_dwn[i]}"
		CxV_SC_nn="${CxV_SC_nn},${i}"
	done
	CxV_SC_up_ll="${CxV_SC_up_ll:1}"; CxV_SC_dwn_ll="${CxV_SC_dwn_ll:1}"; CxV_SC_nn="${CxV_SC_nn:1}"
	
	VxC_SC_up_ll=""; VxC_SC_dwn_ll=""; VxC_SC_nn=""
	for ((i=1;i<=${#VxC_SC_up[@]};++i)); do
		VxC_SC_up_ll="${VxC_SC_up_ll},${VxC_SC_up[i]}"
		VxC_SC_dwn_ll="${VxC_SC_dwn_ll},${VxC_SC_dwn[i]}"
		VxC_SC_nn="${VxC_SC_nn},${i}"
	done
	VxC_SC_up_ll="${VxC_SC_up_ll:1}"; VxC_SC_dwn_ll="${VxC_SC_dwn_ll:1}"; VxC_SC_nn="${VxC_SC_nn:1}"

	# dot gene_overlaps_dotplot plots comparing the CxV and VxC endosperm up/dwn genes to each other, and the seedcoat genes (Fig. S5)
	mkdir -p "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots"
	
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_up_ll" "$VxC_endo_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_upreg_CxV_vs_VxC" --namesA $CxV_endo_nn --namesB $VxC_endo_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 7
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_dwn_ll" "$VxC_endo_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_dwnreg_CxV_vs_VxC" --namesA $CxV_endo_nn --namesB $VxC_endo_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 7
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_SC_up_ll" "$VxC_SC_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_upreg_CxV_vs_VxC" --namesA $CxV_SC_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 7 --height 5.1
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_SC_dwn_ll" "$VxC_SC_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_dwnreg_CxV_vs_VxC" --namesA $CxV_SC_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 7 --height 5.1

	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_up_ll" "$VxC_endo_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_up_CxV_vs_dwn_VxC" --namesA $CxV_endo_nn --namesB $VxC_endo_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 7
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_dwn_ll" "$VxC_endo_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_dwn_CxV_vs_up_VxC" --namesA $CxV_endo_nn --namesB $VxC_endo_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 7
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_SC_up_ll" "$VxC_SC_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_up_CxV_vs_dwn_VxC" --namesA $CxV_SC_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 7
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_SC_dwn_ll" "$VxC_SC_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_dwn_CxV_vs_up_VxC" --namesA $CxV_SC_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 7

	# also compare all CxV endo clusters to themselves, etc.
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_up_ll" "$CxV_endo_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_upreg_CxV_vs_CxV" --namesA $CxV_endo_nn --namesB $CxV_endo_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 9
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_endo_up_ll" "$VxC_endo_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_upreg_VxC_vs_VxC" --namesA $VxC_endo_nn --namesB $VxC_endo_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 10 --height 7
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_SC_up_ll" "$CxV_SC_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_upreg_CxV_vs_CxV" --namesA $CxV_SC_nn --namesB $CxV_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 7 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_SC_up_ll" "$VxC_SC_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_upreg_VxC_vs_VxC" --namesA $VxC_SC_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 8.1 --height 5.1
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_dwn_ll" "$CxV_endo_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_dwnreg_CxV_vs_CxV" --namesA $CxV_endo_nn --namesB $CxV_endo_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 9
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_endo_dwn_ll" "$VxC_endo_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_dwnreg_VxC_vs_VxC" --namesA $VxC_endo_nn --namesB $VxC_endo_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 10 --height 7
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_SC_dwn_ll" "$CxV_SC_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_dwnreg_CxV_vs_CxV" --namesA $CxV_SC_nn --namesB $CxV_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 7 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_SC_dwn_ll" "$VxC_SC_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_dwnreg_VxC_vs_VxC" --namesA $VxC_SC_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 8.1 --height 5.1
	
	# also compare the endo and SC samples within a cross
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_up_ll" "$CxV_SC_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/CxV_endo_vs_SC_upreg" --namesA $CxV_endo_nn --namesB $CxV_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_endo_up_ll" "$VxC_SC_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/VxC_endo_vs_SC_upreg" --namesA $VxC_endo_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 10 --height 5.1
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_dwn_ll" "$CxV_SC_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/CxV_endo_vs_SC_dwnreg" --namesA $CxV_endo_nn --namesB $CxV_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_endo_dwn_ll" "$VxC_SC_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/VxC_endo_vs_SC_dwnreg" --namesA $VxC_endo_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 10 --height 5.1

	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_up_ll" "$CxV_SC_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/CxV_endo_up_vs_SC_dwn" --namesA $CxV_endo_nn --namesB $CxV_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_endo_up_ll" "$VxC_SC_dwn_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/VxC_endo_up_vs_SC_dwn" --namesA $VxC_endo_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 10 --height 5.1
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_dwn_ll" "$CxV_SC_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/CxV_endo_dwn_vs_SC_up" --namesA $CxV_endo_nn --namesB $CxV_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 12 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_endo_dwn_ll" "$VxC_SC_up_ll" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/VxC_endo_dwn_vs_SC_up" --namesA $VxC_endo_nn --namesB $VxC_SC_nn --popsize 4500 --scoretype hgeo --sizeupper 100 --width 10 --height 5.1

	# also compare to the Belmonte data
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_up_ll" "$belmontefiles" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_upreg_CxV_vs_Bel" --namesA $CxV_endo_nn --namesB $belmontenames --popsize 29427 --scoretype hgeo --sizeupper 100 --width 12 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_endo_dwn_ll" "$belmontefiles" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_dwnreg_CxV_vs_Bel" --namesA $CxV_endo_nn --namesB $belmontenames --popsize 29427 --scoretype hgeo --sizeupper 100 --width 12 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_SC_up_ll" "$belmontefiles" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_upreg_CxV_vs_Bel" --namesA $CxV_SC_nn --namesB $belmontenames --popsize 29427 --scoretype hgeo --sizeupper 100 --width 7 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$CxV_SC_dwn_ll" "$belmontefiles" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_dwnreg_CxV_vs_Bel" --namesA $CxV_SC_nn --namesB $belmontenames --popsize 29427 --scoretype hgeo --sizeupper 100 --width 7 --height 4
	
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_endo_up_ll" "$belmontefiles" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_upreg_VxC_vs_Bel" --namesA $VxC_endo_nn --namesB $belmontenames --popsize 29427 --scoretype hgeo --sizeupper 100 --width 10 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_endo_dwn_ll" "$belmontefiles" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/endo_dwnreg_VxC_vs_Bel" --namesA $VxC_endo_nn --namesB $belmontenames --popsize 29427 --scoretype hgeo --sizeupper 100 --width 10 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_SC_up_ll" "$belmontefiles" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_upreg_VxC_vs_Bel" --namesA $VxC_SC_nn --namesB $belmontenames --popsize 29427 --scoretype hgeo --sizeupper 100 --width 8.1 --height 4
	$path_to_scripts/gene_overlaps_dotplot.R "$VxC_SC_dwn_ll" "$belmontefiles" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/SC_dwnreg_VxC_vs_Bel" --namesA $VxC_SC_nn --namesB $belmontenames --popsize 29427 --scoretype hgeo --sizeupper 100 --width 8.1 --height 4
	
	# finally, compare overlap of the nuclei in each cluster to the cell cycle data	(Fig. S6)
	mkdir -p "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists"
	
	sed '1s/geneID/name/' $outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt > $outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_fxheader.txt
	$path_to_scripts/merge_by_column.R $outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_fxheader.txt $path_to_files/trajectory_nuclei_all_stats.txt name - --tokeep merged | cut -f1,13 > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_endo_CC.txt
	
	sed '1s/geneID/name/' $outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt > $outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_fxheader.txt
	$path_to_scripts/merge_by_column.R $outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_fxheader.txt $path_to_files/trajectory_nuclei_all_stats.txt name - --tokeep merged | cut -f1,13 > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_endo_CC.txt
	
	sed '1s/geneID/name/' $outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt > $outdir/clustering/SC3_subsets/CxV_SC_4DAP_final_clusters_fxheader.txt
	$path_to_scripts/merge_by_column.R $outdir/clustering/SC3_subsets/CxV_SC_4DAP_final_clusters_fxheader.txt $path_to_files/trajectory_nuclei_all_stats.txt name - --tokeep merged | cut -f1,13 > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_SC_CC.txt
	
	sed '1s/geneID/name/' $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt > $outdir/clustering/SC3_subsets/VxC_SC_4DAP_final_clusters_fxheader.txt
	$path_to_scripts/merge_by_column.R $outdir/clustering/SC3_subsets/VxC_SC_4DAP_final_clusters_fxheader.txt $path_to_files/trajectory_nuclei_all_stats.txt name - --tokeep merged | cut -f1,13 > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_SC_CC.txt
	
	CCnames="G0,G1,G1toS,S,G2,M"
	CCfiles_CxVendo=""
	CCfiles_CxVSC=""
	CCfiles_VxCendo=""
	CCfiles_VxCSC=""
	for cc in G0 G1 G1toS S G2 M; do
		awk -F$'\t' -v a="$cc" '$2 == a' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_endo_CC.txt | cut -f1 > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_endo_CC_${cc}.txt
		awk -F$'\t' -v a="$cc" '$2 == a' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_endo_CC.txt | cut -f1 > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_endo_CC_${cc}.txt
		awk -F$'\t' -v a="$cc" '$2 == a' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_SC_CC.txt | cut -f1 > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_SC_CC_${cc}.txt
		awk -F$'\t' -v a="$cc" '$2 == a' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_SC_CC.txt | cut -f1 > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_SC_CC_${cc}.txt
		CCfiles_CxVendo="${CCfiles_CxVendo},$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_endo_CC_${cc}.txt"
		CCfiles_VxCendo="${CCfiles_VxCendo},$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_endo_CC_${cc}.txt"
		CCfiles_CxVSC="${CCfiles_CxVSC},$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_SC_CC_${cc}.txt"
		CCfiles_VxCSC="${CCfiles_VxCSC},$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_SC_CC_${cc}.txt"
	done
	CCfiles_CxVendo="${CCfiles_CxVendo:1}"
	CCfiles_VxCendo="${CCfiles_VxCendo:1}"
	CCfiles_CxVSC="${CCfiles_CxVSC:1}"
	CCfiles_VxCSC="${CCfiles_VxCSC:1}"
	
	# CxV endo
	clusfileslist=""; clusnames=""
	for ((i=1;i<=14;++i)); do
		clusfileslist="${clusfileslist},${outdir}/clustering/repeat_noCC/gene_overlaps/CxV_endo_4DAP/cluster${i}_wCC.txt"; clusnames="${clusnames},${i}"
	done
	clusfileslist="${clusfileslist:1}"; clusnames="${clusnames:1}"	
	totnuc=$( cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_endo_CC.txt | wc -l )
	totnuc=$(( $totnuc - 1 ))
	$path_to_scripts/gene_overlaps_dotplot.R "$clusfileslist" "$CCfiles_CxVendo" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/CxV_endo_nuc_vs_CC" --namesA $clusnames --namesB $CCnames --popsize $totnuc --scoretype hgeo --width 12 --height 4 --sizeupper 15

	# VxC endo
	clusfileslist=""; clusnames=""
	for ((i=1;i<=11;++i)); do
		clusfileslist="${clusfileslist},${outdir}/clustering/repeat_noCC/gene_overlaps/VxC_endo_4DAP/cluster${i}_wCC.txt"; clusnames="${clusnames},${i}"
	done
	clusfileslist="${clusfileslist:1}"; clusnames="${clusnames:1}"	
	totnuc=$( cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_endo_CC.txt | wc -l )
	totnuc=$(( $totnuc - 1 ))
	$path_to_scripts/gene_overlaps_dotplot.R "$clusfileslist" "$CCfiles_VxCendo" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/VxC_endo_nuc_vs_CC" --namesA $clusnames --namesB $CCnames --popsize $totnuc --scoretype hgeo --width 10 --height 4 --sizeupper 15

	# CxV SC
	mkdir -p $outdir/clustering/repeat_noCC/gene_overlaps/CxV_seedcoat_4DAP
	clusfileslist=""; clusnames=""
	for ((i=1;i<=6;++i)); do
		awk -v a="$i" '$2 == a {print $1}' "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt" > "$outdir/clustering/repeat_noCC/gene_overlaps/CxV_seedcoat_4DAP/cluster${i}_wCC.txt"
		clusfileslist="${clusfileslist},${outdir}/clustering/repeat_noCC/gene_overlaps/CxV_seedcoat_4DAP/cluster${i}_wCC.txt"; clusnames="${clusnames},${i}"
	done
	clusfileslist="${clusfileslist:1}"; clusnames="${clusnames:1}"	
	totnuc=$( cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/CxV_SC_CC.txt | wc -l )
	totnuc=$(( $totnuc - 1 ))
	$path_to_scripts/gene_overlaps_dotplot.R "$clusfileslist" "$CCfiles_CxVSC" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/CxV_SC_nuc_vs_CC" --namesA $clusnames --namesB $CCnames --popsize $totnuc --scoretype hgeo --width 7 --height 4 --sizeupper 15

	# VxC SC
	mkdir -p $outdir/clustering/repeat_noCC/gene_overlaps/VxC_seedcoat_4DAP
	clusfileslist=""; clusnames=""
	for ((i=1;i<=8;++i)); do
		awk -v a="$i" '$2 == a {print $1}' "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt" > "$outdir/clustering/repeat_noCC/gene_overlaps/VxC_seedcoat_4DAP/cluster${i}_wCC.txt"
		clusfileslist="${clusfileslist},${outdir}/clustering/repeat_noCC/gene_overlaps/VxC_seedcoat_4DAP/cluster${i}_wCC.txt"; clusnames="${clusnames},${i}"
	done
	clusfileslist="${clusfileslist:1}"; clusnames="${clusnames:1}"	
	totnuc=$( cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/CC_lists/VxC_SC_CC.txt | wc -l )
	totnuc=$(( $totnuc - 1 ))
	$path_to_scripts/gene_overlaps_dotplot.R "$clusfileslist" "$CCfiles_VxCSC" "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/plots/VxC_SC_nuc_vs_CC" --namesA $clusnames --namesB $CCnames --popsize $totnuc --scoretype hgeo --width 8.1 --height 4 --sizeupper 15
			
		
	# (1d) also characterize the clusters by running GO-analysis on all DE genes (both up + down) (Fig. 1d, Ext. data Figs 4 and 5)
	# ----------------
	mkdir -p "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO"
	
	awk -F$'\t' '{OFS=FS} {print $2,$1}' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up.txt > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up_flipcol.txt
	awk -F$'\t' '{OFS=FS} $2 ~ /endo/ {print $0}' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up_flipcol.txt > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up_flipcol_endo_only.txt
	awk -F$'\t' '{OFS=FS} $2 ~ /seedcoat/ {print $0}' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up_flipcol.txt > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up_flipcol_seedcoat_only.txt
	
	tail -n+2 $outdir/all_detected_genes.txt > $outdir/all_detected_genes_noheader.txt
	background="$outdir/all_detected_genes_noheader.txt"
	
	echo "Analyzing GO-terms in all_endo_up, p < 0.005"
	mkdir -p $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_endo_up_p005
	genelist="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up_flipcol_endo_only.txt"
	$path_to_scripts/run_topGO.R "$genelist" $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_endo_up_p005/all_endo_up_p005 --background "$background" --expr_data $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_zscores_fxheader.txt --topN 5 --pval 0.005 --pvalmax 6 --fillupper 4 --filllower -4 --useallGOtermgenes > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_endo_up_p005/log.txt
	
	echo "Analyzing GO-terms in all_seedcoat_up, p < 0.005"
	mkdir -p $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_seedcoat_up_p005
	genelist="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_up_flipcol_seedcoat_only.txt"
	$path_to_scripts/run_topGO.R "$genelist" $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_seedcoat_up_p005/all_seedcoat_up_p005 --background "$background" --expr_data $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_zscores_fxheader.txt --topN 5 --pval 0.005 --pvalmax 6 --fillupper 4 --filllower -4 --useallGOtermgenes > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_seedcoat_up_p005/log.txt
		
	# repeat for the downreg. genes
	awk -F$'\t' '{OFS=FS} {print $2,$1}' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn.txt > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn_flipcol.txt
	awk -F$'\t' '{OFS=FS} $2 ~ /endo/ {print $0}' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn_flipcol.txt > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn_flipcol_endo_only.txt
	awk -F$'\t' '{OFS=FS} $2 ~ /seedcoat/ {print $0}' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn_flipcol.txt > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn_flipcol_seedcoat_only.txt
	
	echo "Analyzing GO-terms in all_endo_dwn, p < 0.005"
	mkdir -p $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_endo_dwn_p005
	genelist="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn_flipcol_endo_only.txt"
	$path_to_scripts/run_topGO.R "$genelist" $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_endo_dwn_p005/all_endo_dwn_p005 --background "$background" --expr_data $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_zscores_fxheader.txt --pval 0.005 --topN 5 --pvalmax 6 --fillupper 4 --filllower -4 --useallGOtermgenes > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_endo_dwn_p005/log.txt
	
	echo "Analyzing GO-terms in all_seedcoat_dwn, p < 0.005"
	mkdir -p $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_seedcoat_dwn_p005
	genelist="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/all_sig_gene_lists/all_clusters_dwn_flipcol_seedcoat_only.txt"
	$path_to_scripts/run_topGO.R "$genelist" $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_seedcoat_dwn_p005/all_seedcoat_dwn_p005 --background "$background" --expr_data $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_zscores_fxheader.txt --pval 0.005 --topN 5 --pvalmax 6 --fillupper 4 --filllower -4 --useallGOtermgenes > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_seedcoat_dwn_p005/log.txt
		
	# repeat for the 60 clusters
	awk '$2 != "NA"' "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_heatmap_k60_zscores_roworder.txt" > "$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_k60_clusters_final.txt"
	genelist="$outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_k60_clusters_final.txt"

	echo "Analyzing GO-terms in the 60 gene clusters, p < 0.005"
	mkdir -p $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_cluster_DEGs_p005
	$path_to_scripts/run_topGO.R "$genelist" $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_cluster_DEGs_p005/all_cluster_DEGs_p005 --pval 0.005 --pvalmax 6 --background "$background" --topN 5 --sampleorder 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60 > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/topGO/all_cluster_DEGs_p005/log.txt
				
	
	# (2) how does allelic expression (mat/pat) vary in our data across the nuclei clusters, cell cycle, and factors?
	# ----------------
	sed 's/^\t/locus_name\t/1' $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_mcounts_norm.txt > $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_mcounts_norm_fxheader.txt
	sed 's/^\t/locus_name\t/1' $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_mcounts_norm.txt > $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_mcounts_norm_fxheader.txt
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_mcounts_norm_fxheader.txt $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_mcounts_norm_fxheader.txt locus_name $outdir/imprinting_analysis/all_mcounts_merged.txt

	sed 's/^\t/locus_name\t/1' $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_pcounts_norm.txt > $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_pcounts_norm_fxheader.txt
	sed 's/^\t/locus_name\t/1' $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_pcounts_norm.txt > $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_pcounts_norm_fxheader.txt
	$path_to_scripts/merge_by_column.R $outdir/imprinting_analysis/ASE_CxV/ASE_CxV_pcounts_norm_fxheader.txt $outdir/imprinting_analysis/ASE_VxC/ASE_VxC_pcounts_norm_fxheader.txt locus_name $outdir/imprinting_analysis/all_pcounts_merged.txt

	mkdir -p "$outdir/imprinting_variability"
		
	# (2a) repeat permutation analysis in 'imprinting mode' - 'separate', over the same sets of genes and factors
	# for each analysis, do normal run + run where genes/rows are automatically ordered according to matched run from above (matchH matches hierarchical clustering above, matchK matches k-means)
	pids=(); logfiles=()
	for ((i=0;i<${#listname[@]};++i)); do
		echo "Analyzing allelic expression pattern of genes in list: ${listname[i]}s"
		touse=$( echo "${listname[i]}" | sed 's/ /_/g' )
	
		# endo nuclei over clusters, controlling for cross and wash
		factype="endo_nuc_control_cross_wash"; clustype="expr_var_over_clusters"
		cluslist="$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_clusters.txt"
		factors="$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross_wash.txt"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK"
		kroworder=$( find $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse -regex '.*_kmeans_heatmap_k[0-9]+_zscores_roworder.txt' -regextype posix-extended )
		kcolorder=$( find $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse -regex '.*_kmeans_heatmap_k[0-9]+_zscores_colorder.txt' -regextype posix-extended )
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --kmeans --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_colorder.txt --genefile $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_roworder.txt --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd4="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
				
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_log.txt" -K "$cmd4" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_log.txt" )
		if [ "${listname[i]}" = "CC DEG" ] || [ "${listname[i]}" = "cluster DEG" ]; then
			kroworder="$outdir/expression_variability/${clustype}/all_nuc_control_cross_tissue/$touse/${touse}_kmeans60_heatmap_k60_zscores_roworder.txt"
			kcolorder="$outdir/expression_variability/${clustype}/all_nuc_control_cross_tissue/$touse/${touse}_kmeans60_heatmap_k60_zscores_colorder.txt"
			
			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60"
			awk '$0 !~ /seedcoat/' $kcolorder > $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60/${touse}_colorder_touse.txt
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60/${touse}_colorder_touse.txt --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_log.txt" )

			# also match row order of endo control cross
			kroworder="$outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_kmeans60_heatmap_k60_zscores_roworder.txt"
			kcolorder="$outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_kmeans60_heatmap_k60_zscores_colorder.txt"
			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo"
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo/${touse} --method separate --customyorder $kcolorder --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_log.txt" )
		fi
	
		# expression patterns in endo nuclei, by all factors
		factype="endo_nuc_control_cross_wash_ploidy"; clustype="expr_var_over_factors"
		factors="$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross_wash_ploidy.txt"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier"
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier/$touse --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --allowmissinggenes --plotupperZ 6 --clustersep 4 --kmax 20 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" )

		# endo nuclei over cell cycle, controlling for cross only
		factype="endo_nuc_control_cross"; clustype="expr_var_over_CC"
		cluslist="$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt"
		factors="$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross.txt"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_max3"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_max3"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_max3"
		kroworder=$( find $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse -regex '.*_kmeans_heatmap_k[0-9]+_zscores_roworder.txt' -regextype posix-extended )
		kcolorder=$( find $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse -regex '.*_kmeans_heatmap_k[0-9]+_zscores_colorder.txt' -regextype posix-extended )
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --kmeans --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_colorder.txt --genefile $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_roworder.txt --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd4="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_log.txt" -K "$cmd4" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_log.txt" )

		# repeat, setting Zmax to 3 and Dmax to 3
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_max3/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --kmeans --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_max3/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_colorder.txt --genefile $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_roworder.txt --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd4="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_max3/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_max3_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_max3_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_max3_log.txt" -K "$cmd4" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_log.txt" )

		if [ "${listname[i]}" = "CC DEG" ] || [ "${listname[i]}" = "cluster DEG" ]; then		
			kroworder="$outdir/expression_variability/${clustype}/all_nuc_control_cross_tissue/$touse/${touse}_kmeans60_heatmap_k60_zscores_roworder.txt"
			kcolorder="$outdir/expression_variability/${clustype}/all_nuc_control_cross_tissue/$touse/${touse}_kmeans60_heatmap_k60_zscores_colorder.txt"
			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60"
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_log.txt" )

			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_max3"
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_max3/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_max3_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_max3_log.txt" )

			# also match row order of endo control cross
			kroworder="$outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_kmeans60_heatmap_k60_zscores_roworder.txt"
			kcolorder="$outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_kmeans60_heatmap_k60_zscores_colorder.txt"
			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo"
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_log.txt" )

			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_max3"
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_max3/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_max3_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_max3_log.txt" )
		fi
	
		# endo nuclei over cell cycle, controlling for cross and wash
		factype="endo_nuc_control_cross_wash"; clustype="expr_var_over_CC"
		cluslist="$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt"
		factors="$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross_wash.txt"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_max3"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_max3"
		mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_max3"
		kroworder=$( find $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse -regex '.*_kmeans_heatmap_k[0-9]+_zscores_roworder.txt' -regextype posix-extended )
		kcolorder=$( find $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse -regex '.*_kmeans_heatmap_k[0-9]+_zscores_colorder.txt' -regextype posix-extended )
		
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --kmeans --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_colorder.txt --genefile $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_roworder.txt --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd4="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_log.txt" -K "$cmd4" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_log.txt" )

		# repeat, setting Zmax to 3 and Dmax to 3
		cmd1="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd2="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_max3/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --kmeans --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd3="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_max3/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_colorder.txt --genefile $outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_hier_heatmap_zscores_roworder.txt --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		cmd4="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_max3/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3_log.txt" -K "$cmd1" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_max3_log.txt" -K "$cmd2" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_max3_log.txt" -K "$cmd3" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_log.txt" )
		bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_max3_log.txt" -K "$cmd4" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK_log.txt" )

		if [ "${listname[i]}" = "CC DEG" ] || [ "${listname[i]}" = "cluster DEG" ]; then		
			kroworder="$outdir/expression_variability/${clustype}/all_nuc_control_cross_tissue/$touse/${touse}_kmeans60_heatmap_k60_zscores_roworder.txt"
			kcolorder="$outdir/expression_variability/${clustype}/all_nuc_control_cross_tissue/$touse/${touse}_kmeans60_heatmap_k60_zscores_colorder.txt"
			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60"
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_log.txt" )

			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_max3"
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_max3/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_max3_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_max3_log.txt" )

			# also match row order of endo control cross
			kroworder="$outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_kmeans60_heatmap_k60_zscores_roworder.txt"
			kcolorder="$outdir/expression_variability/${clustype}/endo_nuc_control_cross/$touse/${touse}_kmeans60_heatmap_k60_zscores_colorder.txt"
			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo"
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_log.txt" )

			mkdir "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_max3"
			cmd5="$path_to_scripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_max3/${touse} --method separate --customxorder --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --customyorder $kcolorder --genefile $kroworder --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
			bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_max3_log.txt" -K "$cmd5" & pids+=( $! ); logfiles+=( "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchK60_endo_max3_log.txt" )
		fi
	
	done

	# wait for all jobs to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "cluster_gene_expression for assessing imprinting variability failed, see logfile ${logfiles[i]}" "$log"
	done


	# reviewer mentioned looking at cell cycle effects within chalazal, do that here using only chalazal nuclei - CxV E12,13,14 and VxC E1,6
	# control for cross only or for cross & wash (Fig. S16e)
	# ----------------
	awk -F$'\t' '{OFS=FS} NR == 1 {print "locus_name","cluster"} $2 == "CxV_endo_12" || $2 == "CxV_endo_13" || $2 == "CxV_endo_14" || $2 == "VxC_endo_1" || $2 == "VxC_endo_6" {print $0}' "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_clusters.txt" > "$outdir/expression_variability/nuclei_lists/chalazal_4DAP_nuclei_clusters_wheader.txt"
	awk -F$'\t' '{OFS=FS} NR == 1 {print "locus_name","cell_cycle"} {print $0}' "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" > "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle_wheader.txt"
	awk -F$'\t' '{OFS=FS} NR == 1 {print "locus_name","cross","wash"} {print $0}' "$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross_wash.txt" > "$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross_wash_wheader.txt"

	$path_to_scripts/merge_by_column.R "$outdir/expression_variability/nuclei_lists/chalazal_4DAP_nuclei_clusters_wheader.txt" "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle_wheader.txt" locus_name "$outdir/expression_variability/nuclei_lists/tmp.txt" --tokeep merged
	$path_to_scripts/merge_by_column.R "$outdir/expression_variability/nuclei_lists/tmp.txt" "$outdir/expression_variability/factorlists/endo_4DAP_nuclei_cross_wash_wheader.txt" locus_name "tmp.txt" --tokeep merged

	cut -f1,3 tmp.txt | tail -n+2 > "$outdir/expression_variability/nuclei_lists/chalazal_4DAP_nuclei_cell_cycle.txt"
	cut -f1,4,5 tmp.txt | tail -n+2 > "$outdir/expression_variability/factorlists/chalazal_4DAP_nuclei_cross_wash.txt"
	cut -f1,4 tmp.txt | tail -n+2 > "$outdir/expression_variability/factorlists/chalazal_4DAP_nuclei_cross.txt"

	awk -F$'\t' '{OFS=FS} NR==1 {print "locus_name"} $2 == "chalazal_PEG" {print $1}' "$outdir/expression_variability/expr_var_over_clusters/endo_nuc_control_cross/all_PEG/all_PEG_hier_chalazal_vs_nonchalazal.txt" > "$outdir/expression_variability/genelists/chalazal_PEGs_wheader.txt"
	awk -F$'\t' '{OFS=FS} NR==1 {print "locus_name"} $2 == "nonchalazal_PEG" {print $1}' "$outdir/expression_variability/expr_var_over_clusters/endo_nuc_control_cross/all_PEG/all_PEG_hier_chalazal_vs_nonchalazal.txt" > "$outdir/expression_variability/genelists/nonchalazal_PEGs_wheader.txt"
	awk -F$'\t' '{OFS=FS} NR==1 {print "locus_name","ptype"} {print $0}' "$outdir/expression_variability/genelists/all_PEG_list.txt" > "$outdir/expression_variability/genelists/all_PEG_list_wheader.txt"
	$path_to_scripts/merge_by_column.R "$outdir/expression_variability/genelists/chalazal_PEGs_wheader.txt" "$outdir/expression_variability/genelists/all_PEG_list_wheader.txt" locus_name - --tokeep merged | tail -n+2 > "$outdir/expression_variability/genelists/chalazal_PEGs.txt"
	$path_to_scripts/merge_by_column.R "$outdir/expression_variability/genelists/nonchalazal_PEGs_wheader.txt" "$outdir/expression_variability/genelists/all_PEG_list_wheader.txt" locus_name - --tokeep merged | tail -n+2 > "$outdir/expression_variability/genelists/nonchalazal_PEGs.txt"

	listname=(); listlist=()
	listname+=( "endo gr SC" ); listlist+=( "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/intersect_genes_endo_gr_seedcoat.txt" )
	listname+=( "endo le SC" ); listlist+=( "$outdir/DE_analysis/other/seedcoat_vs_endo_all_CV_4DAP/intersect_genes_endo_le_seedcoat.txt" )
	listname+=( "endo MEG" ); listlist+=( "$outdir/expression_variability/genelists/endo_MEGs.txt" )
	listname+=( "endo MEG v2" ); listlist+=( "$outdir/expression_variability/genelists/endo_MEGs_v2.txt" )
	listname+=( "chromatin related" ); listlist+=( "$outdir/chromdb_loci_touse.txt" )
	listname+=( "no bias subset" ); listlist+=( "$outdir/expression_variability/genelists/no_bias_subset.txt" )
	listname+=( "cluster DEG" ); listlist+=( "$outdir/expression_variability/genelists/genes_pass_cluster_nobatch.txt" )
	listname+=( "CC DEG" ); listlist+=( "$outdir/expression_variability/genelists/CC_sig_genes.txt" )
	listname+=( "all PEG" ); listlist+=( "$outdir/expression_variability/genelists/all_PEG_list.txt" )
	listname+=( "all MEG" ); listlist+=( "$outdir/expression_variability/genelists/all_MEG_list.txt" )
	listname+=( "chalazal PEG" ); listlist+=( "$outdir/expression_variability/genelists/chalazal_PEGs.txt" )
	listname+=( "nonchalazal PEG" ); listlist+=( "$outdir/expression_variability/genelists/nonchalazal_PEGs.txt" )

	for ((i=0;i<${#listname[@]};++i)); do
	echo "Analyzing allelic expression pattern in cell cycle (chalazal nuclei only) of genes in list: ${listname[i]}s"
	touse=$( echo "${listname[i]}" | sed 's/ /_/g' )

	factype="chalazal_nuc_control_cross"; clustype="expr_var_over_CC"
	cluslist="$outdir/expression_variability/nuclei_lists/chalazal_4DAP_nuclei_cell_cycle.txt"
	factors="$outdir/expression_variability/factorlists/chalazal_4DAP_nuclei_cross.txt"
	mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3"
	mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_max3"
	mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH"; mkdir -p "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_matchH_max3"

	cmd1="$path_to_wipscripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
	cmd2="$path_to_wipscripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --plotupperD 5 --allowmissinggenes --hmethod average --kmeans --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"

	bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_log.txt" "$cmd1"
	bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_log.txt" "$cmd2"

	# repeat, setting Zmax to 3 and Dmax to 3
	cmd1="$path_to_wipscripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"
	cmd2="$path_to_wipscripts/cluster_gene_expression.R $outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_max3/${touse} --method separate --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --clustersamples --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --plotupperD 3 --allowmissinggenes --hmethod average --kmeans --clustersep 4 --kmax 60 --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C --colorsD \#0097E1,\#66CDF6,\#C8FAFD,\#EFFDFF,white,\#EDE6F8,\#EAB9F6,\#EC65E8,\#DC2FF0"

	bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_hier_max3_log.txt" "$cmd1"
	bsub -o "$outdir/imprinting_variability/${clustype}/${factype}/$touse/${touse}_kmeans_max3_log.txt" "$cmd2"
	done

	# also do total expression analysis
	# expression patterns in endo nuclei, over cell cycle, controlling for cross only
	for ((i=0;i<${#listname[@]};++i)); do
	echo "Analyzing allelic expression pattern in cell cycle (chalazal nuclei only) of genes in list: ${listname[i]}s"
	touse=$( echo "${listname[i]}" | sed 's/ /_/g' )

	factype="chalazal_nuc_control_cross"; clustype="expr_var_over_CC"
	cluslist="$outdir/expression_variability/nuclei_lists/chalazal_4DAP_nuclei_cell_cycle.txt"
	factors="$outdir/expression_variability/factorlists/chalazal_4DAP_nuclei_cross.txt"
	mkdir -p "$outdir/expression_variability/$clustype/$factype/$touse"
	mkdir -p "$outdir/expression_variability/$clustype/$factype/${touse}_max3"

	cmd1="$path_to_wipscripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
	cmd2="$path_to_wipscripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/$touse/${touse}_kmeans --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 6 --kmeans --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
	bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_hier_log.txt" "$cmd1"
	bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_kmeans_log.txt" "$cmd2"

	# repeat with max Z = 3
	cmd1="$path_to_wipscripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/${touse}_max3/${touse}_hier --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
	cmd2="$path_to_wipscripts/cluster_gene_expression.R $outdir/expression_variability/${clustype}/${factype}/${touse}_max3/${touse}_kmeans --expr_matrix $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile $cluslist --genefile ${listlist[i]} --factors $factors --showrownames --NAcolor grey50 --seed 123456 --nreps 1000 --NAreplacedist 0 --plotupperZ 3 --kmeans --clustersep 4 --kmax 60 --colorder G0,G1,G1toS,S,G2,M --colorsZ \#3D42B6,\#4F74B0,\#80ABCD,\#B4D8E7,\#E4F2F7,white,\#FFFEC6,\#F9E19B,\#F2B26E,\#E3754F,\#D1382C"
	bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_hier_log.txt" "$cmd1"
	bsub -o "$outdir/expression_variability/${clustype}/${factype}/${touse}_max3_kmeans_log.txt" "$cmd2"
	done

	# some examples...
	mkdir -p "$outdir/expression_variability/example_plots"
	gene="AT4G16810"
	$path_to_wipscripts/single_cell_RNAseq_plots.R lin $outdir/expression_variability/example_plots/line_${gene}_CC_CZE_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/chalazal_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M

	gene="AT4G11400"
	$path_to_wipscripts/single_cell_RNAseq_plots.R lin $outdir/expression_variability/example_plots/line_${gene}_CC_CZE_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/chalazal_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M
	
	gene="AT5G28850"
	$path_to_wipscripts/single_cell_RNAseq_plots.R lin $outdir/expression_variability/example_plots/line_${gene}_CC_CZE_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/chalazal_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M	
	
fi


# ----------------------
# Step 9: any additional analyses/plots
# ----------------------
if "$skip9"; then printf "\nSkipping step 9 by user request...\n" | tee -a "$log"
else printf "\nStep 9: any additional analyses/plots\n" | tee -a "$log"; ts=$(date +%s)

	# GO term analyses over MEGs and PEGs
	# ----------------
	mkdir -p $outdir/imprinting_analysis/topGO
	
	# background is set of all genes that could be eval as imprinted
	awk -F$'\t' '{OFS=FS} NR != 1 && $8 != "no data" {print $1}' $outdir/imprinting_analysis/final_status.txt > $outdir/imprinting_analysis/topGO/background.txt
	background="$outdir/imprinting_analysis/topGO/background.txt"
	awk -F$'\t' '{OFS=FS} NR != 1 && ($2 == "no bias" || $2 ~ /MEG/ || $2 ~ /PEG/) {print $0}' $outdir/imprinting_analysis/final_status_only.txt > $outdir/imprinting_analysis/topGO/impr_status_list.txt
	awk -F$'\t' '{OFS=FS} NR != 1 && ($2 ~ /MEG/ || $2 ~ /PEG/) {print $0}' $outdir/imprinting_analysis/final_status_only.txt > $outdir/imprinting_analysis/topGO/impr_status_list_MEGs_PEGs_only.txt

	$path_to_scripts/run_topGO.R "$outdir/imprinting_analysis/topGO/impr_status_list.txt" $outdir/imprinting_analysis/topGO/topGO_impr --background "$background" --topN 5 --pval 0.0001 --pvalmax 6 #> $outdir/imprinting_analysis/topGO/topGO_impr_log.txt
	
	# add endo MEGs as own category...
	awk -F$'\t' '{OFS=FS} {print $1,"endo MEG"}' $outdir/expression_variability/genelists/endo_MEGs.txt > $outdir/imprinting_analysis/topGO/impr_status_list_w_endoMEG.txt
	cat $outdir/imprinting_analysis/topGO/impr_status_list.txt >> $outdir/imprinting_analysis/topGO/impr_status_list_w_endoMEG.txt
	sort -u -k1,1 $outdir/imprinting_analysis/topGO/impr_status_list_w_endoMEG.txt | uniq | awk -F$'\t' '$2 == "endo MEG"' > $outdir/imprinting_analysis/topGO/impr_status_list_w_endoMEG_touse.txt
	rm $outdir/imprinting_analysis/topGO/impr_status_list_w_endoMEG.txt

	$path_to_scripts/run_topGO.R "$outdir/imprinting_analysis/topGO/impr_status_list_w_endoMEG.txt" $outdir/imprinting_analysis/topGO/topGO_impr_w_endoMEG --background "$background" --topN 5 --pval 0.0001 --pvalmax 6 #> $outdir/imprinting_analysis/topGO/topGO_impr_w_endoMEG_log.txt

	# do all PEGs, chalazal PEGs and non-chalazal PEGs instead (note - made file by hand based off clustering in fig. 3)
	$path_to_scripts/run_topGO.R "$outdir/expression_variability/expr_var_over_clusters/endo_nuc_control_cross/all_PEG/all_PEG_hier_chalazal_vs_nonchalazal.txt" $outdir/imprinting_analysis/topGO/topGO_PEGs_wchalazal --background "$background" --topN 5 --pval 0.001

	# show MEG, PEG, non-chalazal PEG, and chalazal PEG GO term enrichment
	awk -F$'\t' '$2 ~ /MEG/' $outdir/imprinting_analysis/topGO/impr_status_list_w_endoMEG.txt | sed 's/strong MEG/MEG/g' | sed 's/weak MEG/MEG/g' >| $outdir/imprinting_analysis/topGO/MEGs_w_endoMEG.txt
	cat $outdir/expression_variability/expr_var_over_clusters/endo_nuc_control_cross/all_PEG/all_PEG_hier_chalazal_vs_nonchalazal.txt $outdir/imprinting_analysis/topGO/MEGs_w_endoMEG.txt >| $outdir/imprinting_analysis/topGO/MEGs_PEGs_byregion.txt
	$path_to_wipscripts/run_topGO.R "$outdir/imprinting_analysis/topGO/MEGs_PEGs_byregion.txt" $outdir/imprinting_analysis/topGO/MEGs_PEGs_byregion --background "$background" --topN 5 --pval 0.001 --pvalmax 6


	# Make "saturation analysis" type plots, for any reads (0) or at least 5 reads (5) - Fig. S1
	# ----------------
	mkdir -p "$outdir/gene_detection"
	rm -rf "$outdir/gene_detection/per_lib_summary.txt"
	
	for mm in "0" "5" "10"; do
		for allele in "all" "Col" "Ler" "Cvi"; do
			echo ""
			echo "Saturation analysis for genes with at least $mm coverage, for ${allele} counts"
			[ "$mm" = "5" ] && ss="min5" || ss="min0"
			[ "$mm" = "10" ] && ss="min10"
			echo "sample	tot_genes_detected" > "$outdir/gene_detection/saturation_analysis_${allele}_${ss}.txt"
			echo "lib	${allele}_min${mm}" > "$outdir/gene_detection/tmp_summary.txt"
			for ((i=0;i<${#stubname[@]};++i)); do
				echo "Processing ${stubname[i]}"
				loc="$outdir/${project[i]}_map/${stubname[i]}_map"
				if [ -f "$loc/htseq_count/${stubname[i]}_counts_dedup_${allele}_genes.txt" ]; then
					loc="$outdir/${project[i]}_map/${stubname[i]}_map"
					echo "locus_name" > "$outdir/gene_detection/tmp.txt"					
					awk -F$'\t' -v a="${mm}" '{OFS=FS} $2>a && $1 !~ /__/ && $1 !~ /ATM/ && $1 !~ /ATC/ {print $1}' "$loc/htseq_count/${stubname[i]}_counts_dedup_${allele}_genes.txt" >> "$outdir/gene_detection/tmp.txt"

					nn=$( wc -l "$outdir/gene_detection/tmp.txt" | awk '{print $1}' )
					echo "${stubname[i]}	$nn" >> "$outdir/gene_detection/tmp_summary.txt"
#					echo "${stubname[i]} has $nn detectable genes"

					if [ ! -z "$outdir/gene_detection/tmp.txt" ]; then
						[ -f "$outdir/gene_detection/all_${allele}_${ss}_counts.txt" ] && { $path_to_scripts/merge_by_column.R "$outdir/gene_detection/all_${allele}_${ss}_counts.txt" "$outdir/gene_detection/tmp.txt" locus_name "$outdir/gene_detection/all_${allele}_${ss}_counts.txt" > /dev/null; }
						[ ! -f "$outdir/gene_detection/all_${allele}_${ss}_counts.txt" ] && { cat "$outdir/gene_detection/tmp.txt" > "$outdir/gene_detection/all_${allele}_${ss}_counts.txt"; }
					fi
					if [ -f "$outdir/gene_detection/all_${allele}_${ss}_counts.txt" ]; then
						numdet=$( wc -l "$outdir/gene_detection/all_${allele}_${ss}_counts.txt" | awk '{print $1}' )
#						echo "A total of $numdet genes have been detected so far"
						echo "${stubname[i]}	$numdet" >> "$outdir/gene_detection/saturation_analysis_${allele}_${ss}.txt"
					fi
				fi
			done
			[ ! -f "$outdir/gene_detection/per_lib_summary.txt" ] && mv "$outdir/gene_detection/tmp_summary.txt" "$outdir/gene_detection/per_lib_summary.txt" || { $path_to_scripts/merge_by_column.R "$outdir/gene_detection/per_lib_summary.txt" "$outdir/gene_detection/tmp_summary.txt" lib "$outdir/gene_detection/per_lib_summary.txt" > /dev/null; }
			
		done
	done
	
	
	# make metaplots of RNA-seq coverage over gene bodies to show 3' bias, contrast this with data from
	# the Col-Cvi RILs paper (use all Col and Cvi replicates) - Fig. S2
	# ----------------
	mkdir -p "$outdir/bamCoverage/bamCoverage_dedup_all"
	mkdir -p "$outdir/bamCoverage/bamCoverage_dedup_mat"
	mkdir -p "$outdir/bamCoverage/bamCoverage_dedup_pat"
	rm -rf "$outdir/bamCoverage/lsf_logs"
	mkdir -p "$outdir/bamCoverage/lsf_logs"
	rm -rf "$outdir/bamCoverage/scripts"
	mkdir -p "$outdir/bamCoverage/scripts"
	
	pids=(); logfiles=(); i=0; scriptno=1
	echo "#!/bin/bash" > "$outdir/bamCoverage/scripts/all_script${scriptno}.sh"
	echo "" >> "$outdir/bamCoverage/scripts/all_script${scriptno}.sh"
	
	for ((i=0;i<${#stubname[@]};++i)); do
		ff="$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments_dedup.bam"
		j=$(( $j + 1 )); k=$(( $i + 1 ))
		cmd="bamCoverage -b $ff -o $outdir/bamCoverage/bamCoverage_dedup_all/${stubname[i]}_cov.bw -bs 1 --normalizeUsing CPM"
		if [ "$j" = "20" ] || [ "$k" = "${#stubname[@]}" ]; then
			# submit last set of jobs
			echo "$cmd" >> "$outdir/bamCoverage/scripts/all_script${scriptno}.sh"
			echo "" >> "$outdir/bamCoverage/scripts/all_script${scriptno}.sh"
			chmod 755 "$outdir/bamCoverage/scripts/all_script${scriptno}.sh"
			bsub -o "$outdir/bamCoverage/lsf_logs/all_script${scriptno}_log.txt" -K "$outdir/bamCoverage/scripts/all_script${scriptno}.sh" & pids+=( $! ); logfiles+=( "$outdir/bamCoverage/lsf_logs/all_script${scriptno}_log.txt" )
		
			# start new script (if we haven't reached end of array)
			if [ "$k" != "${#stubname[@]}" ]; then
				scriptno=$(( $scriptno + 1 ))
				j=0
				echo "#!/bin/bash" > "$outdir/bamCoverage/scripts/all_script${scriptno}.sh"
				echo "" >> "$outdir/bamCoverage/scripts/all_script${scriptno}.sh"
			fi
		else
			echo "$cmd" >> "$outdir/bamCoverage/scripts/all_script${scriptno}.sh"
		fi											
	done
		
	# wait for all jobs to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "bamCoverage failed for all counts, see log file ${logfiles[i]}" "$log"
	done
	
					
	# Repeat for maternal counts
	pids=(); logfiles=(); i=0; scriptno=1
	echo "#!/bin/bash" > "$outdir/bamCoverage/scripts/mat_script${scriptno}.sh"
	echo "" >> "$outdir/bamCoverage/scripts/mat_script${scriptno}.sh"
	
	for ((i=0;i<${#stubname[@]};++i)); do
		if [[ "${cross[i]}" = "Cx"* ]]; then
			ff="$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments_to_Col_dedup.bam"
		elif [[ "${cross[i]}" = "Lx"* ]]; then
			ff="$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments_to_Ler_dedup.bam"
		else
			ff="$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments_to_Cvi_dedup.bam"
		fi
		j=$(( $j + 1 )); k=$(( $i + 1 ))
		cmd="bamCoverage -b $ff -o $outdir/bamCoverage/bamCoverage_dedup_mat/${stubname[i]}_mat_cov.bw -bs 1 --normalizeUsing CPM"
		if [ "$j" = "20" ] || [ "$k" = "${#stubname[@]}" ]; then
			# submit last set of jobs
			echo "$cmd" >> "$outdir/bamCoverage/scripts/mat_script${scriptno}.sh"
			echo "" >> "$outdir/bamCoverage/scripts/mat_script${scriptno}.sh"
			chmod 755 "$outdir/bamCoverage/scripts/mat_script${scriptno}.sh"
			bsub -o "$outdir/bamCoverage/lsf_logs/mat_script${scriptno}_log.txt" -K "$outdir/bamCoverage/scripts/mat_script${scriptno}.sh" & pids+=( $! ); logfiles+=( "$outdir/bamCoverage/lsf_logs/mat_script${scriptno}_log.txt" )
		
			# start new script (if we haven't reached end of array)
			if [ "$k" != "${#stubname[@]}" ]; then
				scriptno=$(( $scriptno + 1 ))
				j=0
				echo "#!/bin/bash" > "$outdir/bamCoverage/scripts/mat_script${scriptno}.sh"
				echo "" >> "$outdir/bamCoverage/scripts/mat_script${scriptno}.sh"
			fi
		else
			echo "$cmd" >> "$outdir/bamCoverage/scripts/mat_script${scriptno}.sh"
		fi											
	done
		
	# wait for all jobs to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "bamCoverage failed for mat counts, see log file ${logfiles[i]}" "$log"
	done
	
	
	# Repeat for paternal counts
	pids=(); logfiles=(); i=0; scriptno=1
	echo "#!/bin/bash" > "$outdir/bamCoverage/scripts/pat_script${scriptno}.sh"
	echo "" >> "$outdir/bamCoverage/scripts/pat_script${scriptno}.sh"
	
	for ((i=0;i<${#stubname[@]};++i)); do
		if [[ "${cross[i]}" = *"xC" ]]; then
			ff="$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments_to_Col_dedup.bam"
		elif [[ "${cross[i]}" = *"xL" ]]; then
			ff="$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments_to_Ler_dedup.bam"
		else
			ff="$outdir/${project[i]}_map/${stubname[i]}_map/STAR/${stubname[i]}_unique_alignments_to_Cvi_dedup.bam"
		fi
		j=$(( $j + 1 )); k=$(( $i + 1 ))
		cmd="bamCoverage -b $ff -o $outdir/bamCoverage/bamCoverage_dedup_pat/${stubname[i]}_pat_cov.bw -bs 1 --normalizeUsing CPM"
		if [ "$j" = "20" ] || [ "$k" = "${#stubname[@]}" ]; then
			# submit last set of jobs
			echo "$cmd" >> "$outdir/bamCoverage/scripts/pat_script${scriptno}.sh"
			echo "" >> "$outdir/bamCoverage/scripts/pat_script${scriptno}.sh"
			chmod 755 "$outdir/bamCoverage/scripts/pat_script${scriptno}.sh"
			bsub -o "$outdir/bamCoverage/lsf_logs/pat_script${scriptno}_log.txt" -K "$outdir/bamCoverage/scripts/pat_script${scriptno}.sh" & pids+=( $! ); logfiles+=( "$outdir/bamCoverage/lsf_logs/pat_script${scriptno}_log.txt" )
		
			# start new script (if we haven't reached end of array)
			if [ "$k" != "${#stubname[@]}" ]; then
				scriptno=$(( $scriptno + 1 ))
				j=0
				echo "#!/bin/bash" > "$outdir/bamCoverage/scripts/pat_script${scriptno}.sh"
				echo "" >> "$outdir/bamCoverage/scripts/pat_script${scriptno}.sh"
			fi
		else
			echo "$cmd" >> "$outdir/bamCoverage/scripts/pat_script${scriptno}.sh"
		fi											
	done
		
	# wait for all jobs to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "bamCoverage failed for pat counts, see log file ${logfiles[i]}" "$log"
	done
	
	# repeat for the Col, Cvi replicates from RIL paper
	pids=(); logfiles=()
	for samp in col_1 col_2 col_3 cvi_1 cvi_2 cvi_3; do
		ff="/archive/gehring/2017.01.24-15807/solexa_gehring/colette/RILs/rnaseq_analysis_new/${samp}_map/${samp}_all_unique_alignments_dedup.bam"
		cmd="bamCoverage -b $ff -o $outdir/bamCoverage/${samp}_RILs_cov.bw -bs 1 --normalizeUsing CPM"
		bsub -o "$outdir/bamCoverage/lsf_logs/${samp}_FILs_log.txt" -K "$cmd" & pids+=( $! ); logfiles+=( "$outdir/bamCoverage/lsf_logs/${samp}_FILs_log.txt" )
	done
	
	# wait for all jobs to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "bamCoverage failed for pat counts, see log file ${logfiles[i]}" "$log"
	done
	
	# now use deeptools to make coverage plots
	# get BED file first
	cat $gtf_genes | tr ';' '\t' | awk -F$'\t' '{OFS=FS} {if ($3 == "exon" && $9 ~ /gene_id/) {print $1,$4,$5,$7,$9} else if ($3 == "exon" && $10 ~ /gene_id/) {print $1,$4,$5,$7,$10}}' | sed 's/gene_id "//g' | sed 's/"$//g' | awk -F$'\t' '{OFS=FS} {if ($5 in s) {if ($2 < s[$5]) {s[$5]=$2}; if ($3 > e[$5]) {e[$5]=$3}} else {s[$5]=$2; e[$5]=$3; c[$5]=$1; ss[$5]=$4}} END {for (a in e) {print c[a],s[a],e[a],a,"0",ss[a]}}' | sort -k1,1 -k2n,2 | sed 's/ //g' > $outdir/TAIR10_plus_araport11_nonoverlapping.bed
	
	# drop all ChrC/M
	awk '$1 != "ChrM" && $1 != "ChrC"' $outdir/TAIR10_plus_araport11_nonoverlapping.bed > $outdir/TAIR10_plus_araport11_nonoverlapping_noChrCM.bed
	
	# use introns from TAIR10 (where nearly all PCGs are from anyway); add back strand info
	printf "chr\tstart\tend\tlocus_name\n" >| $outdir/tmp.txt
	awk '$1 != "ChrM" && $1 != "ChrC"' /lab/solexa_gehring/genomes/A_thaliana/annotations/TAIR10/introns.bed | sed 's/ //g' | cut -f1 -d'.' >> $outdir/tmp.txt
	
	printf "locus_name\tstrand\n" >| $outdir/tmp2.txt
	cut -f4,6 /lab/solexa_gehring/genomes/A_thaliana/annotations/TAIR10/genes.bed >> $outdir/tmp2.txt
	cut -f4,6 /lab/solexa_gehring/genomes/A_thaliana/annotations/TAIR10/pseudogenes.bed >> $outdir/tmp2.txt
	cut -f4,6 /lab/solexa_gehring/genomes/A_thaliana/annotations/TAIR10/TEgenes.bed >> $outdir/tmp2.txt

	$path_to_scripts/merge_by_column.R $outdir/tmp.txt $outdir/tmp2.txt locus_name - --manyto1 --tokeep allx | awk -F$'\t' '{OFS=FS} NR!=1 {print $1,$2,$3,$4,"0",$5}' > TAIR10_introns_noChrCM.bed
	
	bwfiles_all=""; bwfiles_mat=""; bwfiles_pat=""; j=1; k=1
	for ((i=0;i<${#stubname[@]};++i)); do
		echo "Processing ${j}th sample..."
		bwfiles_all="${bwfiles_all} $outdir/bamCoverage/bamCoverage_dedup_all/${stubname[i]}_cov.bw"
		
		if [ $(( $j % 100 )) -eq 0 ]; then
			echo "Running computeMatrix..."
			cmd="computeMatrix scale-regions -S $bwfiles_all -R $outdir/TAIR10_plus_araport11_nonoverlapping_noChrCM.bed -a 500 -b 500 -m 2000 -bs 50 -o $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_part${k}.txt"
			bsub -o "$outdir/bamCoverage/lsf_logs/computeMatrix_batch${k}.txt" -K "$cmd" & pids+=( $! ); logfiles+=( "$outdir/bamCoverage/lsf_logs/computeMatrix_batch${k}.txt" )
			cmd="computeMatrix scale-regions -S $bwfiles_all -R TAIR10_introns_noChrCM.bed -a 200 -b 200 -m 200 -bs 10 -o $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_part${k}.txt"
			bsub -o "$outdir/bamCoverage/lsf_logs/computeMatrix_introns_batch${k}.txt" -K "$cmd" & pids+=( $! ); logfiles+=( "$outdir/bamCoverage/lsf_logs/computeMatrix_introns_batch${k}.txt" )

			k=$(( $k + 1 ))
			bwfiles_all=""; bwfiles_mat=""; bwfiles_pat=""
		fi

		j=$(( $j + 1 ))
	done
	
	# submit last job
	echo "Running computeMatrix..."
	cmd="computeMatrix scale-regions -S $bwfiles_all -R $outdir/TAIR10_plus_araport11_nonoverlapping_noChrCM.bed -a 500 -b 500 -m 2000 -bs 50 -o $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_part${k}.txt"
	bsub -o "$outdir/bamCoverage/lsf_logs/computeMatrix_batch${k}.txt" -K "$cmd" & pids+=( $! ); logfiles+=( "$outdir/bamCoverage/lsf_logs/computeMatrix_batch${k}.txt" )
	cmd="computeMatrix scale-regions -S $bwfiles_all -R TAIR10_introns_noChrCM.bed -a 200 -b 200 -m 200 -bs 10 -o $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_part${k}.txt"
	bsub -o "$outdir/bamCoverage/lsf_logs/computeMatrix_introns_batch${k}.txt" -K "$cmd" & pids+=( $! ); logfiles+=( "$outdir/bamCoverage/lsf_logs/computeMatrix_introns_batch${k}.txt" )
	
	# wait for all jobs to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait ${pids[i]} || err_msg "computeMatrix failed, see log file ${logfiles[i]}" "$log"
	done
	
	# repeat for the Col/Cvi RIL samples
	bwfiles=""
	for samp in col_1 col_2 col_3 cvi_1 cvi_2 cvi_3; do
		bwfiles="${bwfiles} $outdir/bamCoverage/${samp}_RILs_cov.bw"
	done		
	
	computeMatrix scale-regions -S $bwfiles -R $outdir/TAIR10_plus_araport11_nonoverlapping_noChrCM.bed -a 500 -b 500 -m 2000 -bs 50 -o $outdir/bamCoverage/bamCoverage_dedup_all_RIL_samples.txt
	computeMatrix scale-regions -S $bwfiles -R TAIR10_introns_noChrCM.bed -a 200 -b 200 -m 200 -bs 10 -o $outdir/bamCoverage/bamCoverage_introns_dedup_all_RIL_samples.txt
	
	# now run plotprofile on each matrix, make it output final values also
	plotProfile --matrixFile $outdir/bamCoverage/bamCoverage_dedup_all_RIL_samples.txt --outFileName $outdir/bamCoverage/bamCoverage_dedup_all_RIL_samples.pdf --outFileNameData $outdir/bamCoverage/bamCoverage_dedup_all_RIL_samples_metaplot.txt --perGroup 
	plotProfile --matrixFile $outdir/bamCoverage/bamCoverage_introns_dedup_all_RIL_samples.txt --outFileName $outdir/bamCoverage/bamCoverage_introns_dedup_all_RIL_samples.pdf --outFileNameData $outdir/bamCoverage/bamCoverage_introns_dedup_all_RIL_samples_metaplot.txt --perGroup 

	for idx in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do
		fflist1="${fflist1} $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_part${idx}.txt"
		fflist2="${fflist2} $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_part${idx}.txt"
		echo "Running plotProfile on idx $idx"
		plotProfile --matrixFile $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_part${idx}.txt --outFileName $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_part${idx}.pdf --outFileNameData $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_part${idx}_metaplot.txt --perGroup
		plotProfile --matrixFile $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_part${idx}.txt --outFileName $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_part${idx}.pdf --outFileNameData $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_part${idx}_metaplot.txt --perGroup
	done
		
	# also make version with avg. and sd. of each column only
	rm $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_allsamples_metaplot.txt
	for ((i=1;i<=15;++i)); do
		tail -n+3 $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_part${i}_metaplot.txt >> $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_allsamples_metaplot.txt
	done
	printf "bin\tsnRNAseq_mean\tsnRNAseq_sd\n" > $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_allsamples_mean_sd.txt
	cat $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_allsamples_metaplot.txt | awk -F$'\t' '{OFS=FS} {for (i=3;i<=NF;i++){S[i]+=$i; SS[i]+=($i)^2}} END {for (i in S){print i-2,S[i]/NR,sqrt(SS[i]/NR - (S[i]/NR)^2)}}' >> $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_allsamples_mean_sd.txt

	# repeat for the RIL samples
	printf "bin\tRIL_mean\tRIL_sd\n" > $outdir/bamCoverage/bamCoverage_dedup_all_RIL_samples_mean_sd.txt
	tail -n+3 $outdir/bamCoverage/bamCoverage_dedup_all_RIL_samples_metaplot.txt | awk -F$'\t' '{OFS=FS} {for (i=3;i<=NF;i++){S[i]+=$i; SS[i]+=($i)^2}} END {for (i in S){print i-2,S[i]/NR,sqrt(SS[i]/NR - (S[i]/NR)^2)}}' >> $outdir/bamCoverage/bamCoverage_dedup_all_RIL_samples_mean_sd.txt
	
	paste $outdir/bamCoverage/bamCoverage_dedup_all_snRNAseq_allsamples_mean_sd.txt <( cut -f2-3 $outdir/bamCoverage/bamCoverage_dedup_all_RIL_samples_mean_sd.txt ) > $outdir/bamCoverage/bamCoverage_dedup_all_samples_mean_sd.txt
	
	# repeat for introns
	rm $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_allsamples_metaplot.txt
	for ((i=1;i<=15;++i)); do
		tail -n+3 $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_part${i}_metaplot.txt >> $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_allsamples_metaplot.txt
	done
	printf "bin\tsnRNAseq_mean\tsnRNAseq_sd\n" > $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_allsamples_mean_sd.txt
	cat $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_allsamples_metaplot.txt | awk -F$'\t' '{OFS=FS} {for (i=3;i<=NF;i++){S[i]+=$i; SS[i]+=($i)^2}} END {for (i in S){print i-2,S[i]/NR,sqrt(SS[i]/NR - (S[i]/NR)^2)}}' >> $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_allsamples_mean_sd.txt

	# repeat for the RIL samples
	printf "bin\tRIL_mean\tRIL_sd\n" > $outdir/bamCoverage/bamCoverage_introns_dedup_all_RIL_samples_mean_sd.txt
	tail -n+3 $outdir/bamCoverage/bamCoverage_introns_dedup_all_RIL_samples_metaplot.txt | awk -F$'\t' '{OFS=FS} {for (i=3;i<=NF;i++){S[i]+=$i; SS[i]+=($i)^2}} END {for (i in S){print i-2,S[i]/NR,sqrt(SS[i]/NR - (S[i]/NR)^2)}}' >> $outdir/bamCoverage/bamCoverage_introns_dedup_all_RIL_samples_mean_sd.txt
	
	paste $outdir/bamCoverage/bamCoverage_introns_dedup_all_snRNAseq_allsamples_mean_sd.txt <( cut -f2-3 $outdir/bamCoverage/bamCoverage_introns_dedup_all_RIL_samples_mean_sd.txt ) > $outdir/bamCoverage/bamCoverage_introns_dedup_all_samples_mean_sd.txt
		

	# for Table S2, get summary of significant endo/seedcoat enrichment
	# ----------------
	printf "locus_name\tendo\tseedcoat\n" > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_factor2_pscores_fxheader.txt
	tail -n+2 $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_factor2_pscores.txt >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_factor2_pscores_fxheader.txt
		
	printf "locus_name\tDEGG\n" > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_heatmap_k60_factor2_pscores_roworder_fxheader_noNAs.txt
	awk -F$'\t' '$2 != "NA"' $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_heatmap_k60_factor2_pscores_roworder.txt >> $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_heatmap_k60_factor2_pscores_roworder_fxheader_noNAs.txt
	
	$path_to_scripts/merge_by_column.R $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_factor2_pscores_fxheader.txt $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/cluster_DEG_kmeans60_heatmap_k60_factor2_pscores_roworder_fxheader_noNAs.txt locus_name - | awk -F$'\t' '{OFS=FS} NR!=1 {ss_e[$4]+=$2; ss_s[$4]+=$3; tt[$4]+=1} END {for (a in ss_e) {print a,ss_e[a]/tt[a],ss_s[a]/tt[a]}}' > $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/endo_seedcoat_avg_biases.txt
	
	cat $outdir/expression_variability/expr_var_over_clusters/all_nuc_control_cross_tissue/cluster_DEG/endo_seedcoat_avg_biases.txt
	
	
	# make plots of expression over the clusters with Becky's probes (Fig. 2a and Ext. data Fig. 6)
	# ----------------
	# version with max at 12
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/misc_plots/probe_expr_plots/all_probes_CxV_endo_4DAP_max12" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Becky_insitu_genes_list.txt" --yorder AT4G11080,AT1G09380,AT5G10440,AT1G44090,AT2G44240,AT4G13380 --sizeupper 1 --dotsize 15 --fillupper 12 --xorder 4,3,8,10,11,9,12,13,14,1,7,2,6,5 --outputdatawide
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/misc_plots/probe_expr_plots/all_probes_CxV_seedcoat_4DAP_max12" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Becky_insitu_genes_list.txt" --yorder AT4G11080,AT1G09380,AT5G10440,AT1G44090,AT2G44240,AT4G13380 --sizeupper 1 --dotsize 15 --fillupper 12 --xorder 2,1,5,4,3,6 --outputdatawide
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/misc_plots/probe_expr_plots/all_probes_VxC_endo_4DAP_max12" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Becky_insitu_genes_list.txt" --yorder AT4G11080,AT1G09380,AT5G10440,AT1G44090,AT2G44240,AT4G13380 --sizeupper 1 --dotsize 15 --fillupper 12 --xorder 2,4,3,1,6,8,10,9,7,11,5 --outputdatawide
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/misc_plots/probe_expr_plots/all_probes_VxC_seedcoat_4DAP_max12" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Becky_insitu_genes_list.txt" --yorder AT4G11080,AT1G09380,AT5G10440,AT1G44090,AT2G44240,AT4G13380 --sizeupper 1 --dotsize 15 --fillupper 12 --xorder 2,1,5,4,7,8,6,3 --outputdatawide

	# for the sup probes, decrease color scale upper limit a bit to 10	
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/misc_plots/probe_expr_plots/sup_probes_CxV_endo_4DAP" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Becky_insitu_genes_list_sup.txt" --yorder AT1G44090,AT4G13380 --sizeupper 1 --fillupper 10 --dotsize 15 --xorder 4,3,8,10,11,9,12,13,14,1,7,2,6,5 --outputdatawide
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/misc_plots/probe_expr_plots/sup_probes_CxV_seedcoat_4DAP" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Becky_insitu_genes_list_sup.txt" --yorder AT1G44090,AT4G13380 --sizeupper 1 --fillupper 10 --dotsize 15 --xorder 2,1,5,4,3,6 --outputdatawide
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/misc_plots/probe_expr_plots/sup_probes_VxC_endo_4DAP" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Becky_insitu_genes_list_sup.txt" --yorder AT1G44090,AT4G13380 --sizeupper 1 --fillupper 10 --dotsize 15 --xorder 2,4,3,1,6,8,10,9,7,11,5 --outputdatawide
	$path_to_scripts/single_cell_RNAseq_plots.R 'dot' "$outdir/misc_plots/probe_expr_plots/sup_probes_VxC_seedcoat_4DAP" --CPM "$outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt" --sampfile <( tail -n+2 "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters.txt" ) --genefile "/lab/solexa_gehring/colette/single_cell_seq/Becky_insitu_genes_list_sup.txt" --yorder AT1G44090,AT4G13380 --sizeupper 1 --fillupper 10 --dotsize 15 --xorder 2,1,5,4,7,8,6,3 --outputdatawide
	
	
	# make plots of total expression over seed coat related genes (Ext. data Fig. 3 and Fig. S3)
	# ----------------
	mkdir -p $outdir/misc_plots/tot_expr/seed_coat
	mkdir -p $outdir/misc_plots/tot_expr/endo
	
	# MYB5 (epidermis-specific, see https://www.sciencedirect.com/science/article/pii/S0012160608012505)
	gene="AT3G13540"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 12 --includezeros --plotwidth 8
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 12 --includezeros --plotwidth 9
	
	# GL2 (outermost epidermal layer, see refs in http://www.plantphysiol.org/content/139/2/701, Windsor JB, Symonds VV, Mendenhall J, Lloyd AM (2000) Arabidopsis seed coat development: morphological differentiation of the outer integument. Plant J 22: 483–493)
	gene="AT1G79840"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 8 --includezeros --plotwidth 8
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 8 --includezeros --plotwidth 9

	# SCL15 chalazal seed coat (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4507008/)
	gene="AT4G36710"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 3 --includezeros --plotwidth 8
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 3 --includezeros --plotwidth 9
	
	# TT1 (endothelium marker - https://www.frontiersin.org/articles/10.3389/fpls.2019.01801/full)
	gene="AT1G34790"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 3.5 --includezeros --plotwidth 8
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 3.5 --includezeros --plotwidth 9
	
	# TT12 (endothelium marker - http://www.plantcell.org/content/13/4/853)
	gene="AT3G59030"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 10 --includezeros --plotwidth 8
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 10 --includezeros --plotwidth 9

	# TT10 (endothelium marker - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1276023/)
	gene="AT5G48100"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 5.5 --includezeros --plotwidth 8
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 5.5 --includezeros --plotwidth 9

	# deltaVPE (ii2 and ii3 - http://www.plantcell.org/content/17/3/876)
	gene="AT3G20210"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 8 --ymin 9 --ymax 15
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 9 --ymin 9 --ymax 15

	# BAN (ii1 - https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004856)
	gene="AT1G61720"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 8 --ymin 5 --ymax 13
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 9 --ymin 5 --ymax 13

	# TT16 (ii1 and ii2 - https://dev.biologists.org/content/144/8/1490)
	gene="AT5G23260"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 8
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 9

	# AGL69 (chalazal seed coat - https://www.pnas.org/content/pnas/107/18/8063.full.pdf)
	gene="AT5G65070"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 2 --includezeros --plotwidth 8
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --ymin 0 --ymax 2 --includezeros --plotwidth 9

	# SUC2 (chalazal seed coat - http://www.plantphysiol.org/content/139/2/701)
	gene="AT1G22710"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 8 --ymax 0.5
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 9 --ymax 0.5

	# ATT1 (endothelium - https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-313X.2007.03348.x)
	gene="AT4G00360"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 8 --ymin 0 --ymax 10
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 9 --ymin 0 --ymax 10

	# outer seed coat
	gene="AT3G15510"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 8 --ymin 0 --ymax 6
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 9 --ymin 0 --ymax 6
	
	
	# plots of total expression in endosperm
	# ----------------
	# chromatin-related genes in endosperm (Ext. data Fig. 9)
	gene="AT4G11400"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_CxV_v2 --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 13 --noprelog --ymax 275
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_VxC_v2 --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 10.5 --noprelog --ymax 275

	gene="AT4G34060"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_CxV_v2 --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 13 --noprelog --ymax 250
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_VxC_v2 --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 10.5 --noprelog --ymax 250

	gene="AT5G52230"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_CxV_v2 --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 13 --noprelog --ymax 600
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_VxC_v2 --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 10.5 --noprelog --ymax 600

	gene="AT4G31900"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_CxV_v2 --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 13 --noprelog --ymax 600
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_VxC_v2 --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 10.5 --noprelog --ymax 600


	# a few example MEGs comparing seed coat and endo expression (Fig. S13)
	gene="AT1G75900"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 8 --noprelog
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 9 --noprelog
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 13 --noprelog
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 10.5 --noprelog
	
	gene="AT3G23060"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 8 --noprelog
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/seed_coat/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 9 --noprelog
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_CxV --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 13 --noprelog
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/${gene}_VxC --CPM $outdir/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes $gene --dotsize 30 --sizeupper 1 --includezeros --plotwidth 10.5 --noprelog


	# plots of imprinted expression in endosperm (also includes plots of total expression) - Fig. 3
	mkdir -p $outdir/misc_plots/impr_expr

	gene="AT3G60350"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CxV_E --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 14 --includezeros --xorder 12,13,14,10,11,9,6,5,1,7,4,3,2,8 --ymin 0 --ymax 300
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_VxC_E --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 11 --includezeros --xorder 1,6,4,3,9,7,8,10,2,5,11 --ymin 0 --ymax 300
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CxV_E_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 14 --sizeupper 1 --xorder 12,13,14,10,11,9,6,5,1,7,4,3,2,8 --ymax 10
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_VxC_E_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 11 --sizeupper 1 --xorder 1,6,4,3,9,7,8,10,2,5,11 --ymax 10

	gene="AT5G10950"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CxV_E --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 14 --includezeros --xorder 12,13,14,10,11,9,6,5,1,7,4,3,2,8 --ymin 0 --ymax 120
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_VxC_E --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 11 --includezeros --xorder 1,6,4,3,9,7,8,10,2,5,11 --ymin 0 --ymax 120
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CxV_E_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/clustering/SC3_subsets/CxV_endo_4DAP_final_clusters_noheader.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 14 --sizeupper 1 --xorder 12,13,14,10,11,9,6,5,1,7,4,3,2,8 --ymax 8
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_VxC_E_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/clustering/SC3_subsets/VxC_endo_4DAP_final_clusters_noheader.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 11 --sizeupper 1 --xorder 1,6,4,3,9,7,8,10,2,5,11 --ymax 8

	# and over cell cycle - Fig. 3
	gene="AT1G66200"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CC --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 8 --xorder G0,G1,G1toS,S,G2,M --ymin 0 --includezeros
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CC_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M --ymin 0

	gene="AT2G32720"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CC --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 8 --xorder G0,G1,G1toS,S,G2,M --ymin 0 --includezeros
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CC_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M --ymin 0

	gene="AT2G45470"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CC --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 8 --xorder G0,G1,G1toS,S,G2,M --ymin 0 --includezeros
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CC_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M --ymin 0

	gene="AT1G78570"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CC --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 8 --xorder G0,G1,G1toS,S,G2,M --ymin 0 --includezeros
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/${gene}_CC_AS --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M --ymin 0


	# plots over cell cycle (Fig. S16)
	# ----------------
	gene="AT4G12870"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/CC_${gene} --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 8 --xorder G0,G1,G1toS,S,G2,M --ymin 0 --ymax 3500
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/CC_line_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M --ymin 0 --ymax 200
	$path_to_scripts/single_cell_RNAseq_plots.R bar $outdir/misc_plots/impr_expr/CC_bar_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --xorder G0,G1,G1toS,S,G2,M
	
	gene="AT3G53420"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/CC_${gene} --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 8 --xorder G0,G1,G1toS,S,G2,M --ymin 0
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/CC_line_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M --ymin 0
	$path_to_scripts/single_cell_RNAseq_plots.R bar $outdir/misc_plots/impr_expr/CC_bar_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --xorder G0,G1,G1toS,S,G2,M
	
	gene="AT5G48650"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/CC_${gene} --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 8 --xorder G0,G1,G1toS,S,G2,M --ymin 0
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/CC_line_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M --ymin 0
	$path_to_scripts/single_cell_RNAseq_plots.R bar $outdir/misc_plots/impr_expr/CC_bar_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --xorder G0,G1,G1toS,S,G2,M
	
	gene="AT1G57820"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/CC_${gene} --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 8 --xorder G0,G1,G1toS,S,G2,M --ymin 0
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/CC_line_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M --ymin 0
	$path_to_scripts/single_cell_RNAseq_plots.R bar $outdir/misc_plots/impr_expr/CC_bar_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --xorder G0,G1,G1toS,S,G2,M
	
	gene="AT4G02060"
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/tot_expr/endo/CC_${gene} --CPM snRNAseq_analysis/count_matrices/all_cpm_in_genes_dedup_passQC.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes $gene --noprelog --dotsize 30 --sizeupper 1 --plotwidth 8 --xorder G0,G1,G1toS,S,G2,M --ymin 0
	$path_to_scripts/single_cell_RNAseq_plots.R lin $outdir/misc_plots/impr_expr/CC_line_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --dotsize 30 --addsims --noprelog --plotwidth 8 --sizeupper 1 --xorder G0,G1,G1toS,S,G2,M --ymin 0
	$path_to_scripts/single_cell_RNAseq_plots.R bar $outdir/misc_plots/impr_expr/CC_bar_${gene} --mcounts $outdir/imprinting_analysis/all_mcounts_merged.txt --pcounts $outdir/imprinting_analysis/all_pcounts_merged.txt --sampfile "$outdir/expression_variability/nuclei_lists/endo_4DAP_nuclei_cell_cycle.txt" --genes ${gene} --xorder G0,G1,G1toS,S,G2,M
	

	# also look at cell cycle breakdown over those nuclei (which remember are probably mostly 4C) - Ext. data Fig. 3c
	mkdir -p $outdir/misc_plots/cell_cycle/seed_coat
	awk -F$'\t' '{OFS=FS} NR==FNR {a[$1]=$2; next} $1 in a {print $0,a[$1]}' $outdir/clustering/SC3_subsets/CxV_seedcoat_4DAP_final_clusters_noheader.txt $outdir/expression_variability/nuclei_lists/all_4DAP_nuclei_cell_cycle.txt | cut -f2,3 | sort | uniq -c | sed 's/^ *//g' | tr ' ' '\t' > $outdir/misc_plots/cell_cycle/seed_coat/CxV_seedcoat_4DAP_CC.txt
	awk -F$'\t' '{OFS=FS} NR==FNR {a[$1]=$2; next} $1 in a {print $0,a[$1]}' $outdir/clustering/SC3_subsets/VxC_seedcoat_4DAP_final_clusters_noheader.txt $outdir/expression_variability/nuclei_lists/all_4DAP_nuclei_cell_cycle.txt | cut -f2,3 | sort | uniq -c | sed 's/^ *//g' | tr ' ' '\t' > $outdir/misc_plots/cell_cycle/seed_coat/VxC_seedcoat_4DAP_CC.txt
	
	$path_to_scripts/barchart.R $outdir/misc_plots/cell_cycle/seed_coat/CxV_seedcoat_4DAP_CC.txt $outdir/misc_plots/cell_cycle/seed_coat/CxV_seedcoat_4DAP_CC.pdf --xvals 3 --yvals 1 --factor 2 --stack --flevels G0,G1,G1toS,S,G2,M
	$path_to_scripts/barchart.R $outdir/misc_plots/cell_cycle/seed_coat/CxV_seedcoat_4DAP_CC.txt $outdir/misc_plots/cell_cycle/seed_coat/CxV_seedcoat_4DAP_CC_pstack.pdf --xvals 3 --yvals 1 --factor 2 --pstack --flevels G0,G1,G1toS,S,G2,M
	$path_to_scripts/barchart.R $outdir/misc_plots/cell_cycle/seed_coat/VxC_seedcoat_4DAP_CC.txt $outdir/misc_plots/cell_cycle/seed_coat/VxC_seedcoat_4DAP_CC.pdf --xvals 3 --yvals 1 --factor 2 --stack --flevels G0,G1,G1toS,S,G2,M
	$path_to_scripts/barchart.R $outdir/misc_plots/cell_cycle/seed_coat/VxC_seedcoat_4DAP_CC.txt $outdir/misc_plots/cell_cycle/seed_coat/VxC_seedcoat_4DAP_CC_pstack.pdf --xvals 3 --yvals 1 --factor 2 --pstack --flevels G0,G1,G1toS,S,G2,M
	
fi




