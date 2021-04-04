#!/bin/bash

# ------------------------------------------------------------------------------------
# v1.0 by Colette L. Picard
# 11/17/2018
# ------------------------------------------------------------------------------------

# Usage:
# rna_seq_map.sh [options] -1 forward.fq -g genome -o outdir

# -------------------------
# Version history:
# v.1.0: initial build - 11/17/2018
# -------------------------

# Description printed when "help" option specified:
read -d '' usage <<"EOF"
v.1.0 by Colette L Picard, 11/17/2018
-------------------------------
Usage:
rna_seq_map.sh -1 forward.fq -g genome -o outdir [options]
-------------------------------

At its most basic, this is a simple script for the initial processing of RNA-seq data. 
It performs read QC (using fastqc), filtering (using trim_galore/cutadapt) and mapping 
(using STAR).

This version works for both paired-end and single-end data; if providing single-end
data only -1 will be used, if providing paired-end data both -1 and -2 will be used.

OVERVIEW:
(1) Run fastqc on input .fq files to assess read quality
(2) Filter reads based on quality using trim_galore
(3) Align reads to reference genome (or metagenome) using STAR
(4) Remove PCR duplicates and clean up
(5) Summarize mapping results

-------------------------------
ADDITIONAL CONSIDERATIONS:
- Raw read files are provided using -1 (and optionally -2, if paired-end reads) and
should be in the FASTQ format.
- The genome to map reads to is provided with -g and must first be prepped using STAR
in --runMode genomeGenerate (see below for example syntax)
- The location for output files is provided with -o and must be a non-existant or empty
folder, else the script will (by default) refuse to overwrite its contents. Note this can be overridden
with the -r option, which you use at your own risk.
- All other parameters are optional (see list below), but consider setting -n, and looking
at the FASTQC report to set -a and -t. 
- Double-check your FASTQ quality encoding (will be in FASTQC output) and set -3 if
encoding is PHRED+33 instead of PHRED+64 (see https://en.wikipedia.org/wiki/FASTQ_format)
- This script uses a java program (MarkDuplicates.jar, from the picard-tools suite) to remove
PCR duplicates; by default java is given access to 10gb of data; this value may need to be
changed depending on your own setup (use -m).
-------------------------------

Example syntax for generating STAR-ready genome from a genome.fa and set of annotations
in GTF format (see https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf for more info):

STAR --runMode genomeGenerate --genomeDir folder/for/output/files --genomeFastaFiles genome.fa --sjdbGTFfile GTF.gtf --sjdbOverhang 30

NOTE1: genomeGenerate run above, genome annotations should be provided if you have them
(the GTF file). Note that STAR will prep genomes without the GTF, and this can be used as well.
NOTE2: if your genome is really small, --genomeSAindexNbases will have to be set to a smaller
number than the default. See manual for more details.
NOTE3: STAR is used here in --alignEndsType EndToEnd - e.g. soft-clipping of reads is disabled.

EXTRA STUFF:
There are a few additional things this script can handle, beyond the basic analysis noted above:
-------------------------------
(1) The "metagenome" approach

This option is designed to minimize mapping bias in favor of the sequenced 
strain/accession when mapping reads from hybrid tissues. This approach uses a "metagenome" 
containing the genomes of both species or accessions, appended together, so that reads are 
mapped to both genomes simultaneously. 

For example, this could be used if mapping A. thaliana reads from a hybrid of the accessions 
Col and Cvi. Normally, when mapping to the A. thaliana genome (which is derived from Col), 
Col reads will map better than Cvi reads because they will match the provided genome better, 
leading to bias in favor of Col. For the "metagenome" approach, Cvi SNPs are substituted 
into the sequenced genome (Col), and this Cvi "pseudogenome" is appended to the regular genome.
Chromosome names are altered so that Col- and Cvi-chromosomes can be distinguished. Reads 
are then mapped to this "metagenome" instead of the regular Col-derived genome sequence.

WARNING: metagenome mode can only be used with two genomes with identical coordinates! That means
SNPs are ok, but no indels. This is because tools for converting mapped read coordinates
from one genome to another (e.g. LiftOver) are not splice-aware.

In this mode, the user provides a (STAR-prepped) "metagenome" with -g and a mapping of the metagenome 
chromosomes to the "regular" chromosomes with -C. If -C is not provided, the provided genome
is assumed to be a normal genome, not a metagenome. Initial QC and mapping proceeds as normal, after which
reads that mapped uniquely to the metagenome are extracted as potentially allele-specific. 
Ambiguously mapped reads are further filtered to identify those which map exactly once to each of the
two original genomes (in the same location) - these are considered uniquely mapped, but
not specific to either of the two original genomes. Note that if you choose to also provide
a GTF file with -G, you will need to make a "meta-GTF" file with the same chromosomes as the metagenome
instead of using the regular GTF file.

An example of the metagenome chromosomes "mapping" file -C:
chr	Col	Cvi
Chr1	Chr1C	Chr1V
Chr2	Chr2C	Chr2V
Chr3	Chr3C	Chr3V
Chr4	Chr4C	Chr4V
Chr5	Chr5C	Chr5V
ChrC	ChrCC	ChrCV
ChrM	ChrMC	ChrMV

-------------------------------
(2) Support for mapping to spike-in genomes

If using RNA spike-ins like ERCC (ThermoFisher), users will also need to map reads to the 
spike-in "genome". This script allows users to simultaneously map normal and spike-in reads.
To use this feature, append the spike-in genome to the genome being used and run STAR genomeGenerate
on this spike-in appended genome. Provide this genome normally using -g, and use -P 
provide a list of the chromosome names corresponding to all chromosomes/sequences in the 
spike-in genome (the names after the '>' on alternating lines in the fasta file). After 
mapping, any reads mapping to any of the chromosomes in the list provided with -P will be
separated from the other reads into a new file and counted. If providing a GTF file, can append
the spike-in GTF to the regular GTF file. Note that this can also be used with the metagenome
approach (just append the spike-in chromosomes to the metagenome, and provide both -P and -C).
-------------------------------

Required and optional arguments are listed below. Default values, when applicable, are given
at the end of the line in square brackets. Required programs and scripts also listed.

User-specified options:
Required arguments:
	-1 forward : name of file containing single-end reads or forward paired-end reads (FASTQ format, must have extension *.txt, *.fastq, or *.fq)
	-g genome : folder containing STAR-indexed reference genome for alignment (prepped using STAR in --genomeGenerate mode)
	-o outdir : name of directory with all output files; will be created if nonexistant
Required arguments for metagenome approach:
	-C metachrom : list of chromosomes and their names for both strains (3 columns: e.g. "Chr1	ChrC1	ChrV1" where ChrC1 is Col, ChrV1 is Cvi chr1, with header e.g. "chr	Col	Cvi")
Additional options:
	-2 reverse : name of file containing reverse reads if using paired-end data (FASTQ format, must have extension *.txt, *.fastq, or *.fq) []
	-n name : stubname used to name output files (e.g. "$n"_pairs.sam, etc.) ["expt"]
	-t trim : (trim_galore) trims t bases immediately from 5' end of forward (and reverse reads, if applicable) during QC analysis (see trim_galore option -t) [0]
	-q minQ : (trim_galore) quality to use for 3' trimming of reads (see option -q in trim_galore manual) [25]
	-a adapter : (trim_galore) sequence of adapter for adapter trimming (default is universal Illumina adapter) ["AGATCGGAAGAG"]
	-i minIL : (STAR) minimum intron length (see STAR manual, --alignIntronMin option) [70]
	-I maxIL : (STAR) maximum intron length (see STAR manual, --alignIntronMax option) [5000]
	-d dist : (STAR) max allowed distance between mate pairs (see STAR --alignMatesGapMax option) [100000]
	-N maxmismatch : (STAR) max number of mismatches, as a fraction of read length (for PE reads, read1+read2 length) (see STAR -outFilterMismatchNoverReadLmax option) [0.05]
	-A name1 : (metagenome only) name of strain A, e.g. corresponding to column 2 in file from -C ["strain1"]
	-B name2 : (metagenome only) name of strain B, e.g. corresponding to column 3 in file from -C ["strain2"]
	-P spikeins : file containing list of chromosome names corresponding to spike-ins; any reads mapped to chromosomes listed here will be split into own file (append spike-in .fa file to main genome or metagenome) []
	-s path_to_scripts : path to folder containing all required scripts (note - $scriptDir in default is the location of this script) ["$scriptDir/scripts"]
	-m memalloc : memory (in Gb) to allocate to java when running MarkDuplicates [10]
	-T threads : number of threads that STAR will use [1]
Flag options:
	-F : force allow reads from the two input files to not be same length (e.g if already filtered; usually should not be used) [forcelen=false]
	-S : library is stranded (default unstranded) [stranded=false]
	-3 : quality encoding starts at ASCII char 33 (!) instead of char 64 (@) [phred33=false]
	-r : allow overwrite of existing files in provided outdir - default FALSE (WARNING: all files currently in outdir will be deleted!) [overwrite=false]
	-0 : checks that all required programs installed on PATH and all required helper scripts can be located, then exits without running
	-h : prints this version and usage information
	
Required installed on user PATH:
	- fastqc (from Babraham Institute, tested on v0.4.1)
	- trim_galore (from Babraham Institute, tested on v0.4.1)
	- cutadapt (tested on v.1.8 by Marcel Martin)
	- STAR (by Dobin et al. 2013, tested on v2.6.1d)
	- samtools (by Heng Li, tested on v1.7)
	- python (tested on v2.7.6)
	- java (tested on v.1.7.0_181, OpenJDK Runtime Environment (IcedTea 2.6.14) (7u181-2.6.14-0ubuntu0.2), OpenJDK 64-Bit Server VM (build 24.181-b01, mixed mode))
	
Must be in path_to_scripts (if not on default path, specify location with -s):
	- get_uniq_from_meta.py (by Colette L Picard, 17 May 2018)
	- MarkDuplicates.jar (from the picard-tools suite, Broad Institute, downloaded Oct. 2, 2015)

------------------------------------------------------------------------------------
EOF

[[ $# -eq 0 ]] && { printf "%s\n" "$usage"; exit 0; } 		# if no user-supplied arguments, print usage and exit

# ----------------------
# Get user-specified arguments
# ----------------------

# Initiate environment
scriptDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )	# location of this script 
workdir=$( pwd )											# working directory

# Required arguments:
# ----------------------
forward=""							# name of file containing single-end reads or forward paired-end reads (FASTQ format, must have extension *.txt, *.fastq, or *.fq)
genome=""							# reference genome for alignment (prepped using STAR in genomeGenerate mode)
outdir=""							# name of directory with all output files; will be created if nonexistant

# Additional options:
# ----------------------
reverse=""							# name of file containing reverse reads if using paired-end data (FASTQ format, must have extension *.txt, *.fastq, or *.fq)
name="expt"							# stubname used to name output files (e.g. "$n"_pairs.sam, etc.)
trim=0								# (trim_galore) trims t bases immediately from 5' end of both forward and reverse reads during QC analysis (see trim_galore -t)
minQ=25								# (trim_galore) quality to use for 3' trimming of reads (see option -q in trim_galore manual)
adapter="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"							# (trim_galore) sequence adapter for adapter trimming
minIL=70								# (STAR) minimum intron length (see STAR manual, --alignIntronMin)
maxIL=5000								# (STAR) maximum intron length (see STAR manual, --alignIntronMax)
dist=100000								# (STAR) max allowed distance between mate pairs (see STAR --alignMatesGapMax option)
maxmismatch=0.05								# (STAR) max number of mismatches, as a fraction of total read length (of both pairs together if PE) (see STAR --outFilterMismatchNoverReadLmax option)
name1="strain1"							# (metagenome only) name of reference strain, e.g. corresponding to column 2 in file from -C
name2="strain2"							# (metagenome only) name of alternate strain, e.g. corresponding to column 3 in file from -C
spikeins=""							# file containing list of chromosome names corresponding to spike-ins; any reads mapped to chromosomes listed here will be split into own file (append spike-in .fa file to main genome or metagenome)
path_to_scripts="$scriptDir/scripts"							# path to folder containing all required scripts (note - $scriptDir in default is the location of this script)
memalloc=10							# memory (in Gb) to allocate to java when running MarkDuplicates
threads=1						# number of threads that STAR will use

# Flag options:
# ----------------------
stranded=false							# library is stranded (false by default)
forcelen=false							# force allow reads from the two input files to not be same length (e.g if already filtered; usually should not be used)
phred33=false							# quality encoding starts at ASCII char 33 (!) instead of char 64 (@)
overwrite=false							# allow overwrite of existing files in provided outdir - default FALSE (WARNING: all files currently in outdir will be deleted!)

# Required arguments for metagenome approach:
# ----------------------
metachrom=""							# list of chromosomes and their names for both strains (3 columns: e.g. "Chr1	ChrC1	ChrV1" where ChrC1 is Col, ChrV1 is Cvi chr1, with header e.g. "chr	Col	Cvi")

checkdep=false

# ----------------------
while getopts "1:g:o:C:2:n:t:q:a:G:l:i:I:d:N:A:B:P:s:m:T:F3Sr0h" opt; do
	case $opt in
		1)	# name of file containing single-end reads or forward paired-end reads (FASTQ format, must have extension *.txt, *.fastq, or *.fq)
			forward="$OPTARG"
			;;
		g)	# reference genome for alignment (prepped using STAR in genomeGenerate mode)
			genome="$OPTARG"
			;;
		o)	# name of directory with all output files; will be created if nonexistant
			outdir="$OPTARG"
			;;
		C)	# list of chromosomes and their names for both strains (3 columns: e.g. "Chr1	ChrC1	ChrV1" where ChrC1 is Col, ChrV1 is Cvi chr1, with header e.g. "chr	Col	Cvi")
			metachrom="$OPTARG"
			;;
		2)	# name of file containing reverse reads if using paired-end data (FASTQ format, must have extension *.txt, *.fastq, or *.fq)
			reverse="$OPTARG"
			;;
		n)	# stubname used to name output files (e.g. "$n"_pairs.sam, etc.)
			name="$OPTARG"
			;;
		t)	# (trim_galore) trims t bases immediately from 5' end of both forward and reverse reads during QC analysis (see trim_galore -t)
			trim="$OPTARG"
			;;
		q)	# (trim_galore) quality to use for 3' trimming of reads (see option -q in trim_galore manual)
			minQ="$OPTARG"
			;;
		a)	# (trim_galore) sequence of forward adapter for adapter trimming
			adapter="$OPTARG"
			;;
		i)	# (STAR) minimum intron length (see STAR manual, --alignIntronMin option)
			minIL="$OPTARG"
			;;
		I)	# (STAR) maximum intron length (see STAR manual, --alignIntronMax option)
			maxIL="$OPTARG"
			;;
		d)	# (STAR) max allowed distance between mate pairs (see STAR -alignMatesGapMax option)
			dist="$OPTARG"
			;;
		N)	# (STAR) max number of mismatches as a fraction of read length
			maxmismatch="$OPTARG"
			;;
		A)	# (metagenome only) name of reference strain, e.g. corresponding to column 2 in file from -C
			name1="$OPTARG"
			;;
		B)	# (metagenome only) name of alternate strain, e.g. corresponding to column 3 in file from -C
			name2="$OPTARG"
			;;
		P)	# file containing list of chromosome names corresponding to spike-ins; any reads mapped to chromosomes listed here will be split into own file (append spike-in .fa file to main genome or metagenome)
			spikeins="$OPTARG"
			;;
		s)	# path to folder containing all required scripts (note - $scriptDir in default is the location of this script)
			path_to_scripts="$OPTARG"
			;;
		m)	# memory (in Gb) to allocate to java when running MarkDuplicates
			memalloc="$OPTARG"
			;;
		T)	# number of threads that STAR will use
			threads="$OPTARG"
			;;
		F)	# force allow reads from the two input files to not be same length (e.g if already filtered; usually should not be used)
			forcelen=true
			;;
		3)	# quality encoding starts at ASCII char 33 (!) instead of char 64 (@)
			phred33=true
			;;
		S)	# library is stranded
			stranded=true
			;;
		r)	# allow overwrite of existing files in provided outdir - default FALSE (WARNING: all files currently in outdir will be deleted!)
			overwrite=true
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

# ----------------------
# Helper functions for this script
# ----------------------
err_msg ()
# prints an error message both to stdout and to the log file, then exits
# usage: err_msg msg logfile
{
	printf "Error: $1 \n"
	printf "Exited due to error: $1 \n" >> "$2"
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

compress_sam ()
# sorts and compresses a SAM file to BAM and then indexes it; removes the original SAM file
# Usage: compress_sam infile scriptlog
{
	infile="$1"; logfile="$2"
	# grab filename minus extension
	ff="${infile%.*}"
	if [[ -f "$infile" && -s "$infile" ]]; then
		samtools sort "$infile" -o "${ff}.bam" -T "${ff}" > /dev/null 2>&1
		[ $? != 0 ] && err_msg "samtools sort failed for input file $infile" "$logfile"
		samtools index "${ff}.bam"
		[ $? != 0 ] && err_msg "samtools index failed for input file ${ff}.bam" "$logfile"
		rm "$infile"
	fi
}


# ----------------------
# Main code
# ----------------------

# Check that all programs required on PATH are installed
# ----------------------
command -v trim_galore >/dev/null 2>&1 || { echo "Error: trim_galore is required on PATH but was not found"; exit 1; }
command -v STAR >/dev/null 2>&1 || { echo "Error: STAR is required on PATH but was not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools is required on PATH but was not found"; exit 1; }
command -v cutadapt >/dev/null 2>&1 || { echo "Error: cutadapt is required on PATH but was not found"; exit 1; }
command -v python >/dev/null 2>&1 || { echo "Error: python2 is required on PATH but was not found"; exit 1; }
command -v java >/dev/null 2>&1 || { echo "Error: java is required on PATH but was not found"; exit 1; }


# Check that all other required helper scripts can be located
# ----------------------
[ -f "$path_to_scripts/get_uniq_from_meta.py" ] || { echo "Error: required helper script get_uniq_from_meta.py was not found in expected location $path_to_scripts"; exit 1; }
[ -f "$path_to_scripts/MarkDuplicates.jar" ] || { echo "Error: required helper script MarkDuplicates.jar was not found in expected location $path_to_scripts"; exit 1; }

# Check that MarkDuplicates is working
# ----------------------
res=$( java -jar "$path_to_scripts/MarkDuplicates.jar" 2>&1 | head -1 )
[ "$res" = "ERROR: Option 'OUTPUT' is required." ] || { echo "Could not run MarkDuplicates from the picard-tools suite, check if it's compatible with your java installation"; exit 1; }

# Check that all required python modules are installed
# ----------------------
python -c "import argparse"; [ $? != 0 ] && { echo "Error: required python module argparse not installed"; exit 1; }

# Done checking all requirements. Stop here if -0 flagged.
# ----------------------
"$checkdep" && exit 0

# If output folder doesn't exist, make it, else overwrite if user used -r, else error
# ----------------------
if [ ! -d "$outdir" ]; then
	mkdir "$outdir"
else
	if [ "$(ls -A $outdir)" ]; then
		if [ "$overwrite" = "true" ]; then
			echo "Overwriting previous contents of output dir $outdir"
			rm -rf "$outdir"/*	
		else
			echo "Error: provided output directory is not empty. To allow overwrite of non-empty dir, use -r flag. WARNING: all existing files in -outdir- will be deleted. Seriously."
			exit 1
		fi
	fi
fi

# Check all required inputs are provided
# ----------------------
[ -z "$forward" ] && { echo "Error: -1 forward is a required argument (name of file containing single-end reads or forward paired-end reads (FASTQ format, must have extension *.txt, *.fastq, or *.fq))"; exit 1; }
[ -z "$genome" ] && { echo "Error: -g genome is a required argument (reference genome for alignment (prepped using STAR in genomeGenerate mode))"; exit 1; }
[ -z "$outdir" ] && { echo "Error: -o outdir is a required argument (name of directory with all output files; will be created if nonexistant)"; exit 1; }

# Check that all required inputs exist
# ----------------------
[ -f "$forward" ] || { echo "Could not open forward reads file $forward"; exit 1; }
[[ ! -z "$reverse" && ! -f "$reverse" ]] && { echo "Could not open reverse reads file $reverse"; exit 1; }
[ -f "$genome/chrLength.txt" ] || { echo "Could not open STAR genome file $genome/chrLength.txt"; exit 1; }
[ -f "$genome/chrName.txt" ] || { echo "Could not open STAR genome file $genome/chrName.txt"; exit 1; }
[ -f "$genome/chrNameLength.txt" ] || { echo "Could not open STAR genome file $genome/chrNameLength.txt"; exit 1; }

# Check library type was provided
# ----------------------
[[ "$firststrand" = "true" && "$secondstrand" = "true" ]] && { echo "Error: flags -f and -c (library strandedness) are mutually exclusive"; exit 1; }
[[ "$firststrand" = "false" && "$secondstrand" = "false" ]] && libtype="fr-unstranded"
[ "$firststrand" = "true" ] && libtype="fr-firststrand"
[ "$secondstrand" = "true" ] && libtype="fr-secondstrand"

# Check optional input files exist if provided
# ----------------------
[[ ! -z "$metachrom" && ! -f "$metachrom" ]] && { echo "Error: could not open metachrom file $metachrom"; exit 1; }
[[ ! -z "$gtf" && ! -f "$gtf" ]] && { echo "Error: could not open GTF file $gtf"; exit 1; }
[[ ! -z "$spikeins" && ! -f "$spikeins" ]] && { echo "Could not open spike-in chromosome file $spikeins"; exit 1; }



log="$outdir/${name}_log.txt" 	# create log file
time_start=$(date)				# time run was started
time_ss=$(date +%s)				# time run was started (in seconds)

# Output user-derived options to stdout and to log file
# ----------------------
echo "Running rna_seq_map.sh v1.0 (11/17/2018):" | tee "$log"
echo "Run start on: $time_start" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Working directory: $( pwd )" | tee -a "$log"
[ ! -z "$reverse" ] || echo "Reads file: $forward" | tee -a "$log"
[ ! -z "$reverse" ] && echo "Forward reads file: $forward" | tee -a "$log"
[ ! -z "$reverse" ] && echo "Reverse reads file: $reverse" | tee -a "$log"
echo "Output directory: $outdir" | tee -a "$log"
echo "-------------------------" | tee -a "$log"

if [ ! -z "$metachrom" ]; then
	echo "Mapping reads to the $name1 + $name2 metagenome:" | tee -a "$log"
	echo "Metagenome: $genome" | tee -a "$log"
	echo "File with chromosome meta info: $metachrom" | tee -a "$log"
else
	echo "Mapping reads to: $genome" | tee -a "$log"
fi
[ ! -z "$spikeins" ] && echo "Spike-in genome included in main genome, list of chromosomes provided in: $spikeins" | tee -a "$log"

echo "-------------------------" | tee -a "$log"
echo "Additional settings:" | tee -a "$log"
echo "Stubname: $name" | tee -a "$log"
echo "Helper scripts are in: $path_to_scripts" | tee -a "$log"
[ "$phred33" = "true" ] && echo "Quality score encoding is: phred+33" | tee -a "$log" || echo "Quality score encoding is: phred+64" | tee -a "$log"
echo "trim_galore settings:" | tee -a "$log"
echo "   # bases to trim at 5' end of forward and reverse reads: $trim" | tee -a "$log"
echo "   Quality cutoff for trimming at 3' end: $minQ" | tee -a "$log"
echo "   Adapter sequence to trim: $adapter" | tee -a "$log"
echo "STAR settings:" | tee -a "$log"
echo "   Number of threads to use: $threads" | tee -a "$log"
echo "   Library is stranded: $stranded" | tee -a "$log"
echo "   Minimum intron length: $minIL" | tee -a "$log"
echo "   Maximum intron length: $maxIL" | tee -a "$log"
echo "   Max allowed distance between mates (if paired-end): $dist" | tee -a "$log"
echo "   Max number of mismatches as a fraction of total read length: $maxmismatch" | tee -a "$log"
echo "-------------------------" | tee -a "$log"


# ----------------------
# Step 0: check input files ok
# ----------------------
printf "\nStep 0: Checking that all inputs are ok\n" | tee -a "$log"

# Get all chromosomes in provided file
# ----------------------
allchrs=$( cat "$genome/chrName.txt" | tr '\n' ' ' | tr '\r' ' ' | tr '\r\n' ' ' )
allchrsarray=( $allchrs )
numchrs="${#allchrsarray[@]}"
echo "$numchrs chromosomes detected in genome" | tee -a "$log"

# If spike-in chromosome list provided, convert it into space-separated string and ensure those chrs present in genome
# ----------------------
if [ ! -z "$spikeins" ]; then
	spikeinchrs=$( cat "$spikeins" | tr '\n' ' ' | tr '\r' ' ' | tr '\r\n' ' ' )
	spikeinchrsarray=( $spikeinchrs )
	for ((i = 0; i < ${#spikeinchrsarray[@]}; ++i)); do
		[[ "$allchrs" =~ ${spikeinchrsarray[i]} ]] || { echo "Error: at least one chromosome in $spikeins not found in genome $genome"; exit 1; }
	done
	echo " - ${#spikeinchrsarray[@]} chromosomes are from the spike-in genome" | tee -a "$log"
fi

# If metagenome is being used, check that all chromosomes in cols 2 and 3 of metachrom file can be located
# ----------------------
if [ ! -z "$metachrom" ]; then
	# get info from the metachrom file
	normchrs=$(cut -f1 "$metachrom" | tail -n +2); str1chrs=$(cut -f2 "$metachrom" | tail -n +2); str2chrs=$(cut -f3 "$metachrom" | tail -n +2)
	normchrsarray=($normchrs); str1chrsarray=($str1chrs); str2chrsarray=($str2chrs)	
	for ((i = 0; i < ${#str1chrsarray[@]}; ++i)); do
		[[ "$allchrs" =~ ${str1chrsarray[i]} ]] || { echo "Error: at least one chromosome in -C metachrom file in column 2 (${name1}) not found in genome $genome"; exit 1; }
	done
	for ((i = 0; i < ${#str2chrsarray[@]}; ++i)); do
		[[ "$allchrs" =~ ${str2chrsarray[i]} ]] || { echo "Error: at least one chromosome in -C metachrom file in column 3 (${name2}) not found in genome $genome"; exit 1; }
	done
	echo " - ${#str1chrsarray[@]} chromosomes each are from $name1 and $name2 (metagenome approach)" | tee -a "$log"
fi

# If GTF file was provided, check that all chromosomes in the GTF file are in genome (not all chromosomes 
# in genome must appear in GTF - e.g. if genome contains spike-in chrs but GTF doesn't, that's fine)
# ----------------------
if [ ! -z "$gtf" ]; then
	# in GTF format, chromosome is given by first column
	gtfchrs=$( cut -f1 "$gtf" | sort -u | tr '\n' ' ' )
	gtfchrsarray=( $gtfchrs )
	for ((i = 0; i < ${#gtfchrsarray[@]}; ++i)); do
		[[ "$allchrs" =~ ${gtfchrsarray[i]} ]] || { echo "Error: chromosome ${gtfchrsarray[i]} appears in GTF file $gtf but not in genome $genome"; exit 1; }
	done
fi

# Check that input FASTQ file(s)s have correct extensions
# ----------------------
f_basename=$( basename $forward ); ext_f="${f_basename#*.}"
f_dirname=$(dirname $forward)
[ ! -z "$reverse" ] && { r_basename=$( basename $reverse ); ext_r="${r_basename#*.}"; }
[ ! -z "$reverse" ] && r_dirname=$(dirname $reverse)
[[ ! -z "$reverse" && "$ext_f" != "$ext_r" ]] && err_msg "Input files $forward and $reverse do not have the same file extension; note this error may also occur if you have dots (.) in your filenames." "$log"

if [ "$ext_f" != "txt" ] && [ "$ext_f" != "fq" ] && [ "$ext_f" != "fastq" ]; then
	err_msg "Input file(s) -1, -2 have unrecognized extension $ext_f , permitted extensions are .txt, .fq, and .fastq (must be FASTQ format); note this error may also occur if you have dots (.) in your filenames." "$log"
fi

# Get initial number of reads in input file(s); if paired-end then check that # reads same in both files
# ----------------------
f_tot=$( wc -l "$forward" | awk '{print $1}' )								# of lines in forward file
[ ! -z "$reverse" ] && r_tot=$( wc -l "$reverse" | awk '{print $1}' )		# of lines in reverse file

[ $(( $f_tot % 4 )) -ne 0 ] && err_msg "number of lines in -1 fastq file ($f_tot) is not a multiple of four (expected for fastq file)." "$log"
[[ ! -z "$reverse" && $(( $r_tot % 4 )) -ne 0 ]] && err_msg "number of lines in -2 fastq file ($r_tot) is not a multiple of four (expected for fastq file)." "$log"

# Get # of reads by dividing by 4, check these are same
f_tot=$(( $f_tot/4 ))
[ ! -z "$reverse" ] && r_tot=$(( $r_tot/4 ))
[[ ! -z "$reverse" && "$f_tot" -ne "$r_tot" ]] && err_msg "-1 and -2 files do not contain same number of reads.\n # Forward reads = $f_tot,  # Reverse reads = $r_tot" "$log"

# also get read lengths (assume all reads same length, so use length of read in first entry)
readf=$( head -n2 $forward | tail -1 ); readlenf=$( echo ${#readf} )
[ ! -z "$reverse" ] && { readr=$( head -n2 $reverse | tail -1 ); readlenr=$( echo ${#readr} ); }
if [ "$forcelen" = "false" ]; then
	[[ ! -z "$reverse" && "$readlenf" -ne "$readlenr" ]] && err_msg "Reads in -1 and -2 are not the same length (based on first read only); override with -F" "$log"
fi
printf "\nStarting number of reads: %s" "$f_tot" | tee -a "$log"	
[ ! -z "$reverse" ] && printf "\nReads are ${readlenf}x${readlenr}bp\n" | tee -a "$log"	
[ ! -z "$reverse" ] || printf "\nReads are ${readlenf}bp\n" | tee -a "$log"	
echo "-------------------------" | tee -a "$log"


# ----------------------
# Step 1: run fastqc on the forward and reverse reads seperately
# ----------------------
ts=$(date +%s)
printf "\nStep 1: running fastqc to assess initial read quality\n" | tee -a "$log"
mkdir "$outdir/fastqc_prefiltering"

# -----------RUN FASTQC-----------
fastqc "$forward" -o "$outdir/fastqc_prefiltering" 2> /dev/null
[ ! -z "$reverse" ] && fastqc "$reverse" -o "$outdir/fastqc_prefiltering" 2> /dev/null

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"


# ----------------------
# Step 2: filter reads using trim_galore
# ----------------------
ts=$(date +%s)
printf "\nStep 2: filtering and trimming reads using trim_galore\n" | tee -a "$log"
mkdir "$outdir/qc_filtered"; mkdir "$outdir/fastqc_postfiltering"

[ "$phred33" = "true" ] && qualstr="phred33" || qualstr="phred64"
[ "$trim" = "0" ] && trimstr="" || { [ ! -z "$reverse" ] && trimstr=" --clip_R1 $trim --clip_R2 $trim" || trimstr=" --clip_R1 $trim"; }
[ ! -z "$reverse" ] && pairedstr=" --paired -a2 $adapter" || pairedstr=""
[ ! -z "$reverse" ] && infiles="$forward $reverse" || infiles="$forward"

# -----------RUN TRIM_GALORE-----------
trim_galore --fastqc --${qualstr} -a "$adapter" -q "$minQ" --stringency 3${trimstr}${pairedstr} -o "$outdir/qc_filtered" $infiles > "$outdir/qc_filtered/log_trim_galore.txt" 2>&1
[ $? != 0 ] && err_msg "error when running trim_galore, see $outdir/qc_filtered/log_trim_galore.txt" "$log"

fname=$(basename "$forward"); f_ext="${fname#*.}"
[ "$f_ext" = "fq" ] && fname="${fname%.*}"			# if extension of input files is .fq or fastq, trim_galore output will strip the .fq extension before adding _val_1.fq
[ "$f_ext" = "fastq" ] && fname="${fname%.*}"

# summarize results, move FASTQC output files to new folder, rename files for consistency
if [ ! -z "$reverse" ]; then		
	rname=$(basename "$reverse"); r_ext="${rname#*.}"
	[ "$r_ext" = "fq" ] && rname="${rname%.*}"
	[ "$r_ext" = "fastq" ] && rname="${rname%.*}"

	# Output summary of filtering results
	totQC=$( awk '/Total number/ {print $6}' "$outdir/qc_filtered/log_trim_galore.txt" )
	echo "Number of read pairs analyzed: $totQC" | tee -a "$log"
	removedQC=$( awk '/Number of sequence pairs removed/ {print $19}' "$outdir/qc_filtered/log_trim_galore.txt" )
	printf "Number pairs removed (one or both reads < 20bp): $removedQC (%0.1f%%)\n" "$(echo "scale=3; $removedQC / $totQC *100" | bc)" | tee -a "$log"
	remQC=$(( $totQC - $removedQC ))
	printf "Number pairs remaining: $remQC (%0.1f%%)\n" "$(echo "scale=3; $remQC / $totQC *100" | bc)" | tee -a "$log"
	
	[ "$remQC" = "0" ] && err_msg "No reads passed quality filtering step, please check your filtering parameters (in particular, are you using correct PHRED+ encoding?)" "$log"

	# rename output files for clarity
	mv "$outdir/qc_filtered/${fname}_val_1.fq" "$outdir/qc_filtered/${name}_f_pairs.fq"
	mv "$outdir/qc_filtered/${rname}_val_2.fq" "$outdir/qc_filtered/${name}_r_pairs.fq"
	
	# move fastqc files
	mv "$outdir/qc_filtered/${fname}_val_1_fastqc.html" "$outdir/fastqc_postfiltering/${name}_f_fastqc.html"
	mv "$outdir/qc_filtered/${rname}_val_2_fastqc.html" "$outdir/fastqc_postfiltering/${name}_r_fastqc.html"
	mv "$outdir/qc_filtered/${fname}_val_1_fastqc.zip" "$outdir/fastqc_postfiltering/${name}_f_fastqc.zip"
	mv "$outdir/qc_filtered/${rname}_val_2_fastqc.zip" "$outdir/fastqc_postfiltering/${name}_r_fastqc.zip"
else
	# Output summary of filtering results
	totQC=$( awk '/sequences processed in total/ {print $1}' "$outdir/qc_filtered/log_trim_galore.txt" )
	echo "Number of reads analyzed: $totQC" | tee -a "$log"
	removedQC=$( awk '/Sequences removed because they became shorter/ {print $14}' "$outdir/qc_filtered/log_trim_galore.txt" )
	printf "Number reads removed (read < 20bp): $removedQC (%0.1f%%)\n" "$(echo "scale=3; $removedQC / $totQC *100" | bc)" | tee -a "$log"
	remQC=$(( $totQC - $removedQC ))
	printf "Number reads remaining: $remQC (%0.1f%%)\n" "$(echo "scale=3; $remQC / $totQC *100" | bc)" | tee -a "$log"

	[ "$remQC" = "0" ] && err_msg "No reads passed quality filtering step, please check your filtering parameters (in particular, are you using correct PHRED+ encoding?)" "$log"

	mv $( printf '%s\n' $outdir/qc_filtered/${fname}*trimmed.fq ) "$outdir/qc_filtered/${name}_trimmed.fq"
	mv $( printf '%s\n' $outdir/qc_filtered/${fname}*trimming_report.txt ) "$outdir/qc_filtered/${name}_trim_report.txt"

	# move and rename output files
	mv "$outdir/qc_filtered/${fname}_trimmed_fastqc.html" "$outdir/fastqc_postfiltering/${name}_fastqc.html"
	mv "$outdir/qc_filtered/${fname}_trimmed_fastqc.zip" "$outdir/fastqc_postfiltering/${name}_fastqc.zip"
fi

[ "$remQC" = "0" ] && err_msg "No reads passed quality filtering step, please check your filtering parameters (in particular, are you using correct PHRED+ encoding?)" "$log"

te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	


# ----------------------
# Step 3: align filtered reads to the provided reference genome or metagenome
# ----------------------
ts=$(date +%s)
[ -z "$metachrom" ] && printf "\nStep 3: aligning to provided reference genome using STAR\n" | tee -a "$log"
[ ! -z "$metachrom" ] && printf "\nStep 3: aligning reads to provided metagenome using STAR\n" | tee -a "$log"
mkdir "$outdir/STAR"

# get optional parameters to include in STAR command
[ "$stranded" = "true" ] && strstr="" || strstr=" --outSAMstrandField intronMotif"
[ ! -z "$reverse" ] && infiles="$outdir/qc_filtered/${name}_f_pairs.fq $outdir/qc_filtered/${name}_r_pairs.fq" || infiles="$outdir/qc_filtered/${name}_trimmed.fq"
[ ! -z "$metachrom" ] && multimap=2 || multimap=1

# -----------RUN STAR-----------
STAR --alignSoftClipAtReferenceEnds No --alignEndsType EndToEnd${strstr} --runThreadN "$threads" --genomeDir "$genome" --readFilesIn $infiles --outFilterType BySJout --outFilterMultimapNmax "$multimap" --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax "$maxmismatch" --alignIntronMin "$minIL" --alignIntronMax "$maxIL" --alignMatesGapMax "$dist" --outFileNamePrefix "$outdir/STAR/${name}" --outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonical > "$outdir/STAR/${name}_STAR_log.txt"
[ $? != 0 ] && err_msg "error when aligning reads using STAR, see $outdir/STAR/log_STAR.txt" "$log"
echo "Done mapping to provided genome with STAR" | tee -a "$log"

# compress FASTQ files since they are no longer needed
gzip $infiles

# delete temp folder created by STAR
rm -rf "$outdir/STAR/${name}_STARtmp"

# initial STAR alignments will be uniquely mapping if ! metagenome mode; if in metagenome mode, reads
# with exactly two alignments need to grabbed also. MAPQ = 255 for uniquely mapped reads, and 
# int(-10*log10(1-(1/Nmap)) for multi-mapping reads - so for exactly two best alignments, MAPQ==3

# note: STAR is outputting singleton alignments from paired-end data in some datasets at a low rate 
# (ex. 23813 out of 25579958 or 0.09%); this anecdotally may be more likely with longer sequencing reads.
# for now, these are dropped here b/c it's unclear that these are trustable, high-quality alignments
if [ ! -z "$reverse" ]; then
	# save singletons to new file
	samtools view -F 4 -f 8 -bo "$outdir/STAR/${name}_discarded_singletons.bam" "$outdir/STAR/${name}Aligned.sortedByCoord.out.bam"
	num_single=$( samtools flagstat "$outdir/STAR/${name}_discarded_singletons.bam" | awk '$0 ~ /singletons/ {print $1}' )
	if [ "$num_single" -gt 0 ]; then
		# warn user and drop the singletons from the main alignments file before proceeding 
		echo "Warning: $num_single singleton alignments were detected from STAR output; these were censored but were saved to $outdir/STAR/${name}_discarded_singletons.bam" | tee -a "$log"
		samtools view -f 1 -F 12 -bo "$outdir/STAR/${name}Aligned.sortedByCoord.out.tmp.bam" "$outdir/STAR/${name}Aligned.sortedByCoord.out.bam"
		mv "$outdir/STAR/${name}Aligned.sortedByCoord.out.tmp.bam" "$outdir/STAR/${name}Aligned.sortedByCoord.out.bam"
	else
		rm "$outdir/STAR/${name}_discarded_singletons.bam"
	fi
fi

if [ ! -z "$metachrom" ]; then
	echo " - Extracting reads that mapped uniquely to the metagenome" | tee -a "$log"
	mv "$outdir/STAR/${name}Aligned.sortedByCoord.out.bam" "$outdir/STAR/${name}_all_alignments.bam"
	samtools index "$outdir/STAR/${name}_all_alignments.bam"
	samtools view -hq 255 "$outdir/STAR/${name}_all_alignments.bam" > "$outdir/STAR/${name}_unique_alignments.sam" 
	compress_sam "$outdir/STAR/${name}_unique_alignments.sam"
else
	mv "$outdir/STAR/${name}Aligned.sortedByCoord.out.bam" "$outdir/STAR/${name}_unique_alignments.bam"
	samtools index "$outdir/STAR/${name}_unique_alignments.bam"
fi

# if using spike-ins, extract those from the uniquely mapping reads here
# ----------------------
if [ ! -z "$spikeins" ]; then
	echo " - Extracting reads mapped to the spike-in chromosomes" | tee -a "$log"

	mv "$outdir/STAR/${name}_unique_alignments.bam" "$outdir/STAR/${name}_unique_alignments_with_spikeins.bam"
	mv "$outdir/STAR/${name}_unique_alignments.bam.bai" "$outdir/STAR/${name}_unique_alignments_with_spikeins.bam.bai"

	# get list of all chromosomes in the original $genome provided (will be same in header of all SAM/BAM files to date)
	samtools view -H "$outdir/STAR/${name}_all_alignments.bam" | awk '/^@SQ/' | cut -f2 | cut -f2 -d':' > "$outdir/fullchrlist.txt"
	
	# drop all chromosomes in $spikeins (list of spike-in chromosomes) to get all non-spike in chromosomes
	tokeep=$( diff "$outdir/fullchrlist.txt" "$spikeins" --unchanged-group-format='' | tr '\n' ' ' | tr '\r' ' ' | tr '\r\n' ' ' )

	# create new file with only the non-spike-in unique alignments
	samtools view -ho "$outdir/STAR/${name}_unique_alignments_tmp.sam" "$outdir/STAR/${name}_unique_alignments_with_spikeins.bam" $tokeep

	# drop all header lines with spike-in chromosomes and compress
	awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {todrop[$1]=""; next} {if ($1 !~ /^@SQ/) {print $0} else if ($1 ~ /^@SQ/ && !(substr($2,4) in todrop)) {print $0}}' "$spikeins" "$outdir/STAR/${name}_unique_alignments_tmp.sam" > "$outdir/STAR/${name}_unique_alignments.sam"
	rm "$outdir/STAR/${name}_unique_alignments_tmp.sam"
	compress_sam "$outdir/STAR/${name}_unique_alignments.sam" "$log"
	
	# create new file with only the uniquely mapping spike-in reads (non-uniquely mapping spike-in reads not counted)
	samtools view -ho "$outdir/STAR/${name}_unique_spikein_alignments_tmp.sam" "$outdir/STAR/${name}_unique_alignments_with_spikeins.bam" $spikeinchrs

	# drop all header lines with non-spike-in chromosomes and compress
	awk -F$'\t' 'BEGIN{OFS=FS} NR==FNR {tokeep[$1]=""; next} {if ($1 !~ /^@SQ/) {print $0} else if ($1 ~ /^@SQ/ && (substr($2,4) in tokeep)) {print $0}}' "$spikeins" "$outdir/STAR/${name}_unique_spikein_alignments_tmp.sam" > "$outdir/STAR/${name}_unique_spikein_alignments.sam"
	rm "$outdir/STAR/${name}_unique_spikein_alignments_tmp.sam"
	compress_sam "$outdir/STAR/${name}_unique_spikein_alignments.sam" "$log"
	
	# delete unneeded files
	rm "$outdir/fullchrlist.txt" "$outdir/STAR/${name}_unique_alignments_with_spikeins.bam" "$outdir/STAR/${name}_unique_alignments_with_spikeins.bam.bai"
fi

# if using the metagenome approach, get strain-specific and unique non-strain-specific alignments
# ----------------------
if [ ! -z "$metachrom" ]; then
	echo " - Determining strain-of-origin of reads that mapped uniquely to the metagenome" | tee -a "$log"		

	# get SAM file from BAM file of all uniquely mapped (non-spike in, if applicable) alignments
	samtools view -ho "$outdir/STAR/${name}_unique_alignments.sam" "$outdir/STAR/${name}_unique_alignments.bam"
	samtools view -ho "$outdir/STAR/${name}_all_alignments.sam" "$outdir/STAR/${name}_all_alignments.bam"
	
	# separate and tag (SAM tag XP:Z) reads from strain1 and strain2
	awk -F$'\t' -v strain="$name1" 'BEGIN{OFS=FS} NR==FNR {if (NR != 1) {pats[$2]=$1}; next} {if ($3 in pats) {$3=pats[$3]; if ($1 ~ /^@/) {print $0} else {print $0,"XP:Z:"strain}}; if (substr($2,4) in pats) {$2=("SN:"pats[substr($2,4)]); print $0}}' "$metachrom" "$outdir/STAR/${name}_unique_alignments.sam" > "$outdir/STAR/${name}_unique_alignments_to_${name1}.sam"
	awk -F$'\t' -v strain="$name2" 'BEGIN{OFS=FS} NR==FNR {if (NR != 1) {pats[$3]=$1}; next} {if ($3 in pats) {$3=pats[$3]; if ($1 ~ /^@/) {print $0} else {print $0,"XP:Z:"strain}}; if (substr($2,4) in pats) {$2=("SN:"pats[substr($2,4)]); print $0}}' "$metachrom" "$outdir/STAR/${name}_unique_alignments.sam" > "$outdir/STAR/${name}_unique_alignments_to_${name2}.sam"
	compress_sam "$outdir/STAR/${name}_unique_alignments_to_${name1}.sam" "$log"
	compress_sam "$outdir/STAR/${name}_unique_alignments_to_${name2}.sam" "$log"

	# also extract the reads that have exactly two best alignments (MAPQ == 3 with current version of STAR)
	echo " - Getting reads that mapped uniquely, but not specifically to either $name1 or $name2" | tee -a "$log"		
	cat "$outdir/STAR/${name}_all_alignments.sam" | awk '/^@/ || $5 == 3' | samtools sort - -n -o "$outdir/STAR/${name}_metagenome_twoalignments.sam" > /dev/null 2>&1
	
	# keep read pairs and singles that have exactly one unique alignment on each chromosome, in the same position, as "none"
	# (uniquely mapping but not to either genome in metagenome) - spike-ins filtered out by default since they aren't in the 'metachrom' file
	mkdir -p "$outdir/STAR/other_files_and_logs"
	
	$path_to_scripts/get_uniq_from_meta.py "$outdir/STAR/${name}_metagenome_twoalignments.sam" "$metachrom" "$outdir/STAR/${name}_unique_alignments_to_either.sam" > "$outdir/STAR/other_files_and_logs/get_unique_alignments_to_either.txt" 2>&1
	[ $? != 0 ] && err_msg "error when running get_uniq_from_meta.py, see $outdir/STAR/other_files_and_logs/get_unique_alignments_to_either.txt" "$log"
	compress_sam "$outdir/STAR/${name}_unique_alignments_to_either.sam" "$log"
	rm "$outdir/STAR/${name}_metagenome_twoalignments.sam" "$outdir/STAR/${name}_unique_alignments.sam" "$outdir/STAR/${name}_unique_alignments.bam" "$outdir/STAR/${name}_all_alignments.sam"
	
	# merge together all strain-specific and non-strain-specific "uniquely mapped" reads
	samtools merge "$outdir/STAR/${name}_unique_alignments.bam" "$outdir/STAR/${name}_unique_alignments_to_either.bam" "$outdir/STAR/${name}_unique_alignments_to_${name1}.bam" "$outdir/STAR/${name}_unique_alignments_to_${name2}.bam"
	samtools index "$outdir/STAR/${name}_unique_alignments.bam"
fi

te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	

# ----------------------
# Step 4: removing PCR duplicates + cleanup
# ----------------------
ts=$(date +%s)
printf "\nStep 4: creating files with PCR duplicates removed\n" | tee -a "$log"
echo "Note: all files with PCR duplicates removed will have suffix _dedup" | tee -a "$log"

# clean up some minor things
mkdir -p "$outdir/STAR/other_files_and_logs"
mv "$outdir/STAR/${name}_STAR_log.txt" "$outdir/STAR/other_files_and_logs"
mv "$outdir/STAR/${name}Log.final.out" "$outdir/STAR/other_files_and_logs/STAR_log_final_summary.txt"
mv "$outdir/STAR/${name}Log.out" "$outdir/STAR/other_files_and_logs/STAR_log.txt"
mv "$outdir/STAR/${name}Log.progress.out" "$outdir/STAR/other_files_and_logs/STAR_progress_log.txt"
mv "$outdir/STAR/${name}SJ.out.tab" "$outdir/STAR/other_files_and_logs"
[ ! -z "$metachrom" ] && mv "$outdir/STAR/${name}_all_alignments.bam" "$outdir/STAR/other_files_and_logs/${name}_all_STAR_alignments_unmodified.bam"
[ ! -z "$metachrom" ] && mv "$outdir/STAR/${name}_all_alignments.bam.bai" "$outdir/STAR/other_files_and_logs/${name}_all_STAR_alignments_unmodified.bam.bai"

# create new BAM files with PCR duplicates removed for all BAM files in $outdir/STAR
for ff in $outdir/STAR/*.bam; do
	bb="${ff%.*}"			# filename without extension
	base=$(basename "$ff")	# filename without path
	bbase="${base%.*}"		# filename without path or extension
	
	java -Xmx${memalloc}g -jar $path_to_scripts/MarkDuplicates.jar I="$ff" O="${bb}_dedup.bam" METRICS_FILE="$outdir/STAR/other_files_and_logs/MarkDuplicates_metrics_${bbase}.txt" REMOVE_DUPLICATES=true > "$outdir/STAR/other_files_and_logs/MarkDuplicates_log_${bbase}.txt" 2>&1
	[ $? != 0 ] && err_msg "MarkDuplicates failed, see $outdir/STAR/other_files_and_logs/MarkDuplicates_log_${bbase}.txt" "$log"
	
	samtools index "${bb}_dedup.bam"
done

te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	

# ----------------------
# Step 5: summarizing results
# ----------------------
ts=$(date +%s)
printf "\nStep 5: summarizing mapping results\n" | tee -a "$log"
mkdir "$outdir/STAR/summaries"
echo "" | tee -a "$log"	

for ff in $outdir/STAR/*.bam; do
	base=$(basename "$ff")	# filename without path
	bbase="${base%.*}"		# filename without path or extension
	samtools flagstat "$ff" > "$outdir/STAR/summaries/${bbase}_summary.txt"
done

[ -z "$reverse" ] && type="reads" || type="pairs"

echo "Input file contained $totQC $type" | tee -a "$log"	
printf "$remQC (%0.1f%%) $type passed quality filtering steps. Of these:\n" "$(echo "scale=3; $remQC / $totQC *100" | bc)" | tee -a "$log"	
for dd in "" "_dedup"; do
	[ "$dd" = "_dedup" ] && str="uniquely after PCR dedup" || str="uniquely"

	[ -z "$reverse" ] && uniq_all=$( awk '$0 ~ /mapped \(/ {print $1}' "$outdir/STAR/summaries/${name}_unique_alignments${dd}_summary.txt" 2>/dev/null )
	[ ! -z "$reverse" ] && uniq_all=$( awk '$0 ~ /with itself and mate mapped/ {print $1/2}' "$outdir/STAR/summaries/${name}_unique_alignments${dd}_summary.txt" 2>/dev/null )
	[ "$remQC" -gt 0 ] && frac_uniq=$( printf "%0.1f" $(echo "scale=3; $uniq_all / $remQC *100" | bc) ) || frac_spike="NA"

	if [ ! -z "$spikeins" ]; then
		[ -z "$reverse" ] && spikeins=$( awk '$0 ~ /mapped \(/ {print $1}' "$outdir/STAR/summaries/${name}_unique_spikein_alignments${dd}_summary.txt" 2>/dev/null )
		[ ! -z "$reverse" ] && spikeins=$( awk '$0 ~ /with itself and mate mapped/ {print $1/2}' "$outdir/STAR/summaries/${name}_unique_spikein_alignments${dd}_summary.txt" 2>/dev/null )
		uniq_all_wspike=$(( $uniq_all + $spikeins ))
		[ "$uniq_all_wspike" -gt 0 ] && frac_spike_mapped=$( printf "%0.1f" $(echo "scale=3; $spikeins / $uniq_all_wspike *100" | bc) ) || frac_spike=0
		[ "$remQC" -gt 0 ] && frac_spike_all=$( printf "%0.1f" $(echo "scale=3; $spikeins / $remQC *100" | bc) ) || frac_spike=0
		[ "$remQC" -gt 0 ] && frac_uniq=$( printf "%0.1f" $(echo "scale=3; $uniq_all_wspike / $remQC *100" | bc) ) || frac_spike=0
		
		printf "	$uniq_all_wspike (%0.1f%% of QC-filtered ${type}) ${type} mapped ${str}. Of these:\n" "$frac_uniq" | tee -a "$log"	
		printf "		$spikeins (%0.1f%% of uniquely mapped, %0.1f%% of QC-filtered) ${type} mapped ${str} to the spike-in genome\n" "$frac_spike_mapped" "$frac_spike_all" | tee -a "$log"	
	else		
		printf "	$uniq_all (%0.1f%% of QC-filtered ${type}) ${type} mapped $str" "$frac_uniq" | tee -a "$log"	
	fi

	if [ ! -z "$metachrom" ]; then
		[ -z "$spikeins" ] && { printf ". Of these:\n" | tee -a "$log"; }
		[ ! -z "$spikeins" ] && toadd=" non-spike-in" || toadd=""
		if [ -z "$reverse" ]; then
			strain1uniq=$( awk '$0 ~ /mapped \(/ {print $1}' "$outdir/STAR/summaries/${name}_unique_alignments_to_${name1}${dd}_summary.txt" 2>/dev/null )
			strain2uniq=$( awk '$0 ~ /mapped \(/ {print $1}' "$outdir/STAR/summaries/${name}_unique_alignments_to_${name2}${dd}_summary.txt" 2>/dev/null )
			eitheruniq=$( awk '$0 ~ /mapped \(/ {print $1}' "$outdir/STAR/summaries/${name}_unique_alignments_to_either${dd}_summary.txt" 2>/dev/null )
		else
			strain1uniq=$( awk '$0 ~ /with itself and mate mapped/ {print $1/2}' "$outdir/STAR/summaries/${name}_unique_alignments_to_${name1}${dd}_summary.txt" 2>/dev/null )
			strain2uniq=$( awk '$0 ~ /with itself and mate mapped/ {print $1/2}' "$outdir/STAR/summaries/${name}_unique_alignments_to_${name2}${dd}_summary.txt" 2>/dev/null )
			eitheruniq=$( awk '$0 ~ /with itself and mate mapped/ {print $1/2}' "$outdir/STAR/summaries/${name}_unique_alignments_to_either${dd}_summary.txt" 2>/dev/null )
		fi
		totallelic=$(( $strain1uniq + $strain2uniq ))
		[ "$uniq_all" -gt 0 ] && frac_str1=$( printf "%0.1f" $(echo "scale=3; $strain1uniq / $uniq_all *100" | bc) ) || frac_str1=0
		[ "$uniq_all" -gt 0 ] && frac_str2=$( printf "%0.1f" $(echo "scale=3; $strain2uniq / $uniq_all *100" | bc) ) || frac_str2=0
		[ "$uniq_all" -gt 0 ] && frac_both=$( printf "%0.1f" $(echo "scale=3; $eitheruniq / $uniq_all *100" | bc) ) || frac_both=0
		[ "$totallelic" -gt 0 ] && frac_str1_allelic=$( printf "%0.1f" $(echo "scale=3; $strain1uniq / $totallelic *100" | bc) ) || frac_str1=0
		[ "$uniq_all" -gt 0 ] && frac_str2_allelic=$( printf "%0.1f" $(echo "scale=3; $strain2uniq / $totallelic *100" | bc) ) || frac_str2=0
		printf "		$strain1uniq (%0.1f%% of uniquely mapped${toadd} ${type}, %0.1f%% of allelic ${type}) ${type} mapped ${str} to ${name1}\n" "$frac_str1" "$frac_str1_allelic" | tee -a "$log"	
		printf "		$strain2uniq (%0.1f%% of uniquely mapped${toadd} ${type}, %0.1f%% of allelic ${type}) ${type} mapped ${str} to ${name2}\n" "$frac_str2" "$frac_str2_allelic" | tee -a "$log"	
		printf "		$eitheruniq (%0.1f%% of uniquely mapped${toadd} ${type}) ${type} mapped ${str} to both strains\n" "$frac_both" | tee -a "$log"	
	else
		printf "\n" | tee -a "$log"		
	fi
done


te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	

time_end=$(date)	# time run was started
time_es=$(date +%s)	# time run was started
echo "" | tee -a "$log"
echo "Run ended $time_end" | tee -a "$log"
echo "Total time elapsed: $( displaytime $(($time_es - $time_ss)) )" | tee -a "$log"





