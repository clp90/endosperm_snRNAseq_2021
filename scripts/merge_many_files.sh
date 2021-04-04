#!/bin/bash

# ------------------------------------------------------------------------------------
# v1.0 by Colette L. Picard
# 06/23/2018
# ------------------------------------------------------------------------------------

# -------------------------
# Version history:
# v.1.0: initial build - 06/23/2018
# -------------------------

# Description printed when "help" option specified:
read -d '' usage <<"EOF"
v.1.0 by Colette L Picard, 06/23/2018

Simple wrapper script for merge_by_column.R to allow for efficient merging of large
numbers of files, when the file being created will be very large. Normally, when
just merging using merge_by_column.R in a loop over each individual file, the
merged file becomes so big that the I/O and the merging in R become very inefficient,
which is exacerbated when this is repeated many times.

This script instead takes a list of all the files you want to merge (provided in
the first input, the inputfilelist.txt), and splits them into separate batches to
be merged. Calls itself recursively until the number of input files is less than 20.
If number of input files is ≤ 20, this script simply merges them one after the other.
Speed-up is achieved by recursive calls to this script, submitted to LSF so they
can run simultaneously.

inputfilelist.txt = list of files to merge
outfile.txt = name for output file
mergebylist = comma-separated list of columns to merge over (see merge_by_column.R, same bhav)
howtomerge = one of three values (see merge_by_column.R also):
	- "all" = all values kept regardless of merge status
	- "allx" = all values from FIRST file kept regardless of merge status, values from other files
			   that didnt merge to first are discarded
	- "merged" = keep only values that merged across all files
	DEFAULT: "all"

Usage:
merge_many_files.sh [options] inputfilelist.txt outfile.txt mergebylist howtomerge

example:
merge_many_files.sh [options] inputfilelist.txt outfile.txt geneID merge
	
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
scriptDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )	# location of this script 
workdir=$( pwd )											# working directory

# Required arguments:
# ----------------------
[ "$#" -le 2 ] && { echo "Error: at least 3 arguments required (see usage)"; exit 1; }
inputfilelist="$1"
outfile="$2"
mergebylist="$3"
[ "$#" -eq 4 ] && howtomerge="$4" || howtomerge="all"	
[ "$howtomerge" != "all" ] && [ "$howtomerge" != "allx" ] && [ "$howtomerge" != "merged" ] && { echo "Error: value of 4th parameter must be either \"all\", \"allx\", or \"merged\""; exit 1; }

# Location of scripts (so that merge_by_column.R can be called)
# ----------------------
path_to_scripts="$scriptDir"	
[ ! -f "$path_to_scripts/merge_by_column.R" ] && { echo "Could not find script merge_by_column.R in $path_to_scripts"; exit 1; }


# Check all inputs provided
# ----------------------
[ -z "$inputfilelist" ] && { echo "Two arguments required (see usage)"; exit 1; }
[ -z "$outfile" ] && { echo "Two arguments required (see usage)"; exit 1; }
[ ! -f "$inputfilelist" ] && { echo "Could not open input file $inputfilelist"; exit 1; }
[ -f "$outfile" ] && { echo "Output file already exists; exiting"; exit 1; }

# Output basic info to user
# ----------------------
echo "Running merge_many_files.sh v1.0 (06/23/2018):"
echo "-------------------------"
echo "Working directory: $( pwd )"
echo "List of files to merge: $inputfilelist"
echo "Output file name: $outfile"
echo "Merging files by column(s) named: $mergebylist"
echo "Merging method: $howtomerge"
echo "-------------------------"


# Get all files to merge and check that they exist
# ----------------------
inputfilearray=()
while read ll; do
	[ -f "$ll" ] || { echo "Error: could not open input file $ll"; exit 1; }
	inputfilearray+=( $ll )
done < "$inputfilelist"
echo "${#inputfilearray[@]} files to be merged"


# Merge files:
# ----------------------
if [ $( wc -l "$inputfilelist" | awk '{print $1}' ) -le 20 ]; then
	# base case, merge all separately
	cat "${inputfilearray[0]}" > "$outfile"
	for ((i=1;i<${#inputfilearray[@]};++i)); do
		$path_to_scripts/merge_by_column.R "$outfile" "${inputfilearray[i]}" "$mergebylist" "$outfile" --tokeep "$howtomerge" > /dev/null
		[ $? -eq 0 ] || { echo "Error in merge_by_column.R"; echo "Failed command: $path_to_scripts/merge_by_column.R $outfile ${inputfilearray[i]} $mergebylist $outfile --tokeep $howtomerge"; exit 1; }
	done

else
	# split into sqrt(num_input_files) separate jobs and re-call merge_many_files.sh
	baseoutfile="${outfile%.*}"
	
	sqrt=$(echo "sqrt ( ${#inputfilearray[@]} )" | bc -l | cut -f1 -d '.')	
	[ "$sqrt" -le 20 ] && sqrt=20
	echo "Separating the ${#inputfilearray[@]} into $sqrt subsets to be merged separately, then merged back together"
	subsetlist=()
	for ((i=0;i<${#inputfilearray[@]};++i)); do
		subsetfilenum=$( echo "1+ ($i / $sqrt)" | bc -l | cut -f1 -d '.')
		[[ " ${subsetlist[*]} " == *" $subsetfilenum "* ]] || subsetlist+=( $subsetfilenum )
		echo "${inputfilearray[i]}" >> "${baseoutfile}_filelist_tmp${subsetfilenum}.txt"
	done
	
	# call merge_many_files on each subset, wait until finished then merge subsets
	pid=()
	for ((i=0;i<${#subsetlist[@]};++i)); do
		cmd="$scriptDir/merge_many_files.sh ${baseoutfile}_filelist_tmp${subsetlist[i]}.txt ${baseoutfile}_tmp${subsetlist[i]}.txt $mergebylist $howtomerge"
		bsub -o "${baseoutfile}_tmp${subsetlist[i]}_log.txt" -K "$cmd" & pid[i]=$!
	done
	
	for ((i=0;i<${#subsetlist[@]};++i)); do
		wait "${pid[i]}" || { echo "merge_many_files failed, see ${baseoutfile}_tmp${subsetlist[i]}_log.txt"; exit 1; }
	done
	
	# once all jobs completed, merge files back together
	echo "Done. Merging subset files back together..."
	echo "Merging in subset 1 out of ${#subsetlist[@]}"
	cat "${baseoutfile}_tmp${subsetlist[0]}.txt" > "$outfile"
	for ((i=1;i<${#subsetlist[@]};++i)); do
		echo "Merging in subset $(( $i + 1 )) out of ${#subsetlist[@]}"
		$path_to_scripts/merge_by_column.R "$outfile" "${baseoutfile}_tmp${subsetlist[i]}.txt" "$mergebylist" "$outfile" --tokeep "$howtomerge" > /dev/null
		[ $? -eq 0 ] || { echo "Error in merge_by_column.R, exiting"; exit 1; }
	done
	
	# remove all temp files
	for ((i=0;i<${#subsetlist[@]};++i)); do
		rm "${baseoutfile}_tmp${subsetlist[i]}_log.txt"
		rm "${baseoutfile}_tmp${subsetlist[i]}.txt"
		rm "${baseoutfile}_filelist_tmp${subsetlist[i]}.txt"
	done
	
fi
















