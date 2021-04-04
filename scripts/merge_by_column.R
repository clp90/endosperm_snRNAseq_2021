#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# version 1.2 (12/15/2016)
# Given two input files and a column name or list of column names, merges the two
# files together according to that name. One-to-many merges are allowed, in which one
# file has more than one observation with the same value(s) in the mergeByList. If
# doing this, must specify the manyto1 option below, and file 1 must be the one with
# multiple obs per value(s) in the mergeByList. Many-to-many merges are not allowed
# (and never a good idea, anyway).

# Note that for obvious reasons, both input files must have headers containing
# variable names, and these must match the names provided in mergeByList. Files
# must also be tab-delimited.

# Also note that the many-to-one merge, but not the normal (one-to-one) merge requires
# the R package plyr. The one-to-one merge will generally be much faster.

# Usage: merge_by_column.R infile1 infile2 mergeByList outfile method

# Example: merge_by_column.R infile1.txt infile2.txt Alyr107_gene,TAIR10_gene imprinting_summary.txt

# -------------------------
# Version history:
# v.1.0: initial build - allows one-to-one merge over one or more variables
# v.1.1: (06/02/2015) added ability to do a many-to-one merge, added more helpful usage info, and added options
# to specify which entries to keep - all, all in file 1, all in file 2, or merged only (see below)
# v.1.2: (12/15/2016) added option to output directly to command line by specifying output file as "-" (similar to samtools)
# to allow for piping results directly for processing.

options(scipen=999)

option_list <- list(
	make_option("--tokeep", default="all",
		help = "Which entries to keep, must be 'all' (keep all obs in both files regardless of merge), 'merged' (only keep merged obs), 'allx' (keep all obs in file 1), or 'ally' (keep all obs in file 2) - default 'all'"),
	make_option("--manyto1", default=FALSE, action="store_true",
		help = "Perform a many-to-one merge: file1 contains more than one obs per unique value(s) in mergeByList, while file2 has only one obs per unique combination of mergeByList"),
	make_option("--NAto0", default=FALSE, action="store_true",
		help = "Convert all NAs into 0s in the merged output")
	)

args <- commandArgs(trailingOnly = TRUE)	
if (length(args) <= 3) {
	cat("Usage: merge_by_column.R [options] infile1 infile2 mergeByList outfile\n")
	cat("----------------------\n")
	cat("infile1 is any tab-delimited file containing variable names found in the first line (header)\n")
	cat("infile2 is any tab-delimited file containing variable names found in the first line, to be merged with infile1\n")
	cat("mergeByList is a comma-seperated (no spaces!) list of variable names found in both infile1 and infile2, to perform the merge over\n")
	cat("outfile is the name for the resulting merged output file\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--tokeep <'all','merged','allx','ally'> : which observations to keep after the merge (all obs, merged only, all obs in file 1, or all obs in file 2)\n")
	cat("--manyto1 : Perform a many-to-one merge: use if one file contains more than one obs per unique value(s) in mergeByList, while other input file has only one obs per unique combination of mergeByList\n")
	cat("--NAto0 : Convert all NAs into 0s in the merged output\n")
	cat("----------------------\n")
	stop("At least one required argument not provided. See usage above.")
}

arguments = parse_args(OptionParser(option_list=option_list), positional_arguments = 4)
opt = arguments$options
infile1 = arguments$args[1]
infile2 = arguments$args[2]
mergeByList = strsplit(arguments$args[3],',')[[1]]
outfile = arguments$args[4]

# summary of calls to this script
if (outfile != "-")	{
	cat("----------------------\n")
	cat("Merging the following two files according to",arguments$args[3],":\n")
	cat(infile1,"\n")
	cat(infile2,"\n")
	cat("----------------------\n")
	cat("Output file:", outfile,"\n")
	if ( opt$tokeep == "all" ) {
		cat("Keeping all obs regardless of merge status\n")
	} else if ( opt$tokeep == "merged" ) {
		cat("Keeping only obs where fields in mergeByList match in both input files\n")
	} else if ( opt$tokeep == "allx" ) {
		cat("Keeping all obs in",infile1,"regardless of merge status\n")
	} else if ( opt$tokeep == "ally" ) {
		cat("Keeping all obs in",infile2,"regardless of merge status\n")
	} else {
		stop("value '",opt$tokeep,"' is not valid for option --tokeep, allowed values are 'all','allx','ally','merged'")
	}
}

# now read in each file, merge 
alldata = read.table(infile1, header=TRUE, sep="\t",row.names=NULL)
temp = read.table(infile2, header=TRUE, sep="\t",row.names=NULL)

# check both data frames contain the column(s) specified by user
for (i in mergeByList) {
	if (! i %in% colnames(alldata)) {
		stop("Error: no column named \"",i,"\" in input file 1: ",infile1)
	}
	if (! i %in% colnames(temp)) {
		stop("Error: no column named \"",i,"\" in input file 2: ",infile2)
	}
}

# determine whether either input file has more than one obs per combination of mergeByList
numUnique1=unique(alldata[mergeByList])
numUnique2=unique(temp[mergeByList])

if ( length(numUnique1[,1]) != length(alldata[,1]) & length(numUnique2[,1]) != length(temp[,1]) ) {
		stop("Error: both input files have multiple observations per value of varlist: ",mergeByList,", which makes merging behavior ambiguous and is not allowed")
} 

if ( opt$manyto1 == TRUE ) {
	if ( length(numUnique1[,1]) == length(alldata[,1]) & length(numUnique2[,1]) == length(temp[,1]) ) {
		cat("Warning: both input files are unique over varlist:",mergeByList,", --manyto1 option not required\n")
	}
} else {
	if ( length(numUnique1[,1]) != length(alldata[,1]) | length(numUnique2[,1]) != length(temp[,1]) ) {
		stop("Error: one of the input files has non-unique values over varlist: ",mergeByList,", user must specify --manyto1 option")
	}
}

# merge together according to specified vars
if ( opt$manyto1 == TRUE ) {
	suppressPackageStartupMessages(library(plyr))
	
	if ( opt$tokeep == "all" ) {
		alldata = join(alldata,temp, by=mergeByList, type="full", match="all")
	} else if ( opt$tokeep == "merged" ) {
		alldata = join(alldata,temp, by=mergeByList, type="inner", match="all")
	} else if ( opt$tokeep == "allx" ) {
		alldata = join(alldata,temp, by=mergeByList, type="left", match="all")
	} else if ( opt$tokeep == "ally" ) {
		alldata = join(alldata,temp, by=mergeByList, type="right", match="all")
	}
} else {
	if ( opt$tokeep == "all" ) {
		alldata = merge(alldata,temp,by=mergeByList, all = TRUE)
	} else if ( opt$tokeep == "merged" ) {
		alldata = merge(alldata,temp,by=mergeByList)
	} else if ( opt$tokeep == "allx" ) {
		alldata = merge(alldata,temp,by=mergeByList, all.x = TRUE)
	} else if ( opt$tokeep == "ally" ) {
		alldata = merge(alldata,temp,by=mergeByList, all.y = TRUE)
	}
}

# print data frame to file
if ( opt$NAto0 == TRUE )	{
	alldata[is.na(alldata)] <- 0
}

if (outfile == "-")	{
	write.table(alldata, "", sep = "\t", quote = FALSE, row.names = FALSE, na="")
} else {
	write.table(alldata, outfile, sep = "\t", quote = FALSE, row.names = FALSE, na="")
}












