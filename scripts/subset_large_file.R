#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# version 1.0 (06/23/2018)
# Given the following two inputs:
# (1) infile = input matrix, containing column names (a header) and/or rownames
# (2) subsetlist = a file (or string - see last sentence) containing a list of values that 
# 	match the column names (if subsetlist provided with --columns) or row names (if subsetlist 
#	provided with --rows) or both (can provide both --columns and --rows if needed) in infile. 
# 	subsetlist can also be a comma-separated list of column or row names, e.g. --cols A,B (just
#	make sure not to include any spaces before or after the commas). Note that if the value
#	provided to --cols or --rows contains a period (.), the script will assume it is a filename.
#	Don't put periods in your column or row names. Just don't do it.
# this script will output the subset of columns in infile whose names match the values
# in subsetlist.

# For example:

# infile.txt
#	A	B	C
#	0	0	1
#	1	1	1
#	1	0	0

# subsetlist.txt (provided with --cols)
#	A
#	C

# running subset_large_file.R infile.txt outfile.txt --cols subsetlist.txt 
# will result in the following output file:
#	A	C
#	0	1
#	1	1
#	1	0

# running subset_large_file.R infile.txt outfile.txt --cols A,C
# will produce the same result.

# WARNING: there is NO check that the column or row names you are providing are unique
# (e.g. occur only once in the file). So running the following will get you:

# subset_large_file.R infile.txt outfile.txt --rows A,1

# outfile.txt:
#	A	B	C
#	1	1	1
#	1	0	0

# Even though you only provided two values to --rows, one of those occurs in two rows,
# so the output has three rows.

# NOTE: Original sort order is preserved. So if you run:
# subset_large_file.R infile.txt outfile.txt --rows 1,A

# outfile.txt:
#	A	B	C
#	1	1	1
#	1	0	0

# Only the first column of a file given to --cols or --rows will be used. If no values
# match the column/row names in infile, then your output will be blank (but no error
# will be reported).

# Usage: 	subset_large_file.R infile outfile --cols columnsubsetfile.txt --rows rowsubsetfile.txt
# or		subset_large_file.R infile outfile --cols columnslist --rows rowslist
# (either --rows or --cols can be omitted, but at least one must be present)

# -------------------------
# Version history:
# v.1.0: initial build

options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)	
if (length(args) <= 2) {
	cat("Usage: subset_by_column.R infile outfile --cols colslist --rows rowslist\n")
	cat("----------------------\n")
	cat("infile is a file with row names (if using --rows) and/or column names (if using --cols)\n")
	cat("outfile is the output file name\n")
	cat("Must provide one of either --rows or --cols, or can provide both\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--rows : either a filename (containing a list of row names from infile) or comma-separated list (e.g. --rows A,B,D)\n")
	cat("--cols : either a filename (containing a list of column names from infile) or comma-separated list (e.g. --cols A,B,D)\n")
	cat("--sortrows : if file provided with --rows, output will also be sorted according to order in --rows file (ignored if no --rows provided)\n")
	cat("--sortcols : if file provided with --cols, output will also be sorted according to order in --cols file (ignored if no --cols provided)\n")
	cat("----------------------\n")
	stop("At least one required argument not provided. See usage above.")
}

option_list <- list(
	make_option("--cols", help = "either a filename (containing a list of row names from infile) or comma-separated list (e.g. --rows A,B,D)"),
	make_option("--rows", help = "either a filename (containing a list of column names from infile) or comma-separated list (e.g. --cols A,B,D)"),
	make_option("--sortrows", default=FALSE, action="store_true", help = "Sort output according to --rows file"),
	make_option("--sortcols", default=FALSE, action="store_true", help = "Sort output according to --cols file (can be used with --sortrows)")
	)

arguments = parse_args(OptionParser(option_list=option_list), positional_arguments = 2)
opt = arguments$options
infile = arguments$args[1]
outfile = arguments$args[2]

cat("----------------------\n")
cat("Subsetting the following file:",infile,"\n")
cat("Saving result to:",outfile,"\n")
if (length(opt$cols) != 0) {
	if (grepl("\\.",opt$cols)) {
		cat("Outputting all columns matching values provided in file",opt$cols,"\n")
		colslist = read.table(opt$cols, sep="\t", stringsAsFactors = F)$V1
		cat(length(colslist),"column IDs detected\n")
	} else {
		cat("Outputting all columns matching values provided:",opt$cols,"\n")
		colslist = strsplit(opt$cols,',')[[1]]
		cat(length(colslist),"column IDs detected\n")
	}
	if (opt$sortcols == FALSE) {
		cat("Order of columns in",infile,"will be preserved\n")
	} else {
		cat("Columns in",infile,"will be sorted to match",opt$cols,"\n")
	}
}		
if (length(opt$rows) != 0) {
	if (grepl("\\.",opt$rows)) {
		cat("Outputting all rows matching values provided in file",opt$rows,"\n")		
		rowslist = read.table(opt$rows, sep="\t", stringsAsFactors = F)$V1
		cat(length(rowslist),"row IDs detected\n")
	} else {
		cat("Outputting all rows matching values provided:",opt$rows,"\n")
		rowslist = strsplit(opt$rows,',')[[1]]
		cat(length(rowslist),"row IDs detected\n")
	}
	if (opt$sortrows == FALSE) {
		cat("Order of rows in",infile,"will be preserved\n")
	} else {
		cat("Rows in",infile,"will be sorted to match",opt$rows,"\n")
	}
}
cat("----------------------\n")


mydata = read.table(infile, sep="\t", stringsAsFactors = F)
cat("Input file has",nrow(mydata),"rows and",ncol(mydata),"columns\n")

if (length(opt$rows) != 0 && length(opt$cols) != 0) {
	colslist = c(mydata[1,1],colslist)
	rowslist = c(mydata[1,1],rowslist)
	mydatasubset = mydata[mydata[,1] %in% rowslist, mydata[1,] %in% colslist]
} else if (length(opt$cols) != 0) {
	colslist = c(mydata[1,1],colslist)
	mydatasubset = mydata[,mydata[1,] %in% colslist]
} else {
	rowslist = c(mydata[1,1],rowslist)
	mydatasubset = mydata[mydata[,1] %in% rowslist,]
}	

# sort if requested
if (opt$sortrows == TRUE) {
	mydatasubset = mydatasubset[match(rowslist, mydatasubset[,1]),]
}

if (opt$sortcols == TRUE) {
	mydatasubset = t(mydatasubset)
	mydatasubset = mydatasubset[match(colslist, mydatasubset[,1]),]
	mydatasubset = t(mydatasubset)	
}

cat("Subset has",nrow(mydatasubset),"rows and",ncol(mydatasubset),"columns\n")
write.table(mydatasubset, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na="")






















