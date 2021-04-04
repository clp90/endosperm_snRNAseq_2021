#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))

# version 1.0 (11/09/2015)
# -------------------------
# Version history:
# v.1.0: initial build - 11/09/2015
# -------------------------

# Description:
# given a tab-separated file containing at least one numeric column, makes a bar chart of
# the data. User can specify which column to plot, as well as a "color" column for color-coding
# the bars and a "label" column to label the x-axis.

# Usage: barchart.R [options] infile.txt outfile.pdf

parser = ArgumentParser()
parser$add_argument("infile", nargs=1, help = "tab-delimited text file containing 2 numeric columns")
parser$add_argument("outfile", nargs=1, help = "name for output file (will be in pdf format)")
parser$add_argument("--yvals", type="integer", default = 2, help = "Column number containing values to use to make chart - default 2")
parser$add_argument("--xvals", type="integer", default = 1, help = "Column number containing values to use to make chart - default 1")
parser$add_argument("--factor", type="integer", help = "Column number containing factor over which to make grouped bar chart - default none")
parser$add_argument("--colors", type="character", help = "Comma-separated list of colors (number == to levels of factor, or 1, if no factor provided)")
parser$add_argument("--flevels", type="character", help = "Comma-separated list of levels of factor (only used if --factor provided) - plots these in order on plot")
parser$add_argument("--xlevels", type="character", help = "Comma-separated list of levels of values on x-axis (use this to change the ordering, default alphabetic)")
parser$add_argument("--xtitle", type="character", help = "Title for x-axis")
parser$add_argument("--ytitle", type="character", help = "Title for y-axis")
parser$add_argument("--yupper", type="double", help = "Upper value for y-axis")
parser$add_argument("--ylower", type="double", help = "Lower value for y-axis")
parser$add_argument("--xdim", default=8, type="double", help = "Width of plot (in inches)")
parser$add_argument("--ydim", default=6, type="double", help = "Height of plot (in inches)")
parser$add_argument("--title", type="character", default="", help = "Title for plot")
parser$add_argument("--header", default=FALSE, action="store_true", help = "Input file contains a header")
parser$add_argument("--vertical", default=FALSE, action="store_true", help = "Plot x axis labels vertically")
parser$add_argument("--median", default=FALSE, action="store_true", help = "Plot median of values over xval and factor instead of mean (default mean)")
parser$add_argument("--stack", default=FALSE, action="store_true", help = "Make stacked barchart over levels of factor")
parser$add_argument("--pstack", default=FALSE, action="store_true", help = "Make 100% stacked barchart over levels of factor")
parser$add_argument("--yline", help = "Plot horizontal line at y = this value (default none)")
parser$add_argument("--spacing", type="double", default=0.9, help = "Spacing between bars")

opt <- parser$parse_args()

args <- commandArgs(trailingOnly = TRUE)	
if (length(opt) <= 1) {
	cat("Usage: piechart.R [options] infile.txt outfile.pdf\n")
	cat("----------------------\n")
	cat("infile.txt contains at least two numeric columns\n")
	cat("outfile.pdf is the name for the plot (in pdf format)\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--val <int> : index of column containing the data to plot in pie chart (default column 1)\n")
	cat("--color <int> : index of column containing colors of pie slices (default will use default color scheme)\n")
	cat("--labels <int> : index of column containing labels for pie slices\n")
	cat("--title <char> : title for plot\n")
	cat("--header : input file contains a header\n")
	cat("----------------------\n")
}

infile = opt$infile
outfile = opt$outfile

if (opt$stack == TRUE && opt$pstack == TRUE) stop("must use either --stack or --pstack")

alldata = read.table(infile, header=opt$header, sep="\t", stringsAsFactors = TRUE)
levelsofx = unique(alldata[,opt$xvals])

if (opt$median == TRUE) {
	ff = "median"
} else {
	ff = "mean"
}

if (length(opt$factor) == 0) {
	meandata = aggregate(alldata[,opt$yvals] ~ alldata[,opt$xvals], FUN= ff )
	colnames(meandata) = c(colnames(alldata)[opt$xvals],colnames(alldata)[opt$yvals])
} else {
	meandata = aggregate(alldata[,opt$yvals] ~ alldata[,opt$xvals] + alldata[,opt$factor], data = alldata, FUN= ff )
	colnames(meandata) = c(colnames(alldata)[opt$xvals],colnames(alldata)[opt$factor],colnames(alldata)[opt$yvals])
}

if (length(opt$xtitle) == 0) {
	xtitle = ""
} else {
	xtitle = opt$xtitle
}
if (length(opt$ytitle) == 0) {
	ytitle = ""
} else {
	ytitle = opt$ytitle
}

if (length(opt$yupper) == 0 && length(opt$ylower) != 0) {
	cat("Warning: value of --ylower will be ignored (to use, must specify both --yupper and --ylower)\n")
} else if (length(opt$yupper) != 0 && length(opt$ylower) == 0) {
	cat("Warning: value of --yupper will be ignored (to use, must specify both --yupper and --ylower)\n")
}

if (length(opt$colors) != 0) {
	colors = strsplit(opt$colors,",")[[1]]
	if (length(opt$factor) == 0) {
		if (length(colors) > 1) {
			stop("If no --factor is specified, can only specify one color")
		}
	} else {
		factorlevels = levels(meandata[,colnames(alldata)[opt$factor]])
		if (length(colors) != length(factorlevels)) {
			stop("If --factor is specified, must provide as many colors with --color as there are levels of factor (there are ",length(factorlevels)," levels of factor)")
		}
	}
}

if (length(opt$flevels) != 0) {
	lvls = strsplit(opt$flevels,",")[[1]]
	if (length(opt$factor) != 0) {
		factorlevels = levels(meandata[,colnames(alldata)[opt$factor]])
		if (length(lvls) != length(factorlevels)) {
			stop("If --factor is specified, must provide as many values with --levels as there are levels of factor (there are ",length(factorlevels)," levels of factor)")
		} else {
			for (i in lvls) {
				if (! i %in% factorlevels) {
					stop("--level ",i," not found in factor vector")
				}
			}		
		}
	}
	meandata[,colnames(alldata)[opt$factor]] <- factor(meandata[,colnames(alldata)[opt$factor]], levels = lvls)
}

if (length(opt$xlevels) != 0) {
	xlvls = strsplit(opt$xlevels,",")[[1]]
	if (length(xlvls) != length(levelsofx)) {
		stop("--xlevels must specify as many values as there are levels of column --xvals in data (there are ",length(levelsofx)," levels in provided --xvals)")
	} else {
		for (i in xlvls) {
			if (! i %in% levelsofx) {
				stop("--xlevel ",i," not found in --xvals column")
			}
		}
	}
	meandata[,colnames(alldata)[opt$xvals]] <- factor(meandata[,colnames(alldata)[opt$xvals]], levels = xlvls)
}

# make the plot
pdf(outfile, width = opt$xdim, height = opt$ydim)
if (length(opt$factor) != 0) {
	a = ggplot(data=meandata, aes_string(x=colnames(alldata)[opt$xvals], y=colnames(alldata)[opt$yvals], fill=colnames(alldata)[opt$factor], width = opt$spacing))
} else {
	a = ggplot(data=meandata, aes_string(x=colnames(alldata)[opt$xvals], y=colnames(alldata)[opt$yvals], width = opt$spacing))
}
if (opt$stack == TRUE) {
	a = a + geom_bar(stat="identity", color="black", position = "stack")
} else if (opt$pstack == TRUE) {
	a = a + geom_bar(stat="identity", color="black", position = "fill")
} else {
	a = a + geom_bar(stat="identity", color="black", position = "dodge")
}
a = a + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab(xtitle) + ylab(ytitle)
if (length(opt$colors) != 0) {
	a = a + scale_fill_manual(values=colors)
}
if (length(opt$yupper) != 0 && length(opt$ylower) != 0) {
	a = a + ylim(opt$ylower, opt$yupper)
}
if (opt$vertical == TRUE) {
	a = a + theme(axis.text.x = element_text(angle = 90, hjust = 1))		# flip x-axis labels vertically
}
if (length(opt$yline) != 0) {
	ylines = strsplit(opt$yline,",")[[1]]
#COM	a = a + geom_hline(aes(yintercept=as.numeric(ylines[1])), color="gray", linetype="dashed") + geom_hline(aes(yintercept=as.numeric(ylines[2])), color="gray", linetype="dashed") + geom_hline(aes(yintercept=as.numeric(ylines[3])), color="gray", linetype="dashed")
	
	for (i in 1:(length(ylines))) {
		cat("Adding line y=",ylines[i],"\n")
		a = a + geom_hline(aes(yintercept=as.numeric(ylines[i])), color="gray", linetype="dashed")
	}
}
print(a)
graphics.off()







