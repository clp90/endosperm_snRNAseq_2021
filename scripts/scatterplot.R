#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
library(RColorBrewer)

# version 1.0 (11/01/2015)
# -------------------------
# Version history:
# v.1.0: initial build - 11/01/2015
# v.1.1: 12/06/2015
# 	- removed --xfactor option (generally unused, weird behavior)
# 	- fixed issue where min/max was not calculated properly if --lims not set (needed to set na.rm=TRUE)
# v.1.2: 08/17/2016
# 	- added option --PDF to output as PDF instead of PNG
# -------------------------

# Description:
# given a tab-separated file containing at least two numeric columns, makes a scatterplot of
# the data. User can specify which columns to plot, as well as a "color" column for color-coding
# the points. User can also specify plot title and axes, as well as horizontal and vertical dashed
# lines for reference.

# Usage: scatterplot.R [options] infile.txt outfile.png

parser = ArgumentParser()
parser$add_argument("infile", nargs=1, help = "tab-delimited text file containing 2 numeric columns")
parser$add_argument("outfile", nargs=1, help = "name for output file (will be in PNG format)")
parser$add_argument("--xval", type="integer", default = 1, help = "Column number containing values to plot along x-axis")
parser$add_argument("--yval", type="integer", default = 2, help = "Column number containing values to plot along y-axis")
parser$add_argument("--color", type="integer", help = "Column number containing colors to use for points, default none (all will be black)")
parser$add_argument("--size", type="integer", help = "Column number containing sizes to use for points, default none (default all points 0.25)")
parser$add_argument("--sizenum", type="double", help = "Size to use for points, overruled if --size is also defined (default 0.25)")
parser$add_argument("--shape", default=20, type="integer", help = "Shape for points (see pch symbols)")
parser$add_argument("--xline", help = "Plot vertical line at x = this value (default none)")
parser$add_argument("--yline", help = "Plot horizontal line at y = this value (default none)")
parser$add_argument("--lcolor", default="gray", help = "Color for line(s) [not used for --fitline line] (default \"gray\")")
parser$add_argument("--aline", type="double", help = "Plot line y = a + bx with this value of a")
parser$add_argument("--bline", type="double", help = "Plot line y = a + bx with this value of b")
parser$add_argument("--xlimlo", type="double", help = "Lower limit for x-axis - default none")
parser$add_argument("--xlimhi", type="double", help = "Upper limit for x-axis - default none")
parser$add_argument("--ylimlo", type="double", help = "Lower limit for y-axis - default none")
parser$add_argument("--ylimhi", type="double", help = "Upper limit for y-axis - default none")
parser$add_argument("--xlabel", type="character", help = "Label for x-axis (default name of --xval row)")
parser$add_argument("--ylabel", type="character", help = "Label for y-axis (default name of --yval row)")
parser$add_argument("--xdim", default=6, type="double", help = "Width of plot (in inches)")
parser$add_argument("--ydim", default=6.5, type="double", help = "Height of plot (in inches)")
parser$add_argument("--nbin", type="integer", help = "For smoothscatter - number of bins along each axis to create (e.g. 50 creates 50x50 grid)")
parser$add_argument("--header", default=FALSE, action="store_true", help = "Input file contains a header")
parser$add_argument("--noxticks", default=FALSE, action="store_true", help = "Remove x-axis ticks in plot")
parser$add_argument("--fitline", default=FALSE, action="store_true", help = "Fit a linear regression line to all points")
parser$add_argument("--addeq", default=FALSE, action="store_true", help = "If --fitline is set, add the equation of the regression to the plot")
parser$add_argument("--smooth", default=FALSE, action="store_true", help = "Use 'smoothscatter' function to plot points")
parser$add_argument("--PDF", default=FALSE, action="store_true", help = "Output plot as PDF instead of PNG (default PNG)")
parser$add_argument("--colorramp", type="character", help = "If using --smooth, use these values in color ramp (comma-separated list), if just one element, it should be the name of an RColorBrewer ramp (e.g. YlOrRd)")

opt <- parser$parse_args()

args <- commandArgs(trailingOnly = TRUE)	
if (length(opt) <= 1) {
	cat("Usage: scatterplot.R [options] infile.txt outfile.png\n")
	cat("----------------------\n")
	cat("infile.txt contains at least two numeric columns\n")
	cat("outfile.png is the name for the plot (in png format)\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--xval <int> : column number for values to be plotted along x-axis (default column 1)\n")
	cat("--yval <int> : column number for values to be plotted along y-axis (default column 2)\n")
	cat("--color <int> : column number containing colors of points to be plotted (default none)\n")
	cat("--size <int> : column number containing sizes of points to be plotted (default none)\n")
	cat("--sizenum <int> : size to use for points, overruled if --size is also defined (default 0.25)\n")
	cat("--shape <int> : pch symbol for marker shape (default 20)\n")
	cat("--xline <float> : plot dotted vertical line at x = this value\n")
	cat("--yline <float> : plot dotted vertical line at y = this value\n")
	cat("--aline <float> : plot line y = a+bx with this value of a (must also specify --bline or line won't be plotted)\n")
	cat("--bline <float> : plot line y = a+bx with this value of b (must also specify --aline or line won't be plotted)\n")
	cat("--xlimlo <float> : desired lower limit of x-axis (default none)\n")
	cat("--xlimhi <float> : desired upper limit of x-axis (default none)\n")
	cat("--ylimlo <float> : desired lower limit of y-axis (default none)\n")
	cat("--ylimhi <float> : desired upper limit of y-axis (default none)\n")
	cat("--xdim : width of plot in inches (default 6)\n")
	cat("--ydim : height of plot in inches (default 6.5)\n")
	cat("--header : input file contains a header\n")
	cat("--noxticks : remove tick marks on x-axis in plot\n")
	cat("--fitline : add line of best fit to points\n")
	cat("--addeq : if --fitline specified, also print equation of regression in plot\n")
	cat("--smooth : create smoothed scatterplot\n")
	cat("--colorramp : colorramp for smoothed scatterplot\n")
	cat("--nbin : For smoothscatter - number of bins along each axis to create (e.g. 50 creates 50x50 grid)\n")
	cat("--PDF : output plot as PDF instead of PNG\n")
	cat("----------------------\n")
}

infile = opt$infile
outfile = opt$outfile

alldata = read.table(infile, header=opt$header, sep="\t", stringsAsFactors = FALSE)
alldata[,opt$xval] = as.numeric(alldata[,opt$xval])
alldata[,opt$yval] = as.numeric(alldata[,opt$yval])

if (length(opt$xlabel) == 0) {
	xlabel = names(alldata)[opt$xval]
} else {
	xlabel = opt$xlabel
}

if (length(opt$ylabel) == 0) {
	ylabel = names(alldata)[opt$yval]
} else {
	ylabel = opt$ylabel
}

if (length(opt$color) == 0) {
	alldata$colorforplot = "gray40"
} else {
	alldata$colorforplot = alldata[,opt$color]
}

if (length(opt$sizenum) == 0) {
	alldata$sizeforplot = 0.25
} else {
	alldata$sizeforplot = opt$sizenum
}

if (length(opt$size) != 0) {
	alldata$sizeforplot = alldata[,opt$size]
}

xlim=c(min(alldata[,opt$xval],na.rm=TRUE),max(alldata[,opt$xval],na.rm=TRUE))
if (length(opt$xlimlo) != 0) {
	xlim[1] = as.numeric(opt$xlimlo)
}
if (length(opt$xlimhi) != 0) {
	xlim[2] = as.numeric(opt$xlimhi)
}

ylim=c(min(alldata[,opt$yval],na.rm=TRUE),max(alldata[,opt$yval],na.rm=TRUE))
if (length(opt$ylimlo) != 0) {
	ylim[1] = as.numeric(opt$ylimlo)
}
if (length(opt$ylimhi) != 0) {
	ylim[2] = as.numeric(opt$ylimhi)
}

if(typeof(xlim) == "character" | typeof(ylim) == "character") {
	print(xlim)
	print(ylim)
	stop("Non-numeric characters found in x or y column, does your input file have a header? If so, you must use the --header option.\n")
}


# make scatterplot
if (opt$PDF == TRUE) {
	pdf(outfile, width = opt$xdim, height = opt$ydim, useDingbats=FALSE)
} else {
	png(outfile, width = opt$xdim, height = opt$ydim, units = 'in', res = 300)
}
if (opt$smooth == TRUE) {
	if (length(opt$colorramp) != 0) {
		colorramp = strsplit(opt$colorramp, split=",")[[1]]
		if (length(colorramp) == 1) {
			colorramp = colorRampPalette(brewer.pal(9,colorramp[1]))
		} else {
			colorramp = colorRampPalette(colorramp)
		}
	} else {
		colorramp = colorRampPalette(c("white", blues9))
	}
		
	if (length(opt$nbin) != 0) {

		if (opt$noxticks == TRUE) {
			a = smoothScatter(alldata[,opt$xval], alldata[,opt$yval], nrpoints=0, nbin=opt$nbin,
				xlab = "", sub = paste(xlabel,"\n"), 
				ylab = "", colramp = colorramp,
				pch = opt$shape, cex=alldata$sizeforplot, ylim=ylim, xlim=xlim, col=alldata$colorforplot, xaxt='n')
			a = a + mtext(paste(ylabel,"\n"), 2, line=1)
		} else {
#			a = smoothScatter(alldata[,opt$xval], alldata[,opt$yval], postPlotHook = fudgeit, nrpoints=0, nbin=opt$nbin,
			a = smoothScatter(alldata[,opt$xval], alldata[,opt$yval], nrpoints=0, nbin=opt$nbin,
				xlab = "", sub = paste(xlabel,"\n"), 
				ylab = "", colramp = colorramp,
				pch = opt$shape, cex=alldata$sizeforplot, ylim=ylim, xlim=xlim, col=alldata$colorforplot)
			a = a + mtext(paste(ylabel,"\n"), 2, line=1)
		}
	} else {

		if (opt$noxticks == TRUE) {
			a = smoothScatter(alldata[,opt$xval], alldata[,opt$yval], nrpoints=0,
				xlab = "", sub = paste(xlabel,"\n"), 
				ylab = "", colramp = colorramp,
				pch = opt$shape, cex=alldata$sizeforplot, ylim=ylim, xlim=xlim, col=alldata$colorforplot, xaxt='n')
			a = a + mtext(paste(ylabel,"\n"), 2, line=1)
		} else {
			a = smoothScatter(alldata[,opt$xval], alldata[,opt$yval], nrpoints=0,
				xlab = "", sub = paste(xlabel,"\n"), 
				ylab = "", colramp = colorramp,
				pch = opt$shape, cex=alldata$sizeforplot, ylim=ylim, xlim=xlim, col=alldata$colorforplot)
			a = a + mtext(paste(ylabel,"\n"), 2, line=1)
		}
	}
} else {
	if (opt$noxticks == TRUE) {
		plot(alldata[,opt$xval], alldata[,opt$yval],
			xlab = "", sub = paste(xlabel,"\n"), 
			ylab = "", mtext(paste(ylabel,"\n"), 2, line=1),
			pch = opt$shape, cex=alldata$sizeforplot, ylim=ylim, xlim=xlim, col=alldata$colorforplot, xaxt='n')
	} else {
		plot(alldata[,opt$xval], alldata[,opt$yval],
			xlab = "", sub = paste(xlabel,"\n"), 
			ylab = "", mtext(paste(ylabel,"\n"), 2, line=1),
			pch = opt$shape, cex=alldata$sizeforplot, ylim=ylim, xlim=xlim, col=alldata$colorforplot)
	}
}
if (length(opt$xline) != 0) {
	xlines = strsplit(opt$xline, split=",")[[1]]
	for (i in xlines) {
		abline(v = as.numeric(i), col = opt$lcolor)
	}
}
if (length(opt$yline) != 0) {
	ylines = strsplit(opt$yline, split=",")[[1]]
	for (i in ylines) {
		abline(h = as.numeric(i), col = opt$lcolor)
	}
}
if (length(opt$aline) != 0 && length(opt$bline) != 0) {
	abline(a = opt$aline, b = opt$bline, col = opt$lcolor, lty=2)
}
if (opt$fitline == TRUE) {
	fit = lm(alldata[,opt$yval]~alldata[,opt$xval])
	abline(coef(fit)[1:2], col="black",lty=2)
	if (opt$addeq == TRUE) {
		cf = round(coef(fit), 2) 
		eq = paste0("y = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")
		mtext(eq, 3, line=-2)
	}
}
graphics.off()











