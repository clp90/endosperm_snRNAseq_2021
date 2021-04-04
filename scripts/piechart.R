#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))

# version 1.0 (11/09/2015)
# -------------------------
# Version history:
# v.1.0: initial build - 11/09/2015
# v.1.1: (08/17/2016) - added option to output as PDF instead of PNG
# -------------------------

# Description:
# given a tab-separated file containing at least one numeric column, makes a pie chart of
# the data. User can specify which column to plot, as well as a "color" column for color-coding
# the points and a "label" column to label slices.

# Usage: piechart.R [options] infile.txt outfile.png

parser = ArgumentParser()
parser$add_argument("infile", nargs=1, help = "tab-delimited text file containing 2 numeric columns")
parser$add_argument("outfile", nargs=1, help = "name for output file (will be in PNG format)")
parser$add_argument("--val", type="integer", default = 1, help = "Column number containing values to use to make chart")
parser$add_argument("--colors", type="integer", help = "Column number containing colors to use for points, default none (all will be black)")
parser$add_argument("--labels", type="integer", help = "Column containing labels for the pie slices")
parser$add_argument("--title", type="character", default="", help = "Title for plot")
parser$add_argument("--header", default=FALSE, action="store_true", help = "Input file contains a header")
parser$add_argument("--PDF", default=FALSE, action="store_true", help = "Output plot as PDF instead of PNG (default PNG)")

opt <- parser$parse_args()

args <- commandArgs(trailingOnly = TRUE)	
if (length(opt) <= 1) {
	cat("Usage: piechart.R [options] infile.txt outfile.png\n")
	cat("----------------------\n")
	cat("infile.txt contains at least two numeric columns\n")
	cat("outfile.png is the name for the plot (in png format)\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--val <int> : index of column containing the data to plot in pie chart (default column 1)\n")
	cat("--color <int> : index of column containing colors of pie slices (default will use default color scheme)\n")
	cat("--labels <int> : index of column containing labels for pie slices\n")
	cat("--title <char> : title for plot\n")
	cat("--header : input file contains a header\n")
	cat("--PDF : output plot as PDF instead of PNG\n")
	cat("----------------------\n")
}

# modified pie chart function that won't draw labels if not provided
mypie <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
									init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
									col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
	if (!is.numeric(x) || any(is.na(x) | x < 0)) 
		stop("'x' values must be positive.")
	if (is.null(labels)) 
		labels <- as.character(seq_along(x))
	else labels <- as.graphicsAnnot(labels)
	x <- c(0, cumsum(x)/sum(x))
	dx <- diff(x)
	nx <- length(dx)
	plot.new()
	pin <- par("pin")
	xlim <- ylim <- c(-1, 1)
	if (pin[1L] > pin[2L]) 
		xlim <- (pin[1L]/pin[2L]) * xlim
	else ylim <- (pin[2L]/pin[1L]) * ylim
	dev.hold()
	on.exit(dev.flush())
	plot.window(xlim, ylim, "", asp = 1)
	if (is.null(col)) 
		col <- if (is.null(density)) 
			c("white", "lightblue", "mistyrose", "lightcyan", 
				"lavender", "cornsilk")
	else par("fg")
	if (!is.null(col)) 
		col <- rep_len(col, nx)
	if (!is.null(border)) 
		border <- rep_len(border, nx)
	if (!is.null(lty)) 
		lty <- rep_len(lty, nx)
	angle <- rep(angle, nx)
	if (!is.null(density)) 
		density <- rep_len(density, nx)
	twopi <- if (clockwise) 
		-2 * pi
	else 2 * pi
	t2xy <- function(t) {
		t2p <- twopi * t + init.angle * pi/180
		list(x = radius * cos(t2p), y = radius * sin(t2p))
	}
	for (i in 1L:nx) {
		n <- max(2, floor(edges * dx[i]))
		P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
		polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
						border = border[i], col = col[i], lty = lty[i])
		P <- t2xy(mean(x[i + 0:1]))
		lab <- as.character(labels[i])
	}
	title(main = main, ...)
	invisible(NULL)
}

infile = opt$infile
outfile = opt$outfile

alldata = read.table(infile, header=opt$header, sep="\t", stringsAsFactors = FALSE)

if (length(opt$colors) == 0) {
	colors = c()
} else {
	colors = alldata[,opt$colors]
}

if(typeof(alldata[,opt$val]) == "character") {
	print(alldata[1:10,opt$val])
	stop("Non-numeric characters found in column ",opt$val," (see above); if the input file has a header, you must specify --header\n")
}

if (opt$PDF == TRUE) {
	pdf(outfile, useDingbats=FALSE)
} else {
	png(outfile, width = 6, height = 6.5, units = 'in', res = 300)
}
if (length(opt$labels) == 0) {
	par(mgp=c(1,0.45,0), tcl=-0.4, mar=c(0,0,0,0))
	mypie(alldata[,opt$val], col=colors, clockwise=TRUE, init.angle=90, main=opt$title, radius = 1)
} else {
	lbls = alldata[,opt$labels]
	par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.1,1.1,1.1,1.1))
	pie(alldata[,opt$val], labels = lbls, col=colors, clockwise=TRUE, init.angle=90, main=opt$title)
}
graphics.off()













