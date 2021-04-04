#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))

# version 1.0 (07/12/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 07/12/2019
# -------------------------

# Description:
# -------------------------
# Simple script to create a dot-plot representing overlaps between one set of gene lists and another.
# Note - originally made for comparing lists of genes, but could theoretically compare any set of strings.

# Required inputs:
# -------------------------
# listA : comma-separated list of files from the first set of conditions, each containing a list of genes
# listB : comma-separated list of files from the second set of conditions, each containing a list of genes
# outprefix : prefix for output files

# Notes:
# -------------------------
# - Any duplicate gene IDs in any input file are automatically ignored

# Usage: gene_overlaps_dotplot.R [options] listA listB outprefix

parser = ArgumentParser()
parser$add_argument("listA", nargs=1, help = "comma-separated list of files from the first set of conditions, each containing a list of genes")
parser$add_argument("listB", nargs=1, help = "comma-separated list of files from the second set of conditions, each containing a list of genes")
parser$add_argument("outprefix", nargs=1, help = "prefix for output files")
parser$add_argument("--namesA", default = "", help = "List of names corresponding to the 'listA' files (for labeling x-axis of plot)")
parser$add_argument("--namesB", default = "", help = "List of names corresponding to the 'listB' files (for labeling y-axis of plot)")
parser$add_argument("--scoretype", default = "max", help = "How to report calculated overlaps - 'num' = raw number of shared genes, 'A' = frac of A genes, 'B' = frac of B genes, 'min' = min of 'A' and 'B', 'max' = max of 'A' and 'B', 'hgeo' = -log(p-value from hypergeometric test)")
parser$add_argument("--popsize", help = "Required if --scoretype 'hgeo' is used: total size of population of genes that the lists A and B were drawn from",type="integer")
parser$add_argument("--sizeupper", help = "Upper limit for size scale (default max in data)",type="double")
parser$add_argument("--sizelower", help = "Lower limit for size scale (default max in data)",type="double")
parser$add_argument("--dotsize", default = 50, help = "(optional, only used if plottype = 'dot') - upper limit for size of points in dot plot (default 10)", type="double")
parser$add_argument("--sizetrans", default="identity", help = "Transformation for size scale (see scale_size_continuous() options in ggplot2; available values: 'asn', 'atanh', 'boxcox', 'exp', 'identity', 'log', 'log10', 'log1p', 'log2', 'logit', 'probability', 'probit', 'reciprocal', 'reverse' and 'sqrt')",type="character")
parser$add_argument("--width", default = 8.5, help = "Width of plot (in inches)", type="double")
parser$add_argument("--height", default = 8, help = "Height of plot (in inches)", type="double")
parser$add_argument("--header", default=FALSE, action="store_true", help = "Input counts matrix is already normalized; do not normalize again")

opt <- parser$parse_args()

outprefix = opt$outprefix

allowed_scoretype = c( "num", "A", "B", "min", "max","hgeo" )
scoretype = opt$scoretype
if (! scoretype %in% allowed_scoretype) {
	stop('--scoretype must be one of either "num", "A", "B", "min", "max", or "hgeo"')
}

if (scoretype == "hgeo" && length(opt$popsize) == 0) {
	stop('If using --scoretype "hgeo" you must also specify --popsize')
}

listA = unlist(strsplit(opt$listA,","))
listB = unlist(strsplit(opt$listB,","))

if (opt$namesA != "") {
	namesA = unlist(strsplit(opt$namesA,","))
} else {
	namesA = paste("A",1:length(listA),sep='')
}
if (length(namesA) != length(listA)) {
	stop("Error, listA contains ",length(listA)," files but ",length(namesA)," names were provided with --namesA")
}

if (opt$namesB != "") {
	namesB = unlist(strsplit(opt$namesB,","))
} else {
	namesB = paste("B",1:length(listB),sep='')
}
if (length(namesB) != length(listB)) {
	stop("Error, listB contains ",length(listB)," files but ",length(namesB)," names were provided with --namesB")
}

res <- c()
for (a in 1:length(listA)) {
	for (b in 1:length(listB)) {	
		
		setA = read.table(listA[a], header=opt$header, sep="\t", stringsAsFactors = FALSE)
		setB = read.table(listB[b], header=opt$header, sep="\t", stringsAsFactors = FALSE)
		
		if (dim(setA)[2] != 1) {
			stop("Error: input file ",listA[a]," has more than one column; only one column of values (gene IDs) allowed")
		}
		
		if (dim(setB)[2] != 1) {
			stop("Error: input file ",listB[b]," has more than one column; only one column of values (gene IDs) allowed")
		}
		
		setA = unique(as.list(setA))[[1]]
		setB = unique(as.list(setB))[[1]]
						
		if (scoretype == "num") {
			val = length(intersect(setA,setB))
		} else if (scoretype == "A") {
			val = length(intersect(setA,setB)) / length(setA)
		} else if (scoretype == "B") {
			val = length(intersect(setA,setB)) / length(setB)
		} else if (scoretype == "min") {
			val = min(length(intersect(setA,setB)) / length(setA), length(intersect(setA,setB)) / length(setB))
		} else if (scoretype == "max") {
			val = max(length(intersect(setA,setB)) / length(setA), length(intersect(setA,setB)) / length(setB))
		} else if (scoretype == "hgeo") {
			int = length(intersect(setA,setB))
			if (int == 0) {
				val = 0
			} else {
				minv = min(length(setA),length(setB))
				maxv = max(length(setA),length(setB))
				val = -1*log10(sum(dhyper(int:minv,minv,opt$popsize-minv,maxv)))
			}
		}
		
		obs = c(namesA[a],namesB[b],a,b,val)
		res = rbind(res,obs)	
	}
}

rownames(res) = NULL
res = as.data.frame(res,stringsAsFactors=FALSE)
colnames(res) = c("nameA","nameB","a","b","value")
res$value = as.numeric(res$value)
res$a = as.numeric(res$a)
res$b = as.numeric(res$b)

if (length(opt$sizeupper) != 0) {
	sizeupper = opt$sizeupper
} else {
	sizeupper = max(res$value[is.finite(res$value)])
}

if (length(opt$sizelower) != 0) {
	sizelower = opt$sizelower
} else {
	sizelower = min(res$value)
}

xbreaks = c(unique(res$a))
ybreaks = c(unique(res$b))

if (scoretype == "num") {
	lname = "Number"
} else if (scoretype == "A") {
	lname = "Frac. of A (row)"
} else if (scoretype == "B") {
	lname = "Frac. of B (col)"
} else if (scoretype == "min") {
	lname = "Min. frac overlap"
} else if (scoretype == "max") {
	lname = "Max frac overlap"
} else if (scoretype == "hgeo") {
	lname = "-log10(pval)"
}

# make plot
theme_set(theme_bw())
res$btouse = rev(res$b)

# replace any Inf values in res$value (which can't be plotted) with sizeupper
# and cap values at sizeupper
res$value_to_use = res$value
res$value_to_use[! is.finite(res$value_to_use)] = sizeupper
res$value_to_use[res$value_to_use > sizeupper] = sizeupper

pdf(paste(outprefix,'.pdf',sep=''), width = opt$width, height = opt$height)
a = ggplot(res, aes(x=a,y=btouse, size=value_to_use)) + geom_point(alpha=0.2)
a = a + scale_size_continuous(range = c(0,opt$dotsize), trans = opt$sizetrans, limits=c(sizelower,sizeupper), name=lname)
a = a + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
a = a + scale_x_continuous(breaks=xbreaks, labels=namesA) + scale_y_continuous(breaks=ybreaks, labels=rev(namesB))
a = a + coord_cartesian(xlim=c(0.5, max(res$a) + 0.5),ylim=c(0.5, max(res$b) + 0.5))
a = a + labs(y = "B lists", x= "A lists")
print(a)
graphics.off()

write.table(res,file=paste(outprefix,"_overlaps.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')







