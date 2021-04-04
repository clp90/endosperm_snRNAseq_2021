#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(DEsingle))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gplots))

# version 1.0 (03/04/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 03/04/2019 (CLP)
# v.1.1: minor update - -7/09/2019
#	- added filtering options based on total detected counts
# -------------------------

# Description:
# This is a wrapper for DEsingle (Zhun Miao, Ke Deng, Xiaowo Wang, Xuegong Zhang (2018). DEsingle for detecting three types of differential expression in single-cell RNA-seq data. Bioinformatics, bty332.)
# Also see:
# https://bioconductor.org/packages/release/bioc/vignettes/DEsingle/inst/doc/DEsingle.html

# Note: the counts matrix will be subset to only group1, group2 cells, but can contain cells
# that aren't in either group1 or group2 (those will be ignored).

# Usage: run_DEsingle.R [options] counts.txt group1 group2 outprefix

parser = ArgumentParser()
parser$add_argument("counts", nargs=1, help = "matrix of raw expression counts (cols:cells x rows:genes)")
parser$add_argument("group1", nargs=1, help = "list of cells (names must correspond to column names in 'counts') in first condition to be compared; NO HEADER")
parser$add_argument("group2", nargs=1, help = "list of cells (names must correspond to column names in 'counts') in second condition to be compared; NO HEADER")
parser$add_argument("outprefix", nargs=1, help = "prefix for output files")
parser$add_argument("--minreads", default = 5, help = "Genes with fewer than this total value across all cells, are censored. Note that this filter is applied before any transformations (TPM, log2) requested.", type="integer")
parser$add_argument("--mincellsfrac", default = 0.1, help = "Genes detected in fewer than fraction of cells in either group 1 or group 2 are censored.", type="double")
parser$add_argument("--pval", default = 0.05, help = "FDR-adjusted pvalue cutoff for significance", type="double")
parser$add_argument("--pseudocount", default = 0.01, help = "Pseudocount to add before taking log2(TPM)", type="double")

opt <- parser$parse_args()

countsfile = opt$counts
group1file = opt$group1
group2file = opt$group2
outprefix = opt$outprefix
minreads = opt$minreads
mincellsfrac = opt$mincellsfrac
pval = opt$pval
pseudocount = opt$pseudocount

args <- commandArgs(trailingOnly = TRUE)	
if (length(opt) <= 1) {
	cat("Usage: run_DEsingle.R [options] counts.txt group1 group2 outprefix\n")
	cat("----------------------\n")
	cat("TODO\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("TODO\n")
	cat("----------------------\n")
}

cat("\nRunning run_DEsingle.R v.1.0 (03/04/2019)\n")
cat("-----------------------\n")
cat("Current directory:",getwd(),"\n")
cat("Counts matrix:",countsfile,"\n")
cat("List of cells in group 1:",group1file,"\n")
cat("List of cells in group 2:",group2file,"\n")
cat("Prefix for output files:",outprefix,"\n")
cat("-----------------------\n")

counts = data.matrix(read.table(countsfile, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
group1 = data.matrix(read.table(group1file, header=FALSE, row.names=1, sep="\t", stringsAsFactors = FALSE))
group2 = data.matrix(read.table(group2file, header=FALSE, row.names=1, sep="\t", stringsAsFactors = FALSE))

# get subset of count matrix corresponding to the two groups
grpsidx = match(c(rownames(group1),rownames(group2)), colnames(counts))
mydata = counts[,grpsidx]

# set up groups vector
group1size = length(rownames(group1))
group2size = length(rownames(group2))
group = factor(c(rep(1,group1size), rep(2,group2size)))
cat(group1size,"cells detected in group1, and",group2size,"in group2\n")
cat("Original matrix contains",dim(mydata)[1],"genes before filtering\n")

# filter matrix
mydata = mydata[rowSums(mydata)>=minreads, ]
if (ncol(mydata) == 0) {
	stop("Oops, no genes remaining after filtering out all genes with fewer than ",minreads," reads detected across all cells in group1 and group2")
}
tmp1 = mydata[,colnames(mydata) %in% rownames(group1)]
tmp2 = mydata[,colnames(mydata) %in% rownames(group2)]

mydata = mydata[((rowSums(tmp1>0) / length(colnames(tmp1)) > mincellsfrac) | (rowSums(tmp2>0) / length(colnames(tmp2)) > mincellsfrac)), ]
if (ncol(mydata) == 0) {
	stop("Oops, no genes remaining after filtering out all genes expressed in fewer than ",mincells," cells in group1 and group2")
}
totgenes_filt = nrow(mydata)
cat("After filtering out poorly expressed genes,",totgenes_filt,"genes remain\n")

# run DEsingle
results <- DEsingle(counts = mydata, group = group)

results.classified <- DEtype(results = results, threshold = pval)
results.sig <- results.classified[results.classified$pvalue.adj.FDR < pval, ]

# summarize results
degup = rownames(results.sig[results.sig$Type == "DEg" & results.sig$State == "up",])
degdwn = rownames(results.sig[results.sig$Type == "DEg" & results.sig$State == "down",])
deaup = rownames(results.sig[results.sig$Type == "DEa" & results.sig$State == "up",])
deadwn = rownames(results.sig[results.sig$Type == "DEa" & results.sig$State == "down",])
desup = rownames(results.sig[results.sig$Type == "DEs" & results.sig$State == "up",])
desdwn = rownames(results.sig[results.sig$Type == "DEs" & results.sig$State == "down",])

cat("\nSummary of DE analysis results:\n")
cat("Total number of genes assayed:",totgenes_filt,"\n")
cat("Total number of general DE genes (upreg in 1 vs. 2):",length(degup),"\n")
cat("Total number of general DE genes (downreg in 1 vs. 2):",length(degdwn),"\n")
cat("Total number of DE status genes (upreg in 1 vs. 2):",length(desup),"\n")
cat("Total number of DE status genes (downreg in 1 vs. 2):",length(desdwn),"\n")
cat("Total number of DE abundance genes (upreg in 1 vs. 2):",length(deaup),"\n")
cat("Total number of DE abundance genes (downreg in 1 vs. 2):",length(deadwn),"\n")
totde = length(degup)+length(degdwn)+length(desup)+length(desdwn)+length(deaup)+length(deadwn)
cat("Total number of DE genes:",totde,"\n")


# output results
write.table(results.classified,file=paste(outprefix,"_all.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')

if (totde > 0) {
	write.table(results.sig,file=paste(outprefix,"_sig_DE.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')

	if (totde >= 2) {
		# make heatmap of DE gene expression (if any DE genes detected)
		delist = c(degup,degdwn,deadwn,deaup,desup,desdwn)
		x = mydata[rownames(mydata) %in% delist,]
		a = results.sig[order(results.sig$Type,results.sig$State),]
		combos <- with(a, paste(Type, State, sep="."))

		totreads = colSums(x)
		x = t(t(x) / totreads*1000000)
		x = log2(x + pseudocount)
		x = x[rownames(a),]

		hmcolorscale = viridis(30)
		col1 <- brewer.pal(12, "Set3")

		png(paste(outprefix,"_heatmap.png",sep=''), width = 8.5, height = 8, units = 'in', res = 300)
		heatmap.2(x, Colv = FALSE, Rowv=FALSE, dendrogram="none", trace="none", density.info="none", col = hmcolorscale, ColSideColors=col1[group], RowSideColors=col1[as.factor(combos)])
		graphics.off()
	}
}


