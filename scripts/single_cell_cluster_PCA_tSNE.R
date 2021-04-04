#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scater))
library(fpc)
library(princurve)

# version 1.1 (05/15/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 06/24/2018
# v.1.1: 05/15/2019
#	- added option to cluster genes instead of cells (see --genes option)
# -------------------------

# Description:
# Accepts a matrix of raw counts, converts these into a SingleCellExperiment via sceset,
# calculates CPM (counts per million) using library size, and takes log2 of those values
# (this normalization approach helps minimize bias towards highly expressed transcripts). 

# Then performs basic PCA, and uses k-means clustering over the first 3 PCs to learn
# the optimal k, then clusters using k-means.

# Usage: single_cell_cluster_PCA_tSNE.R [options] expr_matrix.txt cell_info.txt outprefix

set.seed(123456)			# because k-means uses random starts, need to set seed in order to be reproducible

parser = ArgumentParser()
parser$add_argument("expr_matrix", nargs=1, help = "tab-delimited text file matrix of raw expression counts (cols:cells x rows:genes)")
parser$add_argument("cell_info", nargs=1, help = "tab-delimited text file containing one row per cell + at least one descriptor (e.g. age, tissue, etc.)")
parser$add_argument("outprefix", nargs=1, help = "prefix for output files")
parser$add_argument("--genelist", help = "use only genes in genelist for analysis")
parser$add_argument("--k", help = "override k", type="integer")
parser$add_argument("--minreads", default = 10, help = "Genes with fewer than this many total reads detected, across all cells, are censored", type="integer")
parser$add_argument("--mincells", default = 5, help = "Genes detected in fewer than this many cells are censored", type="integer")
parser$add_argument("--minexpr", default = 1000, help = "Nuclei detecting fewer than this many genes are censored", type="integer")
parser$add_argument("--maxK", default = 20, help = "Max allowed value of k (# of clusters) that will be tested to find optimum", type="integer")
parser$add_argument("--ncomponent", default = 5, help = "Number of principal components to be plotted together in initial plot", type="integer")
parser$add_argument("--ntop", default = 500, help = "Use the ntop most variably expressed genes for PCA", type="integer")
parser$add_argument("--perplexity", default = "50", help = "Perplexity parameter for t-SNE; if comma-separated list provided (e.g. 10,20,50) will try all values", type="character")
parser$add_argument("--tSNE_PCs", default = "3", help = "Must be ≤ 10. Use this many of the top PCs from regular PCA as input to tSNE; if comma-separated list provided (e.g. 2,3,5,10) will try all values", type="character")
parser$add_argument("--genes", default=FALSE, action="store_true", help = "Cluster genes instead of cells (default cells)")
parser$add_argument("--princurve", default=FALSE, action="store_true", help = "Use R package 'princurve' to fit a principal curve to the tSNE points")

opt <- parser$parse_args()

if (length(opt$k) != 0) {
	if (opt$k > opt$maxK) {
		opt$maxK = opt$k
	}
}

perplexity = as.numeric(strsplit(opt$perplexity,',')[[1]])
tSNE_PCs = as.numeric(strsplit(opt$tSNE_PCs,',')[[1]])

args <- commandArgs(trailingOnly = TRUE)	
if (length(opt) <= 1) {
	cat("Usage: single_cell_cluster_PCA.R [options] expr_matrix.txt cell_info.txt outprefix\n")
	cat("----------------------\n")
	cat("expr_matrix.txt : tab-delimited text file matrix of raw expression counts (cols:cells x rows:genes)\n")
	cat("cell_info.txt : tab-delimited text file containing one row per cell + at least one descriptor (e.g. age, tissue, etc.)\n")
	cat("outprefix : prefix for output files\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--genelist : use only genes in this list for analysis\n")
	cat("--k : override estimate for k\n")
	cat("----------------------\n")
}

fexpr_matrix = opt$expr_matrix
fcell_info = opt$cell_info
outprefix = opt$outprefix

cat("\nRunning single_cell_cluster_PCA.R v.1.0 (06/24/2018)\n")
cat("-----------------------\n")
cat("Current directory:",getwd(),"\n")
cat("Expression matrix:",fexpr_matrix,"\n")
cat("Cell info table:",fcell_info,"\n")
cat("Prefix for output files:",outprefix,"\n")
cat("-----------------------\n")
cat("Quality filtering settings:",fexpr_matrix,"\n")
cat("Keeping only cells expressing at least this many genes:",opt$minexpr,"\n")
cat("Keeping only genes expressed in at least this many cells:",opt$mincells,"\n")
cat("Keeping only genes with at least this many reads across cells:",opt$minreads,"\n")
cat("-----------------------\n")
cat("\n")

expr_matrix = read.table(fexpr_matrix, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
cell_info = read.table(fcell_info, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
features = colnames(cell_info)

if (length(opt$genelist) != 0) {
	genelist = read.table(opt$genelist, header=FALSE, stringsAsFactors=FALSE)	
	expr_matrix_touse = subset(expr_matrix, rownames(expr_matrix) %in% genelist$V1)
	cat("Subsetting matrix to only use genes listed in --genelist;",nrow(expr_matrix_touse),"rows remain.\n")
	cat("Outputting subsetted matrix to ",outprefix,"_subsetgenelist.txt\n",sep='')
	write.table(expr_matrix_touse, file=paste(outprefix,"_subsetgenelist.txt",sep=''), sep='\t', quote=FALSE)
	expr_matrix = expr_matrix_touse
}

# Filter matrix based on quality
cat("Filtering input matrix for low quality cells and very lowly expressed genes\n")
totgenes = nrow(expr_matrix)
totcells = ncol(expr_matrix)
cat("Original matrix contains",totgenes,"genes and",totcells,"cells\n")
expr_matrix = expr_matrix[, colSums(expr_matrix>0)>=opt$minexpr]
if (ncol(expr_matrix) == 0) {
	stop("Oops, no cells remaining after filtering out cells expressing fewer than ",opt$minexpr," genes")
}
totcells_filt = ncol(expr_matrix)
expr_matrix = expr_matrix[rowSums(expr_matrix)>=opt$minreads, ]
if (ncol(expr_matrix) == 0) {
	stop("Oops, no genes remaining after filtering out all genes with fewer than ",opt$minreads," reads detected across all cells")
}
expr_matrix = expr_matrix[rowSums(expr_matrix>0)>=opt$mincells, ]
if (ncol(expr_matrix) == 0) {
	stop("Oops, no genes remaining after filtering out all genes expressed in fewer than ",opt$mincells," cells")
}
totgenes_filt = nrow(expr_matrix)
cat("QC filtered matrix contains",totgenes_filt,"genes and",totcells_filt,"cells\n")

# If some cells were dropped from analysis (or descfile contains extra cells) drop them from the cell_info data too
tokeep = colnames(expr_matrix)
cell_info = cell_info[tokeep, ]

# Make SCE object, normalize reads to counts per million and log transform
cat("Creating SingleCellExperiment object using scater and converting counts to log2(CPM)\n")
expr_matrix_fx = as.matrix(expr_matrix)
gene_df <- DataFrame(Gene = rownames(expr_matrix_fx))
sceset = SingleCellExperiment(assays = list(counts = expr_matrix_fx), colData = cell_info, rowData = gene_df)
logcounts(sceset) <- log2(calculateCPM(sceset, use.size.factors = FALSE) + 1)

# if analysis was requested over genes instead of cells (--genes) flag, get counts matrix from
# filtered, transformed sceset object and make into pseudo-sceset
if (opt$genes == TRUE) {
	expr_matrix_filt = t(counts(sceset))
	cell_df <- DataFrame(Cell = rownames(expr_matrix_filt))
	sceset = SingleCellExperiment(assays = list(logcounts = expr_matrix_filt), rowData = cell_df, colData = gene_df)	
	cell_info = gene_df
	features = colnames(colData(sceset))
} 

# Perform initial PCA
cat("Performing initial PCA and outputting plot over the first",opt$ncomponent,"PCs (tune this with --ncomponent)\n")
cat("Outputting plot to",paste(outprefix,"_PCA.png",sep=''),"\n")
png(paste(outprefix,"_PCA.png",sep=''), width = 6.5, height = 6, units = 'in', res = 300)
sceset = runPCA(sceset, ntop=opt$ntop, ncomponents = 10)		# keep first 3 PCs
sceset = plotPCA(sceset, return_SCE = TRUE)
graphics.off()
pcmat = (reducedDim(sceset)[,1:3])		# grab first 3 PCs

# Learn optimal number of clusters k (allow override if specified by user)
cat("Learning optimal number of clusters k\n")
k.max <- opt$maxK		# highest number of clusters that will be tested
wss <- sapply(1:k.max, function(k){kmeans(pcmat, k, nstart=50,iter.max = 15 )$tot.withinss})
wssd = diff(wss)
cat("Plot of within-cluster sum of squares for different values of k saved to",paste(outprefix,"_PCA_optimalK.png",sep=''),"\n")
png(paste(outprefix,"_PCA_optimalK.png",sep=''), width = 6.5, height = 6, units = 'in', res = 300)
plot(1:k.max, wss, type="b", pch = 19, frame = FALSE,  xlab="Number of clusters K", ylab="Total within-clusters sum of squares")
graphics.off()
cat("Plot of within-cluster sum of squares for (K-(K-1)) saved to",paste(outprefix,"_PCA_optimalK.png",sep=''),"\n")
png(paste(outprefix,"_PCA_optimalK_diff.png",sep=''), width = 6.5, height = 6, units = 'in', res = 300)
plot(2:k.max, wssd, type="b", pch = 19, frame = FALSE,  xlab="Number of clusters K", ylab="Total within-clusters sum of squares for (K-(K-1))")
graphics.off()

# Use crude elbow method to ID # of clusters
mwss = 5*max(wssd)
wssd = wssd[wssd < mwss]
k_opt = length(wssd)-1
cat("Optimal K:",k_opt,"\n")
if (length(opt$k) != 0) {
	cat("Overriding k with user-supplied value:",opt$k,"\n")
	k_opt = opt$k
}

# Perform k-means clustering and outputting plots
cat("Performing k-means clustering for",k_opt,"clusters and outputting results\n")
kkres = kmeans(pcmat, k_opt, nstart=50,iter.max = 15 )
nnames = colnames(colData(sceset))
colData(sceset) = cbind(colData(sceset),kkres$cluster)		# add cluster info to stats
colnames(colData(sceset)) = c(nnames,"cluster")
write.table(cbind(pcmat,kkres$cluster), file=paste(outprefix,"_PCA_clusters.txt",sep=''), sep='\t', quote=FALSE)
cell_info$cluster = kkres$cluster
write.table(cbind(cell_info, pcmat), file=paste(outprefix,"_cellinfo_plus_clusters.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)

# Also make plots over all factors (including cluster)
for (i in 1:length(features)) { 
	if (features[i] != "Gene") {
		ff=paste(outprefix,"_PCA_",features[i],".png",sep='')
		png(ff, width = 10.5, height = 10, units = 'in', res = 300)
		a = plotReducedDim(sceset, use_dimred = "PCA", colour_by = features[i], ncomponents=3)
		print(a)
		graphics.off()

		ff=paste(outprefix,"_PCA_PC1_PC2_",features[i],".png",sep='')
		png(ff, width = 6.5, height = 6, units = 'in', res = 300)
		a = plotReducedDim(sceset, use_dimred = "PCA", colour_by = features[i])
		print(a)
		graphics.off()
	}
}

# Also add clusters and other factors to the subsetgenelist and output
cat("Outputting subsetted matrix + clusters to ",outprefix,"_subsetgenelist_wclusters.txt\n",sep='')
for (i in 1:length(features)) { 
	expr_matrix = rbind(expr_matrix,colData(sceset)[[features[i]]])
}
write.table(expr_matrix, file=paste(outprefix,"_subset_wfeatures.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)


# Also perform t-SNE; this t-SNE uses the PCs we just computed (user specifies how many of the top PCs to use)
# Do t-SNE on all combinations of perplexity param. and # of PCs provided by user
for (i in 1:length(perplexity)) { 
	for (j in 1:length(tSNE_PCs)) { 
		cat("Performing t-SNE over the first",tSNE_PCs[j],"PCs with perplexity =",perplexity[i],"\n")
		sceset = runTSNE(sceset, perplexity=perplexity[i], rand_seed=12345, use_dimred="PCA", n_dimred=tSNE_PCs[j], check_duplicates = FALSE)
		tsnematrix = slot(reducedDims(sceset),"listData")$TSNE[,1:2]
			
		# learn optimal k in this context
		wss <- sapply(1:k.max, function(k){kmeans(tsnematrix, k, nstart=50,iter.max = 15 )$tot.withinss})
		wssd = diff(wss)

		# Use crude elbow method to ID # of clusters
		mwss = 5*max(wssd)
		wssd = wssd[wssd < mwss]
		k_opt = length(wssd)-1
		cat("Optimal K:",k_opt,"\n")
		
		# Use user-provided number if requested
		if (length(opt$k) != 0) {
			cat("Overriding k with user-supplied value:",opt$k,"\n")
			k_opt = opt$k
		}
		kkres = kmeans(tsnematrix, k_opt, nstart=50,iter.max = 15 )
		colData(sceset)$cluster = kkres$cluster
		features=colnames(colData(sceset))
		
		# if requested, fit line using princurve
		if (opt$princurve == TRUE) {
			lfit = principal_curve(tsnematrix)
			x = as.data.frame(lfit$s)
		}
		
		# make plot with no coloring
		ff=paste(outprefix,"_tSNE_perp",perplexity[i],"_",tSNE_PCs[j],"PCs.png",sep='')
		png(ff, width = 8.5, height = 8, units = 'in', res = 300)
		a = plotTSNE(sceset)
		print(a)
		graphics.off()
		
		for (k in 1:length(features)) { 
			if (features[k] != "Gene") {
				ff=paste(outprefix,"_tSNE_perp",perplexity[i],"_",tSNE_PCs[j],"PCs_",features[k],".png",sep='')
				png(ff, width = 8.5, height = 8, units = 'in', res = 300)
				if (opt$princurve == TRUE) {
					a = plotTSNE(sceset, colour_by = features[k]) + geom_point(data = x, mapping = aes(x = V1, y = V2))
				} else {
					a = plotTSNE(sceset, colour_by = features[k])
				}
				print(a)
				graphics.off()
			}
		}
		cell_info$cluster = kkres$cluster
		write.table(cbind(cell_info, pcmat), file=paste(outprefix,"_tSNE_perp",perplexity[i],"_",tSNE_PCs[j],"PCs_clusters.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)
		write.table(tsnematrix, file=paste(outprefix,"_tSNE_perp",perplexity[i],"_",tSNE_PCs[j],"PCs_coordinates.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)
	}
}



















