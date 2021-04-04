#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(SC3))
library(pheatmap)

# version 1.0 (6/24/2018)
# -------------------------
# Version history:
# v.1.0: initial build - 6/24/2018
# -------------------------

# Description:
# Uses SC3 (https://www.nature.com/articles/nmeth.4236/) to cluster cells. Input is
# an expression matrix of raw counts, and a file containing information about the cells/
# nuclei (metadata).

# Usage: single_cell_cluster_SC3.R [options] expr_matrix.txt cell_info.txt outprefix


parser = ArgumentParser()
parser$add_argument("expr_matrix", nargs=1, help = "tab-delimited text file matrix of cells x genes")
parser$add_argument("cell_info", nargs=1, help = "tab-delimited text file containing one row per cell + at least one descriptor (e.g. age, tissue, etc.)")
parser$add_argument("outprefix", nargs=1, help = "prefix for output files")
parser$add_argument("--k", type="integer", help = "Number of clusters to make (if not provided, will estimate ideal number)")
parser$add_argument("--minreads", default = 10, help = "Genes with fewer than this many total reads detected, across all cells, are censored", type="integer")
parser$add_argument("--mincells", default = 5, help = "Genes detected in fewer than this many cells are censored", type="integer")
parser$add_argument("--minexpr", default = 1000, help = "Nuclei detecting fewer than this many genes are censored", type="integer")
parser$add_argument("--maxK", default = 30, help = "Max allowed value of k (# of clusters) that will be tested to find optimum", type="integer")
parser$add_argument("--seed", default = 123456, help = "Set seed for reproducibility", type="integer")
parser$add_argument("--nonorm", default=FALSE, action="store_true", help = "Input counts are already normalized, don't do CPM normalization.")

opt <- parser$parse_args()

set.seed(opt$seed)			# setting seed within SC3 currently bugged, so set seed for R globally instead to hopefully get this reproducible

args <- commandArgs(trailingOnly = TRUE)	
if (length(opt) <= 1) {
	cat("Usage: single_cell_cluster_SC3.R [options] expr_matrix.txt cell_info.txt outprefix\n")
	cat("----------------------\n")
	cat("expr_matrix.txt : tab-delimited text file matrix of cells x genes\n")
	cat("cell_info.txt : tab-delimited text file containing one row per cell + at least one descriptor (e.g. age, tissue, etc.)\n")
	cat("outprefix : prefix for output files\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--k : call this many clusters (default lets SC3 estimate best # of clusters)\n")
	cat("----------------------\n")
}

fexpr_matrix = opt$expr_matrix
fcell_info = opt$cell_info
outprefix = opt$outprefix

cat("\nRunning single_cell_cluster_SC3.R v.1.0 (6/24/2018)\n")
cat("-----------------------\n")
cat("Current directory:",getwd(),"\n")
cat("Expression matrix:",fexpr_matrix,"\n")
cat("Cell info table:",fcell_info,"\n")
cat("Prefix for output files:",outprefix,"\n")
cat("-----------------------\n")
if (length(opt$k) == 0) {
	cat("Ideal number of clusters will be estimated by SC3\n")
} else {
	kk = opt$k
	cat("Looking for",kk,"clusters\n")
}
cat("Quality filtering settings:",fexpr_matrix,"\n")
cat("Keeping only cells expressing at least this many genes:",opt$minexpr,"\n")
cat("Keeping only genes expressed in at least this many cells:",opt$mincells,"\n")
cat("Keeping only genes with at least this many reads across cells:",opt$minreads,"\n")
cat("Other settings:",outprefix,"\n")
cat("Random seed:",opt$seed,"\n")
cat("Max allowed value of k =",opt$maxK," (lower this if script is slow)\n")
cat("-----------------------\n")
cat("\n")

expr_matrix = read.table(fexpr_matrix, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
cell_info = read.table(fcell_info, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
features = colnames(cell_info)
celllist = colnames(expr_matrix)

# Filter matrix based on quality
cat("Filtering input matrix for low quality cells and very lowly expressed genes\n")
totgenes = nrow(expr_matrix)
totcells = ncol(expr_matrix)
expr_matrix = expr_matrix[, colSums(expr_matrix>0)>=opt$minexpr]
totcells_filt = ncol(expr_matrix)
expr_matrix = expr_matrix[rowSums(expr_matrix)>=opt$minreads, ]
expr_matrix = expr_matrix[rowSums(expr_matrix>0)>=opt$mincells, ]
totgenes_filt = nrow(expr_matrix)
cat("Original matrix contained",totgenes,"genes and",totcells,"cells\n")
cat("QC filtered matrix contains",totgenes_filt,"genes and",totcells_filt,"cells\n")

# If any cells dropped, drop them also from the cell_info table
if (totcells != totgenes_filt) {
	cellsremaining = colnames(expr_matrix)
	cell_info = cell_info[rownames(cell_info) %in% cellsremaining,]
}

# Make SCE object
expr_matrix = as.matrix(expr_matrix)
expr_df = DataFrame(Gene = rownames(expr_matrix))
sceset = SingleCellExperiment(assays = list(counts = expr_matrix), colData = cell_info, rowData = expr_df)
if (opt$nonorm == FALSE) {
	logcounts(sceset) <- log2(calculateCPM(sceset, use.size.factors = FALSE) + 1)
} else {
	logcounts(sceset) <- log2(counts(sceset) + 1)
}
rowData(sceset)$feature_symbol <- rownames(sceset)

# Run SC3
sceset = calculateQCMetrics(sceset)
sc3data <- sc3(sceset, ks = 2:opt$maxK, biology = TRUE, rand_seed=opt$seed)

if (length(opt$k) == 0) {
	sc3data = sc3_estimate_k(sc3data)
	kk = metadata(sc3data)$sc3$k_estimation
	cat("Estimated k =",kk,"\n")
	if (kk == 1) {
		cat("When k==1, various parts of this script won't work, so I'm bumping k up to 2, but keep in mind that the program recommends only 1 cluster here!\n")
		kk=2
	}
}
cat("Done running SC3\n")

cat("Making consensus plot\n")
png(paste(outprefix,"_consensus.png",sep=''), width = 7.5, height = 9, units = 'in', res = 300)
sc3_plot_consensus(sc3data, k = kk, show_pdata = colnames(cell_info))
graphics.off()

cat("Getting consensus...\n")
consensus = metadata(sc3data)$sc3$consensus[[as.character(kk)]]$consensus
cat("Outputting consensus matrix:\n")
rownames(consensus) = celllist
colnames(consensus) = celllist

cat("Getting hc...\n")
hc = metadata(sc3data)$sc3$consensus[[as.character(kk)]]$hc
#cat("Order of cells from left to right in heatmap:")
cellorder = rownames(cell_info)[hc$order]

# consensus matrix is in original order of nuclei provided in expr_data; re-sort into
# order in clustered matrix
stconsensus = consensus[,match(cellorder, colnames(consensus))]
stconsensus = stconsensus[match(cellorder, colnames(consensus)),]

write.table(data.frame("nuc_ID"=rownames(stconsensus),stconsensus), file=paste(outprefix,"_",kk,"_consensus_matrix.txt",sep=''), sep='\t', quote=FALSE, col.names=T, row.names=F)
write.table(cbind(order=cellorder), file=paste(outprefix,"_left_to_right_order.txt",sep=''), sep='\t', quote=FALSE, col.names=T, row.names=F)
write.table(colData(sc3data), file=paste(outprefix,"_all_cluster_info.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)
write.table(cbind(geneID=rownames(colData(sc3data)),cluster=colData(sc3data)[[paste("sc3_",kk,"_clusters",sep='')]]), file=paste(outprefix,"_",kk,"_cluster_info.txt",sep=''), sep='\t', quote=FALSE, col.names=T, row.names=F)

cat("Making marker genes plot\n")
png(paste(outprefix,"_markergenes.png",sep=''), width = 7.5, height = 30, units = 'in', res = 300)
sc3_plot_markers(sc3data, k = kk, show_pdata = colnames(cell_info))
graphics.off()


cat("Getting features...\n")
show_pdata = features
ann <- cell_info[, colnames(cell_info) %in% show_pdata]

for (i in grep("_log2_outlier_score", colnames(ann))) {
	if (class(ann[, i]) == "factor") {
		ann[, i] <- as.numeric(levels(ann[, i]))[ann[, i]]
	}
}

colnames(consensus) = rownames(ann)

cat("Making plots...\n")
png(paste(outprefix,"_consensus_rowlabels.png",sep=''), width = 7.5, height = 9, units = 'in', res = 300)
pheatmap(consensus, cluster_rows = hc, cluster_cols = hc, cutree_rows = kk, cutree_cols = kk, show_rownames = FALSE, fontsize_col = 3, annotation_col = cell_info)
graphics.off()


png(paste(outprefix,"_expression.png",sep=''), width = 7.5, height = 10, units = 'in', res = 300)
sc3_plot_expression(sc3data, k = kk, show_pdata = features)
graphics.off()

png(paste(outprefix,"_markers.png",sep=''), width = 7.5, height = 30, units = 'in', res = 300)
sc3_plot_markers(sc3data, k = kk, show_pdata = features)
graphics.off()

cat("Outputting cell order, cluster info, marker gene info\n")
# output info about clusters + ordering of cells in the heatmaps
# and marker genes for each cluster
marker_gene_info = cbind(geneID=rowData(sc3data)$Gene,mean_counts=rowData(sc3data)$mean_counts,log10_mean_counts=rowData(sc3data)$log10_mean_counts,rank_counts=rowData(sc3data)$rank_counts,n_cells_counts=rowData(sc3data)$n_cells_counts,sc3_gene_filter=rowData(sc3data)$sc3_gene_filter,sc3__markers_clusts=rowData(sc3data)[[paste("sc3_",kk,"_markers_clusts",sep='')]],sc3__markers_padj=rowData(sc3data)[[paste("sc3_",kk,"_markers_padj",sep='')]],sc3__markers_auroc=rowData(sc3data)[[paste("sc3_",kk,"_markers_auroc",sep='')]])
write.table(marker_gene_info, file=paste(outprefix,"_",kk,"_marker_gene_info.txt",sep=''), sep='\t', quote=FALSE, col.names=T, row.names=F)




