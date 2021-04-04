#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(dplyr))

# version 1.0 (07/11/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 07/11/2019
# v.1.1: 12/18/2019
#	- added additional optional step at end of the analysis, that 
#	- to hierarchical clustering portion, added code to obtain discrete clusters from the clustering
# -------------------------

# Description:
# -------------------------
# Simple script to cluster genes based on expression across different conditions, 
# given a matrix of rows:genes x cols:conditions. Clustering method is regular k-means
# clustering; script can optimize k or you can provide your own value. Can also work
# to assess variability in mat/pat bias across different conditions/factors/clusters.

# Required inputs (normal mode):
# -------------------------
# expr_matrix : gene expression matrix; rows:genes x cols:conditions. Values can be anything (FPKM, etc.).
# outprefix : prefix for output files

# Required inputs (imprinting mode):
# -------------------------
# mcounts : matrix of maternal expression; rows:genes x cols:conditions. Values should be normalized count data ONLY.
# pcounts : matrix of paternal expression; rows:genes x cols:conditions. Values should be normalized count data ONLY.
# outprefix : prefix for output files
# method : either 'combined' or 'separate'

# Notes:
# -------------------------
# (1) How --sampfile works:
#	- can have 1 or 2 columns
#	- first column is expected to contain column names ('samples') from expr_matrix, and only those columns will be kept for analysis
#	- second column is optional, but if provided, is expected to assign each 'sample' to a cluster, which can be any string
#		- if this second column is provided, expression values for all samples in a cluster will be averaged together to get
#		a per-cluster avg.(expression) value, that will be used for clustering instead of the individual samples
#		- additionally, script will randomly shuffle the cluster labels and identify genes with p < --pval significantly enriched
#		or depleted expression in each cluster

# Example:
#	nuc1	clus1
#	nuc2	clus1
#	nuc3	clus1
#	nuc4	clus2
#	nuc5	clus2
#	nuc6	clus2

# (2) How --factors works:
# 	- can have 2 or more columns
# 	- first column must match the samples (column names of expr_matrix)
#	- additional columns will assign those samples to additional factors; these are taken into account when scaling and bootstrapping
#	- will output additional plots comparing between different factors
# 	- can have factors without clusters and vice-versa

# Example:
#	nuc1	endosperm
#	nuc2	seedcoat
#	nuc3	endosperm
#	nuc4	seedcoat

# (3) different runs are only reproducible if --seed is the same (and probably also requires same version of R & libraries)
#
# (4) 'combined' vs. 'separate' methods for imprinting analysis. 'combined' will look for variation by cluster (sampfile values)
# or factor specifically in the %maternal, while 'separate' will model maternal and paternal counts separately, and look separately
# for variation in maternal or paternal expression across clusters instead. It will then calculate a 'difference score', trying to
# find cases where maternal and paternal counts vary across clusters in inconsistent ways.

# Usage: cluster_gene_expression.R [options] expr_matrix.txt outprefix

parser = ArgumentParser()
# arguments required for all runs
parser$add_argument("outprefix", nargs=1, help = "prefix for output files")

# arguments required if looking at variability in gene expression
parser$add_argument("--expr_matrix", nargs=1, help = "gene expression matrix; rows:genes x cols:samples/cells. Values should be normalized (e.g. normalized counts or TPM).")

# arguments required if looking at variability in imprinting
parser$add_argument("--mcounts", nargs=1, help = "maternal gene expression matrix; rows:genes x cols:samples/cells. Values should be normalized counts.")
parser$add_argument("--pcounts", nargs=1, help = "paternal gene expression matrix; rows:genes x cols:samples/cells. Values should be normalized counts.")

# optional arguments if looking at variability in imprinting
parser$add_argument("--method", default = 'combined', help = "Method to use for analysis, either 'combined' or 'separate'. See summary above. Only required if looking at imprinting.")

# optional settings for factors and scaling
parser$add_argument("--sampfile", help = "(NO HEADER!) list of samples (columns in input matrices) to use in plot (all other columns in expr_data will be omitted). If this file has two columns, second column indicates separate groups; see description above.", type="character")
parser$add_argument("--factors", help = "(NO HEADER!) list of samples with columns indicating additional factors (e.g. tissue). Basically the same as --sampfile, except --sampfile has the factor you're interested in, while --factors has any others you want to control for. Can have as many columns (factors) as you want.", type="character")
parser$add_argument("--genefile", help = "(NO HEADER!) File with a list of genes; if this file has two columns, second column assumed to provide 'gene groups'.", type="character")
parser$add_argument("--nolog", default=FALSE, action="store_true", help = "Do -not- take log2 of expr_matrix values before calculating averages (default is to take log2(x+1)). Ignored if analyzing imprinting with --method 'combined'.")
parser$add_argument("--pseudocount", default=1, help = "If --nolog is not provided, values are log2-transformed; this value is added as a pseudocount to avoid taking log(0)", type="double")
parser$add_argument("--allowmissinggenes", default=FALSE, action="store_true", help = "Let's the script continue if not all genes in --genefile are found in the input matrices")

# optional settings for clustering
parser$add_argument("--kmeans", default=FALSE, action="store_true", help = "Use k-means instead of hierarchical clustering")
parser$add_argument("--k", help = "Number of clusters k to use with k-means algorithm, or clusters to obtain from hierarchical clustering", type="integer")
parser$add_argument("--kmax", default = 100, help = "Upper limit on number of clusters k tested, if optimizing value of k (ignored if --k set)", type="integer")
parser$add_argument("--clustersamples", default=FALSE, action="store_true", help = "Cluster samples/cells (columns) in plot (--breakpoints will be ignored; method will always be hierarchical")
parser$add_argument("--hmethod", default = "complete", help = "Method to use for hierarchical clustering (see hclust() method option)")
parser$add_argument("--hdist", default = "euclidean", help = "Distance metric to use for hierarchical clustering (see dist() method option)")
parser$add_argument("--clustbypval", default=FALSE, action="store_true", help = "(Ignored if --sampfile has only one column, since Z-scores/P-vals won't be calculated) Cluster p-values instead of z-scores")
parser$add_argument("--clustbyobs", default=FALSE, action="store_true", help = "(Ignored if --sampfile has only one column, since Z-scores/P-vals won't be calculated) Cluster raw values instead of z-scores")
parser$add_argument("--NAreplacedist", help = "Replace NaN values in Z-score matrix with this value ***only for computing gene x gene distance matrix*** (will still appear as missing in heatmaps); NaNs occur when a gene is zero in all nuclei/samples within group - set this to whatever default value (e.g. 0 for z-scores, 0.67 if using observed values, etc.) makes sense", type="double")
parser$add_argument("--customxorder", default=FALSE, action="store_true", help = "Use order of rows in --genefile for plots (note that if used, all other clustering will not be performed). Compatible with all 'roworder' files output by this script (blanks included).")
parser$add_argument("--customyorder", help = "Use order of columns in this file for plots (note that if used, all other clustering will not be performed). Compatible with all 'colorder' files output by this script.")

# optional settings for output plots and files
parser$add_argument("--colorder", help = "Desired order of columns, provided as a comma-separated list (only used if --clustersamples and --customyorder not used; use --customyorder to specify desired column order from a file instead)")
parser$add_argument("--colorsZ", default = "blue,white,red", help = "Color scheme for heatmap z-scores and p-values, from lowest to highest value, provided as comma-separated list of 3+ colors.")
parser$add_argument("--colorsM", default = "white,red", help = "Color scheme for heatmap unscaled expression values, from lowest to highest value, provided as comma-separated list of 2 colors.")
parser$add_argument("--colorsD", default = "blue,white,red", help = "Color scheme for heatmap of 'D' scores (only applies if in imprinting mode with method == separate), from lowest to highest value, provided as comma-separated list of 3+ colors.")
parser$add_argument("--breakpoints", help = "Number(s) of one or more columns after which to insert an empty (white) column, as a break between groups")
parser$add_argument("--clustersep", default = 1, help = "Number of blank rows to add between groups of genes in the same cluster, to distinguish clusters", type="integer")
parser$add_argument("--width", default = 8.5, help = "Width of plot (in inches)", type="integer")
parser$add_argument("--height", default = 8, help = "Height of plot (in inches)", type="integer")
parser$add_argument("--plotupperZ", help = "Max abs() value of z-score for plot colors", type="integer")
parser$add_argument("--plotupperM", help = "Max abs() value of absolute expression for plot colors", type="integer")
parser$add_argument("--plotupperD", help = "Max abs() value of difference between mat. and pat. values (difscore plots)", type="integer")
parser$add_argument("--NAcolor", default = "white", help = "Color for missing values (including the rows/columns inserted by this script)", type="character")
parser$add_argument("--showrownames", default=FALSE, action="store_true", help = "Show row (gene) names in plot")
parser$add_argument("--outputmatrices", default=FALSE, action="store_true", help = "Output final scaled and unscaled matrices")

# optional settings for identifying DE genes
parser$add_argument("--nreps", default = 1000, help = "Number of times to shuffle cluster labels (in --sampfile), if applicable", type="integer")
parser$add_argument("--pval", default=0.05, help = "P-value cutoff (based on shuffling labels -nreps- times) for gene to be significantly enriched in a cluster", type="double")
parser$add_argument("--minsd", default=0, help = "Min. std. devation for distribution of values for a given cluster, in order to consider it significant (must also pass p-value cutoff; this is for dropping 'significant but not substantial' cases)", type="double")
parser$add_argument("--mincov", default=0, help = "('imprinting' mode only, only if sampfile provided) Minimum per-cluster coverage required, in order to consider a given cluster significant (must pass all other cutoffs)", type="integer")
parser$add_argument("--censorlowsd", default=FALSE, action="store_true", help = "Replace Z-scores, P-scores and sig. values by 'NA' if s.d. for a particular cell is below --minsd (e.g. if s.d. is so low that no meaningful conclusion can be reached, again this is for dropping 'significant but not substantial' cases)")

# additional options
parser$add_argument("--seed", default = 123456, help = "Set random number generator seed; runs with same seed and other inputs will always produce same results.", type="integer")

opt <- parser$parse_args()

set.seed(opt$seed)

if (! is.null(opt$breakpoints)) breakpoints = unlist(strsplit(opt$breakpoints,","))
colorsZ = unlist(strsplit(opt$colorsZ,","))
colorsM = unlist(strsplit(opt$colorsM,","))
colorsD = unlist(strsplit(opt$colorsD,","))
if (length(colorsZ) < 3) stop("--colorsZ must have 3 or more colors")
if (length(colorsZ) %% 2 == 0) stop("--colorsZ must have an odd number of colors")
if (length(colorsD) < 3) stop("--colorsD must have 3 or more colors")
if (length(colorsD) %% 2 == 0) stop("--colorsD must have an odd number of colors")

if (opt$clustbypval == TRUE && opt$clustbyobs == TRUE) stop("Must choose either --clustbypval or --clustbyobs, but not both")
if (! is.null(opt$customyorder) && opt$clustersamples == TRUE) stop("Must choose either --customyorder or --clustersamples, but not both")


# ------------------
# functions
# ------------------
# this function adapted from stackoverflow user dcarlson (https://stackoverflow.com/questions/57974376/how-to-compute-total-within-sum-of-square-in-hierarchical-clustering)
# x = the matrix, h = hclust() object, k = # of clusters to cut tree into
TSS <- function(x, h, k) {
	rr = cutree(h, k = k)
    sum(aggregate(x, by=list(rr), function(x) sum(scale(x,scale=FALSE)^2, na.rm = TRUE))[, -1], na.rm = TRUE)
}

# takes an NxM matrix Mat and a Mx1 dataframe Fac, where colnames(Mat) == rownames(Fac), and merges
# the two in order to collapse Mat by the factor in Fac. Note that if Fac has 2 columns, 2nd column
# is considered groups over which to shuffle (if shuffling, ignored otherwise).
collapseMat <- function(Mat, Fac, shuffle = FALSE, fun = mean) {

	if (ncol(Fac) != 1 && ncol(Fac) != 2) stop("in collapseMat(), Fac must have 1 or 2 columns")
	
	if (! setequal(rownames(Fac),colnames(Mat))) {
		cat("Error in call to collapseMat(); rownames(Fac) != colnames(Mat).\n")
		cat("head(Fac):\n")
		print(head(Fac))
		cat("Mat[1:5,1:5]:\n")
		print(Mat[1:min(5,nrow(Mat)),1:min(5,ncol(Mat))])
		stop("exiting")
	}
	
	if (shuffle == TRUE) {
		if (ncol(Fac) == 2) {
			Fac$newlbl = ave(Fac[,1], Fac[,2], FUN = sample)
			Fac = Fac[,3,drop=FALSE]	
		} else {
			Fac$newlbl = sample(Fac[,1])
			Fac = Fac[,2,drop=FALSE]	
		}
	} else {
		if (ncol(Fac) == 2) {
#			cat("Warning: in collapseMat(), Fac has 2 columns, but shuffle not requested.\n")
			Fac = Fac[,1,drop=FALSE]
		}
	}
	
	Mat = merge(Fac, t(Mat), by=0)
	rownames(Mat) = Mat$Row.names; Mat$Row.names = NULL
	Mat = aggregate(Mat[, 2:ncol(Mat)], list(Mat[,1]), FUN=fun)
	rownames(Mat) = Mat$Group.1; Mat$Group.1 = NULL; Mat = t(Mat)
	return(Mat)
}

collapseMatImprHelper <- function(m, p, f) {
	t = m + p
	f$ID = ifelse(t == 0,0,f$ID)
	f$newlbl = ave(f[,1], f[,2], FUN = sample)
	f = f[,3,drop=FALSE]
	m = aggregate(m, f, FUN=sum); rownames(m) = m[,1]; m[,1] = NULL
	p = aggregate(p, f, FUN=sum); rownames(p) = p[,1]; p[,1] = NULL
	return(list(m,p))
}

# similar to function above, but specifically for imprinting (accepts two matrices, one of mat. and one 
# of pat. counts, and calculated weighted mean frac. maternal instead of mean total expression)
collapseMatImpr <- function(Mat, Pat, Fac, shuffle = FALSE) {

	if (ncol(Fac) != 1 && ncol(Fac) != 2) stop("in collapseMat(), Fac must have 1 or 2 columns")
	
	if (! setequal(rownames(Mat),rownames(Pat))) stop("Error in collapseMatImpr(); rownames of mat. and pat. matrices don't match")
	if (! setequal(colnames(Mat),colnames(Pat))) stop("Error in collapseMatImpr(); colnames of mat. and pat. matrices don't match")

	if (! setequal(rownames(Fac),colnames(Mat))) {
		cat("Error in call to collapseMat(); rownames(Fac) != colnames(Mat).\n")
		cat("head(Fac):\n")
		print(head(Fac))
		cat("Mat[1:5,1:5]:\n")
		print(Mat[1:min(5,nrow(Mat)),1:min(5,ncol(Mat))])
		stop("exiting")
	}
	
	# make sure Mat, Pat and Fac are all in same order
	Pat = Pat[rownames(Mat), colnames(Mat)]
	Fac = Fac[colnames(Mat),,drop=FALSE]
	
	if (shuffle == TRUE) {
		if (ncol(Fac) == 1) {
			Fac$ID = 1
		}		
		# shuffle only the samples with m + p > 0 		
		# set up using first row
		res = collapseMatImprHelper(Mat[1,],Pat[1,], Fac)
		newMat = as.data.frame(t(res[[1]]), stringsAsFactors=FALSE); rownames(newMat) = rownames(Mat)[1]
		newPat = as.data.frame(t(res[[2]]), stringsAsFactors=FALSE); rownames(newPat) = rownames(Mat)[1]	
			
		# for all remaining rows
		if (nrow(Mat) > 2) {
			for (rr in 2:nrow(Mat)) {
				res = collapseMatImprHelper(Mat[rr,],Pat[rr,], Fac)
				tmpMat = as.data.frame(t(res[[1]])); rownames(tmpMat) = rownames(Mat)[rr]
				tmpPat = as.data.frame(t(res[[2]])); rownames(tmpPat) = rownames(Mat)[rr]
				newMat = rbind(newMat,tmpMat); newPat = rbind(newPat,tmpPat)		
			}
		}

		# calculate % mat
		pmat = newMat / (newMat + newPat); pmat = as.matrix(pmat)		
	} else {
		Mat = aggregate(t(Mat), list(Fac[,1]), FUN=sum); rownames(Mat) = Mat[,1]; Mat[,1] = NULL
		Pat = aggregate(t(Pat), list(Fac[,1]), FUN=sum); rownames(Pat) = Pat[,1]; Pat[,1] = NULL
				
		# calculate % mat
		pmat = Mat / (Mat + Pat); pmat = t(pmat); pmat = as.matrix(pmat)	
	}
	return(pmat)
}

# add blank rows to matrix, to separate the different clusters
# matrix must have column named 'cluster' and be in correct row order for final plot (e.g. all members of a cluster together)
# or you'll get -a lot- of breaks
addhblanks <- function(x, sepval) {

	if (! "cluster" %in% colnames(x)) stop("for addhblanks() function, matrix must have column named 'cluster'")

	# add blank rows between each cluster, corresponding to cluster breaks
	rowlist = rownames(x)
	breakpoints = c(which(diff(x$cluster)!=0))
	newroworder = c(); startp = 1
	for (i in 1:length(breakpoints)) {
		for (j in 1:sepval) {
			x[nrow(x)+1,] <- NA
			rownames(x)[nrow(x)] = paste("blank",i,"_",j,sep='')
		}
		newroworder = c(newroworder,rowlist[startp:breakpoints[i]])
	
		for (j in 1:opt$clustersep) {
			newroworder = c(newroworder,paste("blank",i,"_",j,sep=''))
		}
		startp=as.numeric(breakpoints[i])+1
	}
	newroworder = c(newroworder,rowlist[startp:length(rowlist)])
	x = x[newroworder,]
	
	return(x)
}

# add blank columns to matrix, to separate specific groups of samples
# second variable 'breakpoints' indicates number of one or more columns after which to insert empty column(s)
addvblanks <- function(x, breaks) {
	
	collist = colnames(x)
	newcolorder = c(); startp = 1
	for (i in 1:length(breaks)) {
		x[[paste("blank",i,sep='')]] = NA
		newcolorder = c(newcolorder,collist[startp:breaks[i]],paste("blank",i,sep=''))
		startp=as.numeric(breaks[i])+1
	}
	newcolorder = c(newcolorder,collist[startp:length(collist)])
	x = x[,newcolorder]
	
	return(x)
}

# perform hierarchical clustering of the rows of a matrix, and cut dendogram according to optimal or custom k
hierclust <- function(matforclus, hdist, hmethod, kmax, k_opt = NULL) {
	# calculate distance matrix	
	dd = dist(matforclus, method = hdist)
	if (length(which(!is.finite(as.matrix(dd)))) > 0) {
		stop("NAs in gene x gene distance matrix, caused by genes not varying within factor groups provided. Try reducing factors or dropping genes that are always 0 within one or more levels of a factor. To replace NaNs in Z-score matrix with 0s, use --NAreplacedist (use cautiously).")
	}
	# perform clustering
	myhclust <- hclust(dd, method = hmethod)
	
	# cut hclust dendogram into k distinct clusters using either user-supplied k or 'optimal k' defined by elbow method
	if (is.null(k_opt)) {
		cat("   - Finding optimal k for cutree()...")
		# if k was not provided, find optimal k by testing all values up to kmax
		wss <- sapply(1:kmax, function(k){TSS(matforclus,myhclust,k)})
		wssd = diff(wss)

		# Use crude elbow method to ID # of clusters
		mwss = max(2,min(5,kmax/5))*max(wssd)
		wssd = wssd[wssd < mwss]
		k_opt = length(wssd)-1
		cat(" using k =",k_opt,'\n')
	} else {
		# k was provided, use that value
		cat(" using k =",k_opt,' (user-provided)\n')
	}
	if (k_opt < 2) { cat("   - Adjusting to k = 2 since k < 2 is unuseable\n"); k_opt = 2; }	
	
	# get optimal partition of the tree	
	kkres = cutree(myhclust, k = k_opt)
	cluslist = kkres[match(rownames(matforclus),names(kkres))]
	return(list(myhclust,cluslist))
}

# perform kmeans clustering of the rows of a matrix according to optimal or custom k
kclust <- function(matforclus, hdist, kmax, k_opt = NULL) {

	dd = dist(matforclus, method = hdist)
	if (length(which(!is.finite(as.matrix(dd)))) > 0) {
		stop("NAs in gene x gene distance matrix, caused by genes not varying within factor groups provided. Try reducing factors or dropping genes that are always 0 within one or more levels of a factor. To replace NaNs in Z-score matrix with 0s, use --NAreplacedist (use cautiously).")
	}

	if (is.null(k_opt)) {
		cat("   - Finding optimal k for kmeans()...\n")
		# if k was not provided, find optimal k by testing all values up to kmax
		wss <- sapply(1:kmax, function(k){kmeans(matforclus, k, nstart=50,iter.max = 15)$tot.withinss})
		wssd = diff(wss)

		# Use crude elbow method to ID # of clusters
		mwss = max(2,min(5,kmax/5))*max(wssd)
		wssd = wssd[wssd < mwss]
		k_opt = length(wssd)-1
		cat("   - Using k =",k_opt,'\n')
	} else {
		# k was provided, use that value
		cat("   - Using k =",k_opt,' (user-provided)\n')
	}
	if (k_opt < 2) { cat("   - Adjusting to k = 2 since k < 2 is unuseable\n"); k_opt = 2; }	
	
	# use k-means clustering with k_opt to get clusters
	kkres = kmeans(matforclus, k_opt, nstart=50,iter.max = 150)
	cluslist = kkres$cluster
	cluslist = cluslist[match(rownames(matforclus),names(cluslist))]
	return(cluslist)
}

getClusOrder <- function(matforclus, cluslist, vclust=NULL) {
	tt = aggregate(matforclus, list(cluslist), mean, na.rm=TRUE)
	tt$Group.1 = NULL
	if (! is.null(vclust)) {
		if (class(vclust) == "hclust") {	
			tt = tt[,vclust$order]
		} else if (class(vclust) == "numeric" || class(vclust) == "integer") {
			tt = tt[,vclust]
		}
	}
	
	idxlist = rep(-1,max(cluslist))
	maxlist = rep(-1,max(cluslist))
	for (i in 1:max(cluslist)) {
		idx = which(tt[i,] == max(tt[i,],na.rm=TRUE), arr.ind = TRUE)[1,2]
		idxlist[i] = idx
		maxlist[i] = tt[i,idx]
	}

	sorter = as.data.frame(cbind(idxlist,maxlist))
	sorter$origorder = 1:max(cluslist)
	sorter = sorter[order(sorter$idxlist,-sorter$maxlist),]
	sorter$neworder = 1:max(cluslist)
	clusterorder = sorter$origorder
	return(clusterorder)
}

addhblanksall <- function(M,Z,P=NULL,S=NULL) {
	Z = addhblanks(Z, opt$clustersep)
	M = addhblanks(M, opt$clustersep)		
	if (! is.null(P)) {
		P = addhblanks(P, opt$clustersep)
		S = addhblanks(S, opt$clustersep)
	}
	return(list(M,Z,P,S))
}

combMatPatMat <- function(Mat,Pat) {
	colnames(Mat) = paste(colnames(Mat),"_m",sep="")
	colnames(Pat) = paste(colnames(Pat),"_p",sep="")
	res = cbind(Mat,Pat)
	return(res)
}

sortkclus <- function(matForClus, cluslist, vclusttouse) {
	# get new order for row clusters, sorting by order of columns
	clusterorder = getClusOrder(matForClus, cluslist, vclusttouse)
	
	# update cluslist numbering accordingly
	cluslistupdated = as.data.frame(cluslist)
	cluslistupdated$newclus = match(cluslistupdated$cluslist,clusterorder)
	cluslistupdated = cluslistupdated[order(cluslistupdated$newclus),2,drop=FALSE]
	cluslistfinal = cluslistupdated[,1]
	names(cluslistfinal) = rownames(cluslistupdated)
	return(cluslistfinal)
}


# Main function that runs the first part of the analysis
# ------------------
# Required inputs:
# mtx = 	input matrix of expression values, or (if mm == 'impr') a list of 2 matrices corresponding to mat and pat values
# grp = 	a dataframe with rownames == colnames(mtx), and 2 columns: first equal to the 'true' cluster labels, and second equal to groups to shuffle within
# nn = 		number of iterations of shuffling to perform
# Optional inputs:
# mm = 		mode, either 'expr' or 'impr' (determines whether collapseMat() or collapseMatImpr() is used)
# ping = 	(optionally) list of integers; print a statement when iteration number matches value in ping
# minsd = 	minimum std. deviation that shuffled means must have for a particular value to be significant (must also meet other significance criteria)
# censorlowsd = whether to censor from output matrices any value where the s.d. of the shuffled values is below 'minsd' (see 'minsd' above)
# mincov = 	(only if mm == 'impr') minimum coverage required for value to be significant
# cutoffHi/cutoffLo = cutoffs for finding significant bias in either direction
# Output:
# M = average expression/%mat over the true clusters
# Z = z-score estimates for the value in M, relative to the mean and s.d. obtained when shuffling nn times
# P = p-value estimates, based on fraction out of nn of times shuffled value was above/below true value
# S = [-1,0,1] value indicating significant up/more mat (1) or down/more pat (-1) expression, or not sig. (0)
part1 <- function(mtx, grp, nn, mm='expr', ping = c(), censorlowsd = FALSE, minsd = 0, mincov = 0, cutoffHi = 0.975, cutoffLo = 0.025) {

	# get the true, observed means of expression/%mat over levels of the factor in --sampfile
	if (mm == "expr") { 
		M = collapseMat(mtx, grp)
	} else {
		if (length(mtx) != 2) stop("input for first argument in part1() must be a list of two elements, if in 'combined' imprinting mode")
		Am = mtx[[1]]; Ap = mtx[[2]]
		M = collapseMatImpr(Am, Ap, grp)
		TT = Am + Ap; TT = collapseMat(TT, grp, fun = sum)
	} 

	# now shuffle the cluster labels randomly & re-calculate means over each cluster
	for (i in 1:opt$nreps) {
		if (i %in% ping) cat("     - Running ",i,"th iteration...\n",sep='')
		if (mm == "expr") { 
			T = collapseMat(mtx, grp, shuffle = TRUE)
		} else {
			T = collapseMatImpr(Am, Ap, grp, shuffle = TRUE)		
		} 
		T = T[match(rownames(M),rownames(T)),match(colnames(M),colnames(T))]		# sort according to 'true' matrix (should be ok already, but just in case)

		# calculate significance (how often is the random mean > or < true mean?)
		if (exists("G")) {
			G = G + ifelse(T < M, 1, ifelse(T == M,0.5,0))
		} else {
			G = ifelse(T < M, 1, ifelse(T == M,0.5,0))
		}

		# get sum and sum of squares, to get mean and s.d. later
		if (exists("totS")) {
			totS = totS + ifelse(is.finite(T),T,0)
			totSS = totSS + (ifelse(is.finite(T),T,0))^2
			TC = TC + ifelse(is.finite(T),1,0)
		} else {
			totS = ifelse(is.finite(T),T,0)
			totSS = (ifelse(is.finite(T),T,0))^2
			TC = ifelse(is.finite(T),1,0)
		}	
	}
	
	# calculate Z-scores, p-scores (pval estimates) and significance
	Mmean = totS / TC
	Msd = sqrt((totSS / TC) - (Mmean)^2)
	
	# calculate Z-scores
	Z = (M - Mmean) / Msd
	if (nrow(which(is.infinite(Z),arr.ind=TRUE)) > 0) {
		cat("Warning: some estimates are missing due to incomplete shuffling of the samples; please repeat analysis with a higher --nreps value to reduce number of NAs in plots (current number is: ",nrow(which(is.infinite(Z),arr.ind=TRUE)),").\n",sep="")
		Z[! is.finite(Z)] = NaN			
	}
	
	# censor based on sd if requested
	if (censorlowsd == TRUE) {
		numNAs = sum(ifelse(! is.finite(Z),1,0))
		Z[Msd < minsd] = NA
		numNAsfilt = sum(ifelse(! is.finite(Z),1,0)); dd = numNAsfilt - numNAs
		if (dd > 0) cat(" - After censoring cells with low variance (s.d. < ",opt$minsd,"), ",dd," cells were changed to NAs\n",sep="")		
	}
	
	# get remaining matrices P and S
	P = G / TC																					# value between 0 and 1 that estimates p-value; high value (near 1) == sig. higher expression (or if imprinting, more maternal), near 0 == sig. lower expression / more paternal	
	S = ifelse(P > cutoffHi & Msd >= minsd, 1, ifelse(P < cutoffLo & Msd >= minsd, -1, 0))		# 1 if significantly higher expressed (or if imprinting, maternal), -1 if sig. low expressed/paternal, 0 if not significant, NA if could not eval.

	# censor based on coverage, if requested & in imprinting mode
	if (mm == "impr" && mincov > 0) {
		S = ifelse(TT >= mincov, S, 0)
	}	

	# censor any P and S values that are missing in Z (occurs when value doesn't vary over factor, or from sd cutoff above)
	P = ifelse(! is.finite(Z), NaN, P)	
	S = ifelse(! is.finite(Z), NaN, S)
	
	# make matrices of just the genes that are significant
	siggenelist = rownames(S)[rowSums(abs(S),na.rm = TRUE) > 0]
	sigS = S[rownames(S) %in% siggenelist,,drop=FALSE]
	sigP = P[rownames(P) %in% siggenelist,,drop=FALSE]
	sigZ = Z[rownames(Z) %in% siggenelist,,drop=FALSE]
	sigM = M[rownames(M) %in% siggenelist,,drop=FALSE]
	
	# return the four main matrices
	return(list(M,Z,P,S,sigM,sigZ,sigP,sigS,siggenelist))
}

# Specialized heatmap plotting function
# Required arguments:
# Mtx = matrix of values
# hclus = clustering of rows (genes) - either an hclust object, or a vector of cluster numbers (integer) with names() = row names in Mat
# vclus = clustering of cols (samples/clusters) - either an hclust object, or a vector of cluster numbers (integer) with names() = column names in Mat
# outfilename = name for output file

# Optional arguments:
# colorlist = colors for heatmap
# vbreaks = custom set of breaks to introduce between columns (e.g. c(1,4) will insert a blank column after cols #1 and #4)
# width = plot width
# height = plot height
# rowannot = annotations for rows (genes), if any
# paletteLength = resolution of heatmap color scheme
# scaleupper = max value of scale (if -1, set to max abs. value in Mtx; this is default)
# scaletype = 'sym' -> color scale is symmetric around zero, 'pos' -> color scale goes from 0 to scaleupper, 'imp' -> imprinting %mat scale (goes from [0,1], with middle color at 0.67)
# showRowNames = whether to show row names in plot
# fontsize = fontsize for plot
# nacolor = color for missing values in plot
plotheatmap <- function(Mtx, hclus, vclus, outprefix, colorlist = c('blue','white','red'), vbreaks = NULL, clustersep = 1, orderByCustom = FALSE, width = 8.5, height = 8, rowannot = NULL, paletteLength = 50, scaleupper = -1, scaletype = 'sym', showRowNames = FALSE, fontsize = 10, nacolor = 'white', allowmissinggenes=FALSE) {

	# check inputs ok
	if (! scaletype %in% c('sym','pos','imp')) stop("Error in plotkmeans(), parameter 'scaletype' must have one of these three values: 'sym','pos','imp', but has value",scaletype)
	if (! is.null(vbreaks) && class(vclus) == "hclust") cat("Warning: both an hclust object and vbreaks provided to plotheatmap(), vbreaks will be ignored\n")

	# get color palette
	hmcolorscale = colorRampPalette(colorlist)(paletteLength)
	if (scaleupper > 0) {
		maxval = scaleupper
	} else {
		maxval = max(abs(Mtx),na.rm=TRUE)
	}
	
	if (scaletype == 'sym') {
		myBreaks <- c(seq(-1*maxval, maxval, length.out=paletteLength))
	} else if (scaletype == 'pos') {
		myBreaks <- c(seq(0, maxval, length.out=paletteLength))
	} else {
		myBreaks = c(seq(0,0.6667,length.out = (paletteLength/2)+1),seq(0.6667,1,length.out = paletteLength/2)[2:(paletteLength/2)])
	}
	
	# order rows
	if (class(hclus) == "hclust") {
		# hclus is a dendogram, use it directly
		hclustouse = hclus
		write.table(rownames(Mtx)[hclustouse$order], file=paste(outprefix,"_roworder.txt",sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names = FALSE)
	} else if (class(hclus) == "numeric" || class(hclus) == "integer") {
		# hclus is a vector of cluster assignments
		
		# may contain 'blankX' where X is a number; those are ok (and will produce blanks)
		hclus_true = hclus[! grepl("blank", names(hclus))]	
		if (allowmissinggenes == FALSE) {
			if (length(hclus_true) != nrow(Mtx)) stop("error in plotheatmap(), hclus length != number of rows in Mtx")
			if (! setequal(names(hclus_true),rownames(Mtx))) stop("error in plotheatmap(), hclus values don't match row names in Mtx")
		}
		hclus = hclus[order(hclus)]
		Mtx = Mtx[match(names(hclus),rownames(Mtx)),]
		if (clustersep > 0) {
			Mtx$cluster = hclus
			Mtx = addhblanks(Mtx, clustersep)
			hclustoutput = Mtx$cluster
			names(hclustoutput) = rownames(Mtx)
			Mtx$cluster = NULL
		} else {
			hclustoutput = hclus
		}
		
		write.table(hclustoutput, file=paste(outprefix,"_roworder.txt",sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names = TRUE)
		hclustouse = FALSE
	} else {
		hclustouse = FALSE
		write.table(rownames(Mtx), file=paste(outprefix,"_roworder.txt",sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names = TRUE)
	}
	
	# order columns
	if (class(vclus) == "hclust") {
		# vclus is a dendogram, use it directly
		vclustouse = vclus
		write.table(colnames(Mtx)[vclustouse$order], file=paste(outprefix,"_colorder.txt",sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names = FALSE)
	} else if (class(vclus) == "numeric" || class(vclus) == "integer") {
		# vclus is a vector of cluster assignments
		if (length(vclus) != ncol(Mtx)) stop("error in plotheatmap(), vclus length != number of cols in Mtx")
		if (! setequal(names(vclus),colnames(Mtx))) stop("error in plotheatmap(), vclus values don't match column names in Mtx")
		vclus = vclus[order(vclus)]
		Mtx = Mtx[,match(names(vclus),colnames(Mtx))]
		vclustouse = FALSE
		write.table(colnames(Mtx), file=paste(outprefix,"_colorder.txt",sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names = FALSE)
	} else {
		vclustouse = FALSE
		write.table(colnames(Mtx), file=paste(outprefix,"_colorder.txt",sep=''), sep='\t', quote=FALSE, col.names=FALSE, row.names = FALSE)
	}
	
	# add breaks between columns if requested
	if (! is.null(vbreaks) && class(vclus) != "hclust") {
		Mtx = addvblanks(Mtx,vbreaks)
	}
	
	pdf(paste(outprefix,".pdf",sep=''), width = width, height = height)
	p1 = pheatmap(Mtx, cluster_rows = hclustouse, cluster_cols=vclustouse, annotation_row = rowannot, color = hmcolorscale, breaks=myBreaks, show_rownames=showRowNames, fontsize_row=fontsize, na_col=nacolor)
	show(p1)
	graphics.off()
}


# ------------------
# check inputs and read in data
# ------------------
if (! is.null(opt$expr_matrix) && ! is.null(opt$mcounts)) stop("Error: provide either --expr_matrix or --mcounts/--pcounts, not both.")
if (is.null(opt$expr_matrix) && is.null(opt$mcounts)) stop("Error: must provide either --expr_matrix or --mcounts/--pcounts.")
if (! is.null(opt$mcounts) && is.null(opt$pcounts)) stop("Error: if --mcounts is provided, must also provide --pcounts.")
if (is.null(opt$mcounts) && ! is.null(opt$pcounts)) stop("Error: if --pcounts is provided, must also provide --mcounts.")
if (! is.null(opt$expr_matrix)) { mmode = "expression" } else { mmode = "imprinting" }
if (mmode == "imprinting" && opt$method == "combined") opt$nolog = TRUE
if (mmode == "expression") opt$method = "NA"
if (opt$customxorder == TRUE) opt$clustersep = 0
	
# print summary of params
cat("\n")
cat("\nRunning cluster_gene_expression.R v.2.0 (12/18/2020)\n")
cat("-----------------------\n")
cat("Current directory:",getwd(),"\n")
if (mmode == "expression") { 
	cat("Expression matrix:",opt$expr_matrix,"\n") 
} else { 
	cat("Matrix of (norm.) mat. counts:",opt$mcounts,"\n")
	cat("Matrix of (norm.) pat. counts:",opt$pcounts,"\n")
	cat("Analysis method:",opt$method,"\n")
}
cat("Output directory:",opt$outprefix,"\n")
cat("Random seed:",opt$seed,"\n")
cat("-----------------------\n")
cat("Additional options for factors and scaling:\n")
if (! is.null(opt$sampfile)) cat(" - Samples/clusters provided via --sampfile:",opt$sampfile,"\n")
if (! is.null(opt$factors)) cat(" - Factors provided via --factors:",opt$factors,"\n")
if (! is.null(opt$genefile)) cat(" - List of genes to use provided via --genefile:",opt$genefile,"\n")
cat(" - Taking log2(expr + 1) for each value in input matrix before calculations:",(! opt$nolog),"\n")
if (! opt$nolog) cat(" - Pseudocount for log2() transformation:",opt$pseudocount,"\n")
cat("-----------------------\n")
cat("Additional options for clustering:\n")
if (opt$kmeans == TRUE) { 
	cat(" - Performing k-means clustering\n") 
	if (exists("opt$k")) {
		cat(" - Value of k fixed to:",opt$k,"\n")
	} else {
		cat(" - Value of k will be determined using elbow method, with max allowed k:",opt$kmax,"\n")
	}
} else { 
	cat(" - Performing hierarchical clustering\n")
	cat(" - Distance metric to use for clustering:",opt$hdist,"\n")
	cat(" - Linkage method for clustering:",opt$hmethod,"\n")
}
if (opt$clustersamples == TRUE) cat(" - Also clustering columns (samples/clusters) using hierarchical clustering\n")
cat("-----------------------\n")
cat("Additional options for identifying genes with sig. enriched/depleted expression:\n")
cat(" - Number of times to shuffle labels to estimate distribution of expression vals:",opt$nreps,"\n")
cat(" - P-value cutoff for significance (1/2 of this applied to each tail, so this value is the alpha/false pos rate:",opt$pval,"\n")
cat("-----------------------\n")

# read in and check all input files
cat("\nReading in and checking input files...\n")

# check input matrix/matrices ok
if (mmode == "expression") {
	if(isFALSE(file.exists(opt$expr_matrix))) stop("cannot open expr_matrix ",opt$expr_matrix)	
	A = read.table(opt$expr_matrix, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
	A = as.matrix(A); ngenes = nrow(A); nsamps = ncol(A)
	cat(" - Input expr_matrix has",ngenes,"genes (rows) x",nsamps,"samples (columns)\n")
} else {
	if(isFALSE(file.exists(opt$mcounts))) stop("cannot open mcounts matrix ",opt$mcounts)	
	if(isFALSE(file.exists(opt$pcounts))) stop("cannot open pcounts matrix ",opt$pcounts)	

	Am = read.table(opt$mcounts, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
	Ap = read.table(opt$pcounts, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)

	Am = as.matrix(Am); Ap = as.matrix(Ap); ngenes = nrow(Am); nsamps = ncol(Am)
	if (! setequal(rownames(Am),rownames(Ap))) stop("--mcounts and --pcounts files don't have same rows (genes); these two files must have the same rows (genes) x samples (cols)")
	if (! setequal(colnames(Am),colnames(Ap))) stop("--mcounts and --pcounts files don't have same cols (samples); these two files must have the same rows (genes) x samples (cols)")
	Ap = Ap[rownames(Am),colnames(Am)]
	
	cat(" - Input mcounts/pcounts matrices have",ngenes,"genes (rows) x",nsamps,"samples (columns)\n")
}

# check sampfile ok
if (! is.null(opt$sampfile)) {
	if(isFALSE(file.exists(opt$sampfile))) stop("cannot open --sampfile file ",opt$sampfile)	
	clusters = read.table(opt$sampfile, row.names=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
	if (ncol(clusters) > 1) stop("--sampfile file should have 1 or 2 columns; file",opt$sampfile,"has",ncol(clusters)+1,"columns.")
	if (length(rownames(clusters)) != length(unique(rownames(clusters)))) stop("samples list provided in --sampfile cannot contain repeated values")

	# subset A according to samples in this file
	if (mmode == "expression") {
		A = A[,colnames(A) %in% rownames(clusters)]
		if (ncol(A) < nrow(clusters)) stop("not all samples in --sampfile were found in input expr_matrix, aborting (does your file have a header? if so, removing it may fix this error)")
		if (ncol(A) != nsamps) { cat(" - Retaining only the",ncol(A),"samples present in --sampfile\n"); nsamps = ncol(A); }
	} else {
		Am = Am[,colnames(Am) %in% rownames(clusters)]
		Ap = Ap[,colnames(Ap) %in% rownames(clusters)]
		if (ncol(Am) < nrow(clusters)) stop("not all samples in --sampfile were found in input --mcounts/pcounts files, aborting (does your file have a header? if so, removing it may fix this error)")
		if (ncol(Am) != nsamps) { cat(" - Retaining only the",ncol(Am),"samples present in --sampfile\n"); nsamps = ncol(Am); }
	}
	
	if (! is.null(opt$customyorder)) {
		customyorder = read.table(opt$customyorder, row.names=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
		customyorder = rownames(customyorder)
		if (ncol(clusters) == 1) {
			cluslist = unique(clusters[,1])
			if (! setequal(cluslist,customyorder)) stop("Values in --customyorder file don't match values in --sampfile")
		} else {
			if (! setequal(rownames(clusters),customyorder)) stop("Values in --customyorder file don't match values in --sampfile")
		}
	}

	if (ncol(clusters) == 1) colnames(clusters) = c("cluster")
	if (ncol(clusters) == 0) clusters = NULL
} else {
	clusters = NULL
}

# check genefile ok
if (! is.null(opt$genefile)) {
	if(isFALSE(file.exists(opt$genefile))) stop("cannot open --genefile file ",opt$genefile)	
	genelist = read.table(opt$genefile, row.names=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
	if (length(rownames(genelist)) != length(unique(rownames(genelist)))) stop("samples list provided in --genelist cannot contain repeated values")
	if (opt$customxorder == TRUE) {
		genelist_true = genelist
		genelist = genelist_true[! grepl("blank", rownames(genelist_true)),,drop=FALSE]
	}
		
	# subset A according to genes in this file
	if (mmode == "expression") {
		A = A[rownames(A) %in% rownames(genelist),]
		if (opt$allowmissinggenes == FALSE) {
			if (nrow(A) < nrow(genelist)) stop("not all samples in --genelist were found in input expr_matrix, aborting (does your file have a header? if so, removing it may fix this error); you can also override this error by adding --allowmissinggenes to your command.")
		} else {
			cat("Warning:",nrow(genelist) - nrow(A),"genes in --genelist were not found in expr_matrix and will be ignored, since --allowmissinggenes is set to true.\n")
		}
		if (nrow(A) != ngenes) { cat(" - Retaining only the",nrow(A),"genes present in --genelist\n"); ngenes = nrow(A); }
	} else {
		Am = Am[rownames(Am) %in% rownames(genelist),]
		Ap = Ap[rownames(Ap) %in% rownames(genelist),]
		if (opt$allowmissinggenes == FALSE) {
			if (nrow(Am) < nrow(genelist)) stop("not all samples in --genelist were found in input --mcounts/pcounts files, aborting (does your file have a header? if so, removing it may fix this error); you can also override this error by adding --allowmissinggenes to your command.")
		} else {
			cat("Warning:",nrow(genelist) - nrow(Am),"genes in --genelist were not found in --mcounts/--pcounts and will be ignored, since --allowmissinggenes is set to true.\n")
		}
		if (nrow(Am) != ngenes) { cat(" - Retaining only the",nrow(Am),"genes present in --genelist\n"); ngenes = nrow(Am); }
	}
}

# check factors file ok
if (! is.null(opt$factors)) {
	if(isFALSE(file.exists(opt$factors))) stop("cannot open --factors file ",opt$factors)	
	factors = read.table(opt$factors, row.names=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
	if (length(rownames(factors)) != length(unique(rownames(factors)))) stop("samples list (first column) provided in --factors cannot contain repeated values")
	cat(" -",ncol(factors),"factor(s) provided with --factors\n")
	nsamples = nrow(factors)
	if (! is.null(opt$sampfile)) {
		factors = factors[rownames(factors) %in% rownames(clusters),,drop=FALSE]
		if (nrow(clusters) != nrow(factors)) stop("not all samples in --sampfile found in --factors file\n")
		if (nsamples != nrow(factors)) cat(" - Censoring",nsamples - nrow(factors),"samples in --factors file not present in --sampfile\n")	
	} else {
		# drop any columns in A not in factors
		if (mmode == "expression") {
			A = A[,colnames(A) %in% rownames(factors)]
			if (ncol(A) < nrow(factors)) stop("not all samples in --factors were found in input expr_matrix, aborting (does your file have a header? if so, removing it may fix this error)")
			if (ncol(A) != nsamps) { cat(" - Retaining only the",ncol(A),"samples present in --factors\n"); nsamps = ncol(A); }		
		} else {
			Am = Am[,colnames(Am) %in% rownames(factors)]
			Ap = Ap[,colnames(Ap) %in% rownames(factors)]
			if (ncol(Am) < nrow(factors)) stop("not all samples in --factors were found in input expr_matrix, aborting (does your file have a header? if so, removing it may fix this error)")
			if (ncol(Am) != nsamps) { cat(" - Retaining only the",ncol(Am),"samples present in --factors\n"); nsamps = ncol(Am); }		
		}
	}
	
	factors = factors[do.call(order,factors),,drop=FALSE]
	factors$ID = cumsum(!duplicated(factors[1:ncol(factors)]))		# add column that flags all the unique combos of factors with a different ID #
}

# drop all rows in the matrix that are all missing data
if (mmode == "expression") {
	A = A[rowSums(is.na(A)) < ncol(A),]
	if (nrow(A) != ngenes) { cat("Censored",ngenes-nrow(A),"genes (rows) containing all missing data\n"); ngenes = nrow(A); }
	# drop these also from --genefile if provided
	if (! is.null(opt$genefile)) genelist = genelist[rownames(genelist) %in% rownames(A),,drop=FALSE]
} else {
	S = Am + Ap		# total allelic counts
	S = S[rowSums(is.na(S)) + rowSums(S == 0) < ncol(S),]
	if (nrow(S) != ngenes) { 
		cat("Censored",ngenes-nrow(S),"genes (rows) containing all missing data\n"); ngenes = nrow(S)
		Am = Am[rownames(Am) %in% rownames(S),]
		Ap = Ap[rownames(Ap) %in% rownames(S),]
		# drop these also from --genefile if provided
		if (! is.null(opt$genefile)) genelist = genelist[rownames(genelist) %in% rownames(Am),,drop=FALSE]
	}
}

# lower kmax to 1/2 the number of genes (if actual value of kmax is higher than that)
if (opt$kmax <= ngenes / 2) kmax = opt$kmax
if (opt$kmax > ngenes / 2 && opt$customxorder == FALSE) { kmax = floor(ngenes / 2); cat(" - lowering --kmax to",floor(ngenes / 2),"to ensure clustering algorithms will run\n") }

# for evaluating significance
cutoffLo = opt$pval/2; cutoffHi = 1 - cutoffLo

# fix colors for heatmaps if in imprinting mode
if (mmode == "imprinting" && opt$method == 'combined') {
	if (opt$colorsM == "white,red") {			# the default value
		colorsM = c('blue','white','red')
	} else {
		if (length(colorsM) != 3) stop("if evaluating imprinting (--mcounts/pcounts provided), must provide 3 colors to --colorsM")
	}
}

# ------------------
# main code
# ------------------

# (1) calculate 'scores' assessing expression/imprinting variability
# ----------------

# log-transform data before proceeding unless user requests otherwise
if (opt$nolog == FALSE) {
	cat(" - Taking log2 of data before averaging (e.g. mean will resemble geometric mean instead of arithmetic mean)\n")
	if (mmode == 'expression') {
		A = log2(A + opt$pseudocount)
	} else {
		Am = log2(Am + opt$pseudocount)
		Ap = log2(Ap + opt$pseudocount)
	}
}

# for iterations below, alert user when these iterations have been reached (notify every time 10% progress achieved)
ping = seq(1,opt$nreps,by=(opt$nreps / 10))

# if clusters were provided with --sampfile, get average expression over all nuclei in sample
# and calculate z-score of enrichment compared to if labels were randomly shuffled (within factors, if also provided)
if (! is.null(clusters)) {
	cat("\n")

	# get values to shuffle over
	if (! is.null(opt$factors)) {	
		shufgroups = merge(clusters,factors,by=0)
		rownames(shufgroups) = shufgroups$Row.names; shufgroups = shufgroups[,c(2,ncol(shufgroups))]
	} else {
		shufgroups = clusters
		shufgroups$ID = 1
	}

	# get the true, observed means of expression/%mat over levels of the factor in --sampfile
	if (mmode == "expression") { 
		cat("(1) Assessing expression variability across each cluster...\n")
		res = part1(A, shufgroups, opt$nreps, mm='expr', ping=ping, censorlowsd=opt$censorlowsd, minsd=opt$minsd, mincov=opt$mincov, cutoffHi=cutoffHi, cutoffLo=cutoffLo)
		M = res[[1]]; Z = res[[2]]; P = res[[3]]; S = res[[4]]
		sigM = res[[5]]; sigZ = res[[6]]; sigP = res[[7]]; sigS = res[[8]]; siggenelist = res[[9]]
	} else if (mmode == "imprinting" && opt$method == "combined") {
		cat("(1) Assessing % maternal variability across each cluster ('combined' method)...\n")
		res = part1(list(Am,Ap), shufgroups, opt$nreps, mm='impr', ping=ping, censorlowsd=opt$censorlowsd, minsd=opt$minsd, mincov=opt$mincov, cutoffHi=cutoffHi, cutoffLo=cutoffLo)
		M = res[[1]]; Z = res[[2]]; P = res[[3]]; S = res[[4]]
		sigM = res[[5]]; sigZ = res[[6]]; sigP = res[[7]]; sigS = res[[8]]; siggenelist = res[[9]]
	} else {
		cat("(1a) Assessing expression variability of maternal counts across each cluster ('separate' method)...\n")
		resm = part1(Am, shufgroups, opt$nreps, mm='expr', ping=ping, censorlowsd=opt$censorlowsd, minsd=opt$minsd, mincov=opt$mincov, cutoffHi=cutoffHi, cutoffLo=cutoffLo)
		Mm = resm[[1]]; Zm = resm[[2]]; Pm = resm[[3]]; Sm = resm[[4]]
		cat("(1b) Assessing expression variability of paternal counts across each cluster ('separate' method)...\n")
		resp = part1(Ap, shufgroups, opt$nreps, mm='expr', ping=ping, censorlowsd=opt$censorlowsd, minsd=opt$minsd, mincov=opt$mincov, cutoffHi=cutoffHi, cutoffLo=cutoffLo)
		Mp = resp[[1]]; Zp = resp[[2]]; Pp = resp[[3]]; Sp = resp[[4]]
		
		M = combMatPatMat(Mm,Mp); Z = combMatPatMat(Zm,Zp)
		P = combMatPatMat(Pm,Pp); S = combMatPatMat(Sm,Sp)

		siggenelist = unique(c(resm[[9]],resp[[9]]))		
		sigMm = Mm[rownames(Mm) %in% siggenelist,,drop=FALSE]; sigZm = Zm[rownames(Zm) %in% siggenelist,,drop=FALSE]
		sigPm = Pm[rownames(Pm) %in% siggenelist,,drop=FALSE]; sigSm = Sm[rownames(Sm) %in% siggenelist,,drop=FALSE]
		sigMp = Mp[rownames(Mp) %in% siggenelist,,drop=FALSE]; sigZp = Zp[rownames(Zp) %in% siggenelist,,drop=FALSE]
		sigPp = Pp[rownames(Pp) %in% siggenelist,,drop=FALSE]; sigSp = Sp[rownames(Sp) %in% siggenelist,,drop=FALSE]
		
		sigM = combMatPatMat(sigMm,sigMp); sigZ = combMatPatMat(sigZm,sigZp)
		sigP = combMatPatMat(sigPm,sigPp); sigS = combMatPatMat(sigSm,sigSp)		
	}		
} else {
	# if clusters not provided, just do regular normalization row-wise ((x - mean) / sd) via scale()
	if (mmode == "expression") { 
		cat("(1) Calculating average expression across samples\n")
		M = A
		Z = t(scale(t(M)))
	} else if (mmode == "expression" && opt$method == "combined") {
		cat("(1) Calculating average % maternal across samples\n")
		M = Am / (Am + Ap)
		Z = t(scale(t(M)))
	} else {
		cat("(1) Calculating average maternal and paternal expression across samples\n")
		Mm = Am; Zm = t(scale(t(Mm)))
		Mp = Ap; Zp = t(scale(t(Mp)))
		M = combMatPatMat(Mm,Mp); Z = combMatPatMat(Zm,Zp)
	}	
	siggenelist = c()
}

# write final matrices
cat("   - DONE, writing final scores to .txt files\n")
if (mmode == "expression" || (mmode == "imprinting" && opt$method == "combined")) { 
	write.table(M, file=paste(opt$outprefix,"_unscaled_expr.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
	write.table(Z, file=paste(opt$outprefix,"_zscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
	if (length(siggenelist) > 0) {
		write.table(P, file=paste(opt$outprefix,"_pscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)				
		write.table(S, file=paste(opt$outprefix,"_sig.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)	
	}	
} else {
	write.table(Mm, file=paste(opt$outprefix,"_mat_unscaled_expr.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
	write.table(Mp, file=paste(opt$outprefix,"_pat_unscaled_expr.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
	write.table(Zm, file=paste(opt$outprefix,"_mat_zscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
	write.table(Zp, file=paste(opt$outprefix,"_pat_zscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
	if (length(siggenelist) > 0) {
		write.table(Pm, file=paste(opt$outprefix,"_mat_pscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)				
		write.table(Pp, file=paste(opt$outprefix,"_pat_pscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)				
		write.table(Sm, file=paste(opt$outprefix,"_mat_sig.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
		write.table(Sp, file=paste(opt$outprefix,"_pat_sig.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
	}
}

# if factors provided, also examine variability across levels of each factor, holding all other factors constant (except for the --sampfile factor, which is not kept constant)
factorMs = list(); factorZs = list(); factorPs = list(); factorSs = list()
if (mmode == "imprinting" && opt$method == "separate") { 
	factorMms = list(); factorZms = list(); factorPms = list(); factorSms = list()
	factorMps = list(); factorZps = list(); factorPps = list(); factorSps = list()
}

if (! is.null(opt$factors)) {	
	cat(" - Repeating analysis over each factor provided in --factors file...\n")
			
	for (ff in 1:(ncol(factors)-1)) {	
		# resort order of columns in factors file so that ff-th column is last
		tmp = factors
		tmp$ID = NULL
		neworder = c(setdiff(1:(ncol(factors)-1),ff),ff)
		tmp = tmp[,neworder,drop=FALSE]
		tmp = tmp[do.call(order,tmp),,drop=FALSE]

		# shuffling will only be performed within groups of all other factors (excluding the factor in --sampfile)
		if (ncol(tmp) > 1) {
			tmp$ID = cumsum(!duplicated(tmp[1:(ncol(tmp)-1)]))
		} else {
			tmp$ID = 1
		}
		
		tmp = tmp[,c(ncol(tmp)-1,ncol(tmp))]	# keep only last two columns (the factor & ID)

		cat("   - Comparing factor with levels",paste(unique(tmp[,1]), collapse=", "),"\n")
		if (mmode == "expression") { 
			res = part1(A, tmp, opt$nreps, mm='expr', ping=ping, censorlowsd=opt$censorlowsd, minsd=opt$minsd, mincov=opt$mincov, cutoffHi=cutoffHi, cutoffLo=cutoffLo)
			factorMs[[ff]] = res[[1]]; factorZs[[ff]] = res[[2]]; factorPs[[ff]] = res[[3]]; factorSs[[ff]] = res[[4]]
		} else if (mmode == "imprinting" && opt$method == "combined") {
			res = part1(list(Am,Ap), tmp, opt$nreps, mm='impr', ping=ping, censorlowsd=opt$censorlowsd, minsd=opt$minsd, mincov=opt$mincov, cutoffHi=cutoffHi, cutoffLo=cutoffLo)
			factorMs[[ff]] = res[[1]]; factorZs[[ff]] = res[[2]]; factorPs[[ff]] = res[[3]]; factorSs[[ff]] = res[[4]]
		} else {
			cat("     - Analyzing maternal counts...\n")
			resm = part1(Am, tmp, opt$nreps, mm='expr', ping=ping, censorlowsd=opt$censorlowsd, minsd=opt$minsd, mincov=opt$mincov, cutoffHi=cutoffHi, cutoffLo=cutoffLo)
			factorMms[[ff]] = resm[[1]]; factorZms[[ff]] = resm[[2]]; factorPms[[ff]] = resm[[3]]; factorSms[[ff]] = resm[[4]]
			cat("     - Analyzing paternal counts...\n")
			resp = part1(Ap, tmp, opt$nreps, mm='expr', ping=ping, censorlowsd=opt$censorlowsd, minsd=opt$minsd, mincov=opt$mincov, cutoffHi=cutoffHi, cutoffLo=cutoffLo)
			factorMps[[ff]] = resp[[1]]; factorZps[[ff]] = resp[[2]]; factorPps[[ff]] = resp[[3]]; factorSps[[ff]] = resp[[4]]
			factorMs[[ff]] = combMatPatMat(resm[[1]],resp[[1]]); factorZs[[ff]] = combMatPatMat(resm[[2]],resp[[2]])
			factorPs[[ff]] = combMatPatMat(resm[[3]],resp[[3]]); factorSs[[ff]] = combMatPatMat(resm[[4]],resp[[4]])
		}		

		# write final matrices
		if (mmode == "expression" || (mmode == "imprinting" && opt$method == "combined")) { 
			cat("     - DONE - writing final scores to .txt files\n")
			write.table(factorMs[[ff]], file=paste(opt$outprefix,"_factor",ff,"_unscaled_expr.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
			write.table(factorZs[[ff]], file=paste(opt$outprefix,"_factor",ff,"_zscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
			write.table(factorPs[[ff]], file=paste(opt$outprefix,"_factor",ff,"_pscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)				
			write.table(factorSs[[ff]], file=paste(opt$outprefix,"_factor",ff,"_sig.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
		} else {
			write.table(factorMms[[ff]], file=paste(opt$outprefix,"_factor",ff,"_mat_unscaled_expr.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
			write.table(factorZms[[ff]], file=paste(opt$outprefix,"_factor",ff,"_mat_zscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
			write.table(factorPms[[ff]], file=paste(opt$outprefix,"_factor",ff,"_mat_pscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)				
			write.table(factorSms[[ff]], file=paste(opt$outprefix,"_factor",ff,"_mat_sig.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
			write.table(factorMps[[ff]], file=paste(opt$outprefix,"_factor",ff,"_pat_unscaled_expr.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
			write.table(factorZps[[ff]], file=paste(opt$outprefix,"_factor",ff,"_pat_zscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
			write.table(factorPps[[ff]], file=paste(opt$outprefix,"_factor",ff,"_pat_pscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)				
			write.table(factorSps[[ff]], file=paste(opt$outprefix,"_factor",ff,"_pat_sig.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)		
		}
	}
}

# (1.5) get matrices to use for clustering observations
# ----------------
if (! is.null(clusters)) {
	if (opt$clustbypval == FALSE && opt$clustbyobs == FALSE) { 
		clusbystr = "z-scores"
		matForClus = Z
		if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSig = sigZ; }
	} else if (opt$clustbypval == TRUE) { 
		clusbystr = "p-value estimates"
		matForClus = P
		if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSig = sigP; }
	} else if (opt$clustbyobs == TRUE) {
		clusbystr = "observed values"
		matForClus = M
		if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSig = sigM; }
	}
} else {
	clusbystr = "observed values"; matForClus = M;
}
if (mmode == "imprinting" && opt$method == "separate") {
	if (clusbystr == "z-scores") { 
		D = Zm - Zp
		matForClusD = Zm - Zp
		if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSigD = sigZm - sigZp; sigD = sigZm - sigZp; }
	} else if (clusbystr == "p-value estimates") {
		D = Pm - Pp
		matForClusD = Pm - Pp
		if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSigD = sigPm - sigPp; sigD = sigPm - sigPp; }
	} else {
		D = Mm - Mp
		matForClusD = Mm - Mp
		if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSigD = sigMm - sigMp; sigD = sigMm - sigMp; }
	}
	write.table(D, file=paste(opt$outprefix,"_dscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
	if(length(siggenelist) > 2 && opt$customxorder == FALSE) write.table(sigD, file=paste(opt$outprefix,"_sigonly_dscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
}

# replace missings by provided value, just for calculating distance matrix (will appear as missing in plots though)
if (! is.null(opt$NAreplacedist)) {
	numinf = sum(ifelse(! is.finite(matForClus), 1, 0))
	if (numinf > 0) cat(" - Warning:",numinf,"values in matrix (dimensions",dim(matForClus),") were replaced with",opt$NAreplacedist,"for calculating distance matrix\n")
	matForClus[! is.finite(matForClus)] = opt$NAreplacedist
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
		numinf = sum(ifelse(! is.finite(matForClusSig), 1, 0))
		if (numinf > 0) cat(" - Warning:",numinf,"values in sig. genes matrix (dimensions",dim(matForClusSig),") were replaced with",opt$NAreplacedist,"for calculating distance matrix\n")
		matForClusSig[! is.finite(matForClusSig)] = opt$NAreplacedist	
	}
	if (mmode == "imprinting" && opt$method == "separate") {
		numinf = sum(ifelse(! is.finite(matForClusD), 1, 0))
		if (numinf > 0) cat(" - Warning:",numinf,"values in D matrix (dimensions",dim(matForClusD),") were replaced with",opt$NAreplacedist,"for calculating distance matrix\n")
		matForClusD[! is.finite(matForClusD)] = 0
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
			numinf = sum(ifelse(! is.finite(matForClusSigD), 1, 0))
			if (numinf > 0) cat(" - Warning:",numinf,"values in sig. genes D matrix (dimensions",dim(matForClusSigD),") were replaced with",opt$NAreplacedist,"for calculating distance matrix\n")
			matForClusSigD[! is.finite(matForClusSigD)] = 0
		}
	}
}


# (2) cluster rows (genes)
# ----------------
if (opt$customxorder == FALSE) cat("(2) Clustering rows (genes)\n")
if (opt$customxorder == TRUE) cat("(2) Skipping clustering rows; using custom row (gene) order provided by user\n")

if (opt$customxorder == FALSE) {
	if (opt$kmeans == FALSE) {	
		# cluster rows using hierarchical clustering
		cat(" - Performing hierarchical clustering of rows(genes), using ",clusbystr,", distance metric ",opt$hdist,",and hclust() method ",opt$hmethod,"\n",sep='')	
		res = hierclust(matForClus, opt$hdist, opt$hmethod, kmax, k_opt = opt$k)
		myhclust = res[[1]]
		cluslist = res[[2]]
		
		# save result of clustering
		write.table(cluslist, file=paste(opt$outprefix,"_hier_cutree_clusters.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
		
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
			cat(" - Also clustering the subset of ",length(siggenelist)," rows(genes) with significant p-value\n",sep='')
			kmaxsig = max(2,min(floor(length(siggenelist)/2),kmax))
			if(! is.null(opt$k)) { koptsig = max(2,min(floor(length(siggenelist)/2),opt$k)) } else { koptsig = NULL }
			res = hierclust(matForClusSig, opt$hdist, opt$hmethod, kmaxsig, k_opt = koptsig)
			sigmyhclust = res[[1]]
			sigcluslist = res[[2]]
			write.table(sigcluslist, file=paste(opt$outprefix,"_hier_cutree_clusters_sig.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
		}
		if (mmode == "imprinting" && opt$method == "separate") {
			cat(" - Repeating hierarchical clustering on 'D' matrix (maternal - paternal ",clusbystr,")\n",sep='')	
			res = hierclust(matForClusD, opt$hdist, opt$hmethod, kmax, k_opt = opt$k)
			myhclustD = res[[1]]
			cluslistD = res[[2]]
			write.table(cluslistD, file=paste(opt$outprefix,"_hier_cutree_clusters_Dmat.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)

			if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
				cat(" - Repeating hierarchical clustering on 'D' matrix for the subset of ",length(siggenelist)," rows(genes) with significant p-value\n",sep='')	
				res = hierclust(matForClusSigD, opt$hdist, opt$hmethod, kmaxsig, k_opt = koptsig)
				sigmyhclustD = res[[1]]
				sigcluslistD = res[[2]]
				write.table(sigcluslistD, file=paste(opt$outprefix,"_hier_cutree_clusters_Dmat_sig.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
			}
		}
	} else {
		# cluster rows using k-means
		cat(" - Performing k-means clustering of rows(genes), using ",clusbystr,", and distance metric ",opt$hdist,"\n",sep='')	
		cluslist = kclust(matForClus, opt$hdist, kmax, k_opt = opt$k)	
		myhclust = FALSE
		kopt = max(cluslist)	
		write.table(cluslist, file=paste(opt$outprefix,"_k",kopt,"_clusters.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
				cat(" - Also clustering the subset of ",length(siggenelist)," rows(genes) with significant p-value\n",sep='')
				kmaxsig = max(2,min(floor(length(siggenelist)/2),kmax))
				if(! is.null(opt$k)) { koptsig = max(2,min(floor(length(siggenelist)/2),opt$k)) } else { koptsig = NULL }
				sigcluslist = kclust(matForClusSig, opt$hdist, kmaxsig, k_opt = koptsig)
				sigmyhclust = FALSE
				koptsig = max(sigcluslist)	
				write.table(sigcluslist, file=paste(opt$outprefix,"_k",kopt,"_clusters_sigonly.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
			} else {
				cat(" - Since only ",length(siggenelist)," rows(genes) had significant p-value, skipping clustering for this subset (will still be plotted)\n",sep='')
				sigcluslist = rep(1,length(siggenelist))
				names(sigcluslist) = siggenelist
			}
		}
		if (mmode == "imprinting" && opt$method == "separate") {
			cluslistD = kclust(matForClusD, opt$hdist, kmax, k_opt = opt$k)		
			myhclustD = FALSE
			write.table(cluslistD, file=paste(opt$outprefix,"_k",kopt,"_D_clusters.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
				if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
					sigcluslistD = kclust(matForClusSigD, opt$hdist, kmaxsig, k_opt = koptsig)
					sigmyhclustD = FALSE
					write.table(sigcluslistD, file=paste(opt$outprefix,"_k",kopt,"_D_clusters_sigonly.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
				} else {
					sigcluslistD = rep(1,length(siggenelist))
					names(sigcluslistD) = siggenelist
				}
			}
		}
	}
} else {
	myhclust = FALSE
	sigmyhclust = FALSE
	if (mmode == "imprinting" && opt$method == "separate") {
		myhclustD = FALSE
		sigmyhclustD = FALSE
	}
}


# (3) order columns (samples/clusters) as requested
# ----------------
cat("(3) Ordering columns (samples/clusters)\n")

M = as.data.frame(M); Z = as.data.frame(Z)
if (exists("P")) { P = as.data.frame(P); S = as.data.frame(S); }
if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
	sigM = as.data.frame(sigM); sigZ = as.data.frame(sigZ)
	sigP = as.data.frame(sigP); sigS = as.data.frame(sigS)
}
if (mmode == "imprinting" && opt$method == "separate") {
	D = as.data.frame(D)
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) sigD = as.data.frame(sigD)
}

if (opt$clustersamples == TRUE) {
	cat(" - Clustering columns (samples/clusters) according to",clusbystr,"using hierarchical clustering...\n")
	dd = dist(t(matForClus), method = opt$hdist)
	myvclust <- hclust(dd, method = opt$hmethod)	
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
		dd = dist(t(matForClusSig), method = opt$hdist)
		sigmyvclust <- hclust(dd, method = opt$hmethod)
	}
	if (mmode == "imprinting" && opt$method == "separate") {
		dd = dist(t(matForClusD), method = opt$hdist)
		myvclustD <- hclust(dd, method = opt$hmethod)	
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
			dd = dist(t(matForClusSigD), method = opt$hdist)
			sigmyvclustD <- hclust(dd, method = opt$hmethod)
		}
	}		
} else if (! is.null(opt$colorder)  && opt$clustersamples == FALSE) {
	cat(" - Sorting columns of matrices based on --colorder provided by user...\n")
	colorder = unlist(strsplit(opt$colorder,","))
	myvclust = 1:length(colorder)
	names(myvclust) = colorder
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) { sigmyvclust = myvclust }
	if (mmode == "imprinting" && opt$method == "separate") {
		myvclustD = myvclust
		tmp = as.vector(rbind(paste(names(myvclust),"_m",sep=''),paste(names(myvclust),"_p",sep='')))
		myvclust = 1:length(tmp); names(myvclust) = tmp
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
			sigmyvclustD = sigmyvclust 
			tmp = as.vector(rbind(paste(names(sigmyvclust),"_m",sep=''),paste(names(sigmyvclust),"_p",sep='')))
			sigmyvclust = 1:length(tmp); names(sigmyvclust) = tmp
		}
	}		
} else if (! is.null(opt$customyorder)) {
	cat(" - Sorting columns of matrices based on --customyorder provided by user...\n")
	myvclust = 1:length(customyorder)
	names(myvclust) = customyorder
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) { sigmyvclust = myvclust }
	if (mmode == "imprinting" && opt$method == "separate") {
		myvclustD = myvclust
		tmp = as.vector(rbind(paste(names(myvclust),"_m",sep=''),paste(names(myvclust),"_p",sep='')))
		myvclust = 1:length(tmp); names(myvclust) = tmp
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
			sigmyvclustD = sigmyvclust 
			tmp = as.vector(rbind(paste(names(sigmyvclust),"_m",sep=''),paste(names(sigmyvclust),"_p",sep='')))
			sigmyvclust = 1:length(tmp); names(sigmyvclust) = tmp
		}
	}
} else {
	myvclust = FALSE
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) { sigmyvclust = FALSE }
	if (mmode == "imprinting" && opt$method == "separate") {
		myvclustD = FALSE; sigmyvclustD = FALSE
	}
}
	

# (4) make heatmaps
# ----------------
cat("(4) Making heatmaps of unscaled average expression/%mat, Z-score and p-value estimates, and significance\n")

# get all plotting params
fontsize_all = max(10 - nrow(M) / 25, 3)
if (length(siggenelist) > 2 && opt$customxorder == FALSE) fontsize_sig = max(10 - nrow(sigM) / 25, 3)

if (ncol(genelist) > 0) { rowannot = genelist; } else { rowannot = NULL; }
paletteLength <- 50

if (opt$customxorder == FALSE) {
	cat(" - Making plots of genes divided into",max(cluslist),"clusters...\n")

	# order gene clusters by order of expression in columns
	if (length(myvclust) > 1) { vclusttouse = myvclust } else { vclusttouse = NULL }
	newcluslist = sortkclus(matForClus, cluslist, vclusttouse)
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
		if (length(sigmyvclust) > 1) { vclusttouse = sigmyvclust } else { vclusttouse = NULL }
		newsigcluslist = sortkclus(matForClusSig, sigcluslist, vclusttouse)
	}
	if (mmode == "imprinting" && opt$method == "separate") {
		if (length(myvclustD) > 1) { vclusttouse = myvclustD } else { vclusttouse = NULL }
		newcluslistD = sortkclus(matForClusD, cluslistD, vclusttouse)
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
			if (length(sigmyvclustD) > 1) { vclusttouse = sigmyvclustD } else { vclusttouse = NULL }
			newsigcluslistD = sortkclus(matForClusSigD, sigcluslistD, vclusttouse)
		}
	}
} else {	
	newcluslist = 1:nrow(genelist_true)
	names(newcluslist) = rownames(genelist_true)
	newcluslistD = newcluslist
}
	
# make all plots	
upperM = ifelse(is.null(opt$plotupperM),-1,opt$plotupperM)
upperZ = ifelse(is.null(opt$plotupperZ),-1,opt$plotupperZ)
upperD = ifelse(is.null(opt$plotupperD),-1,opt$plotupperD)

if (is.null(clusters)) {
	matrices = list(M, Z)
	colors = list(colorsM, colorsZ)
	maxvals = c(upperM, upperZ)
	outstr = c('unscaled','zscores')
	if (mmode == "imprinting" && opt$method == "combined") {
		scaletypes = c('imp','sym')
	} else {
		scaletypes = c('pos','sym')
	}
} else {
	matrices = list(M, Z, P, S)
	colors = list(colorsM, colorsZ, colorsZ, colorsZ)
	maxvals = c(upperM, upperZ, 1, 1)
	outstr = c('unscaled','zscores','pscores','sig')
	if (mmode == "imprinting" && opt$method == "combined") {
		scaletypes = c('imp','sym','pos','sym')
	} else {
		scaletypes = c('pos','sym','pos','sym')
	}
}

if (opt$customxorder == TRUE) {
	kstr = 'customx'
	kstrD = 'customx'
} else {
	kstr = paste('k',max(newcluslist),sep='')
	if (mmode == "imprinting" && opt$method == "separate") kstrD = paste('k',max(newcluslistD),sep='')
}

for (i in 1:length(matrices)) {
	plotheatmap(matrices[[i]], newcluslist, myvclust, paste(opt$outprefix,'_heatmap_',kstr,'_',outstr[i],sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
	if (opt$kmeans == FALSE && opt$customxorder == FALSE) plotheatmap(matrices[[i]], myhclust, myvclust, paste(opt$outprefix,'_heatmap_',outstr[i],sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
}

# repeat for sig subset, if exists
if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
	matrices = list(sigM, sigZ, sigP, sigS)
	for (i in 1:length(matrices)) {
		plotheatmap(matrices[[i]], newsigcluslist, sigmyvclust, paste(opt$outprefix,'_heatmap_k',max(newsigcluslist),'_sigonly_',outstr[i],sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_sig, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
		if (opt$kmeans == FALSE) plotheatmap(matrices[[i]], sigmyhclust, sigmyvclust, paste(opt$outprefix,'_heatmap_sigonly_',outstr[i],sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_sig, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
	}
}	
	
# repeat for the 'D' matrix, if exists
if (mmode == "imprinting" && opt$method == "separate") {	

	plotheatmap(D, newcluslistD, myvclustD, paste(opt$outprefix,'_heatmap_',kstrD,'_difscore',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
	if (opt$kmeans == FALSE && opt$customxorder == FALSE) plotheatmap(D, myhclustD, myvclustD, paste(opt$outprefix,'_heatmap_difscore',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
		plotheatmap(sigD, newsigcluslistD, sigmyvclustD, paste(opt$outprefix,'_heatmap_k',max(newsigcluslistD),'_sigonly_difscore',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)	
		if (opt$kmeans == FALSE) plotheatmap(sigD, sigmyhclustD, sigmyvclustD, paste(opt$outprefix,'_heatmap_sigonly_difscore',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)	
	}

	# also make version of plots showing individual mat and pat values, ordered according to D
	if (opt$clustersamples == TRUE) {
		newylbls = as.vector(rbind(paste(myvclustD$labels[myvclustD$order],"_m",sep=''),paste(myvclustD$labels[myvclustD$order],"_p",sep='')))
	} else if (class(myvclustD) != "logical") {
		newylbls = as.vector(rbind(paste(names(myvclustD),"_m",sep=''),paste(names(myvclustD),"_p",sep='')))
	} else {
		newylbls = as.vector(rbind(paste(colnames(D),"_m",sep=''),paste(colnames(D),"_p",sep='')))
	}
	
	if (! setequal(newylbls,colnames(M))) stop("internal error; column names of M != newylbls")
	newylblsclus = 1:length(newylbls)
	names(newylblsclus) = newylbls
	vbreaks = unname(newylblsclus[newylblsclus %% 2 == 0]); vbreaks = vbreaks[1:(length(vbreaks)-1)]		
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
		newylblssig = as.vector(rbind(paste(sigmyvclustD$labels[sigmyvclustD$order],"_m",sep=''),paste(sigmyvclustD$labels[sigmyvclustD$order],"_p",sep='')))
		if (! setequal(newylblssig,colnames(sigM))) stop("internal error; column names of sigM != newylblssig")
		newylblsclussig = 1:length(newylblssig)
		names(newylblsclussig) = newylblssig
		vbreakssig = unname(newylblsclussig[newylblsclussig %% 2 == 0]); vbreakssig = vbreakssig[1:(length(vbreakssig)-1)]
	}				

	if (length(siggenelist) <= 2) {
		matrices = list(M, Z)
	} else {
		matrices = list(M, Z, P, S)
	}
	if (length(siggenelist) <= 2) {
		matrices = matrices[1:2]
	}
	for (i in 1:length(matrices)) {
		plotheatmap(matrices[[i]], newcluslistD, newylblsclus, paste(opt$outprefix,'_heatmap_',kstrD,'_difscore_',outstr[i],sep=''), colorlist = colors[[i]], vbreaks = vbreaks, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
		if (opt$kmeans == FALSE && opt$customxorder == FALSE) plotheatmap(matrices[[i]], myhclustD, newylblsclus, paste(opt$outprefix,'_heatmap_difscore_',outstr[i],sep=''), colorlist = colors[[i]], vbreaks = vbreaks, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
	}
	if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
		matrices = list(sigM, sigZ, sigP, sigS)
		for (i in 1:length(matrices)) {
			plotheatmap(matrices[[i]], newsigcluslistD, newylblsclussig, paste(opt$outprefix,'_heatmap_k',max(newsigcluslistD),'_sigonly_difscore_',outstr[i],sep=''), colorlist = colors[[i]], vbreaks = vbreakssig, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
			if (opt$kmeans == FALSE) plotheatmap(matrices[[i]], sigmyhclustD, newylblsclussig, paste(opt$outprefix,'_heatmap_sigonly_difscore_',outstr[i],sep=''), colorlist = colors[[i]], vbreaks = vbreakssig, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
		}
	}
}

# repeat for all factors
if (! is.null(opt$factors)) {
	for (ff in 1:(ncol(factors)-1)) {
		cat(" - Making plots over levels of factor",ff,"...\n")
		M = as.data.frame(factorMs[[ff]]); Z = as.data.frame(factorZs[[ff]]); P = as.data.frame(factorPs[[ff]]); S = as.data.frame(factorSs[[ff]])
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
			sigS = S[rownames(S) %in% siggenelist,,drop=FALSE]; sigP = P[rownames(P) %in% siggenelist,,drop=FALSE]
			sigZ = Z[rownames(Z) %in% siggenelist,,drop=FALSE]; sigM = M[rownames(M) %in% siggenelist,,drop=FALSE]
		}		
		
		# get matrices for clustering
		if (opt$clustbypval == FALSE && opt$clustbyobs == FALSE) { 
			clusbystr = "z-scores"
			matForClus = Z
			if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSig = sigZ; }
		} else if (opt$clustbypval == TRUE) { 
			clusbystr = "p-value estimates"
			matForClus = P
			if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSig = sigP; }
		} else if (opt$clustbyobs == TRUE) {
			clusbystr = "observed values"
			matForClus = M
			if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSig = sigM; }
		}

		if (mmode == "imprinting" && opt$method == "separate") {
			Mm = as.data.frame(factorMms[[ff]]); Zm = as.data.frame(factorZms[[ff]]); Pm = as.data.frame(factorPms[[ff]]); Sm = as.data.frame(factorSms[[ff]])
			Mp = as.data.frame(factorMps[[ff]]); Zp = as.data.frame(factorZps[[ff]]); Pp = as.data.frame(factorPps[[ff]]); Sp = as.data.frame(factorSps[[ff]])
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
				sigSm = Sm[rownames(Sm) %in% siggenelist,,drop=FALSE]; sigPm = Pm[rownames(Pm) %in% siggenelist,,drop=FALSE]
				sigZm = Zm[rownames(Zm) %in% siggenelist,,drop=FALSE]; sigMm = Mm[rownames(Mm) %in% siggenelist,,drop=FALSE]
				sigSp = Sp[rownames(Sp) %in% siggenelist,,drop=FALSE]; sigPp = Pp[rownames(Pp) %in% siggenelist,,drop=FALSE]
				sigZp = Zp[rownames(Zp) %in% siggenelist,,drop=FALSE]; sigMp = Mp[rownames(Mp) %in% siggenelist,,drop=FALSE]
			}
			if (clusbystr == "z-scores") { 
				D = Zm - Zp
				matForClusD = D
				if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSigD = sigZm - sigZp; sigD = sigZm - sigZp; }
			} else if (clusbystr == "p-value estimates") {
				D = Pm - Pp
				matForClusD = Pm - Pp
				if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSigD = sigPm - sigPp; sigD = sigPm - sigPp; }
			} else {
				D = Mm - Mp
				matForClusD = Mm - Mp
				if(length(siggenelist) > 2 && opt$customxorder == FALSE) { matForClusSigD = sigMm - sigMp; sigD = sigMm - sigMp; }
			}
			write.table(D, file=paste(opt$outprefix,'_factor',ff,"_dscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
			if(length(siggenelist) > 2 && opt$customxorder == FALSE) write.table(sigD, file=paste(opt$outprefix,'_factor',ff,"_sigonly_dscores.txt",sep=''), sep='\t', quote=FALSE, col.names=TRUE, row.names = TRUE)
		}

		# replace missings by provided value, just for calculating distance matrix (will appear as missing in plots though)
		if (! is.null(opt$NAreplacedist)) {
			matForClus[! is.finite(as.matrix(matForClus))] = opt$NAreplacedist
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) matForClusSig[! is.finite(as.matrix(matForClusSig))] = opt$NAreplacedist
			if (mmode == "imprinting" && opt$method == "separate") {
				matForClusD[! is.finite(as.matrix(matForClusD))] = 0
				if (length(siggenelist) > 2 && opt$customxorder == FALSE) matForClusSigD[! is.finite(as.matrix(matForClusSigD))] = 0
			}
		}

		# cluster rows
		dd = dist(matForClus, method = opt$hdist)
		myhclust_f <- hclust(dd, method = opt$hmethod)	
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
			dd = dist(matForClusSig, method = opt$hdist)
			sigmyhclust_f <- hclust(dd, method = opt$hmethod)
		}
		if (mmode == "imprinting" && opt$method == "separate") {
			dd = dist(matForClusD, method = opt$hdist)
			myhclustD_f <- hclust(dd, method = opt$hmethod)	
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
				dd = dist(matForClusSigD, method = opt$hdist)
				sigmyhclustD_f <- hclust(dd, method = opt$hmethod)
			}
		}
		
		# cluster columns, if requested
		if (opt$clustersamples == TRUE) {
			dd = dist(t(matForClus), method = opt$hdist)
			myvclust <- hclust(dd, method = opt$hmethod)	
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
				dd = dist(t(matForClusSig), method = opt$hdist)
				sigmyvclust <- hclust(dd, method = opt$hmethod)
			}
			if (mmode == "imprinting" && opt$method == "separate") {
				dd = dist(t(matForClusD), method = opt$hdist)
				myvclustD <- hclust(dd, method = opt$hmethod)	
				if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
					dd = dist(t(matForClusSigD), method = opt$hdist)
					sigmyvclustD <- hclust(dd, method = opt$hmethod)
				}
			}		
		} else {
			myvclust = FALSE
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) { sigmyvclust = FALSE }
			if (mmode == "imprinting" && opt$method == "separate") {
				myvclustD = FALSE
				sigmyvclustD = FALSE
			}
		}
		
		# convert all matrices to data frames
		M = as.data.frame(M); Z = as.data.frame(Z)
		if (exists("P")) { P = as.data.frame(P); S = as.data.frame(S); }
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
			sigM = as.data.frame(sigM); sigZ = as.data.frame(sigZ)
			sigP = as.data.frame(sigP); sigS = as.data.frame(sigS)
		}
		if (mmode == "imprinting" && opt$method == "separate") {
			D = as.data.frame(D)
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) sigD = as.data.frame(sigD)
		}
				
		# make plots
		matrices = list(M, Z, P, S)
		colors = list(colorsM, colorsZ, colorsZ, colorsZ)
		upperM = ifelse(is.null(opt$plotupperM),-1,opt$plotupperM)
		upperZ = ifelse(is.null(opt$plotupperZ),-1,opt$plotupperZ)
		upperD = ifelse(is.null(opt$plotupperD),-1,opt$plotupperD)
		maxvals = c(upperM, upperZ, 1, 1)
		outstr = c('unscaled','zscores','pscores','sig')
		if (mmode == "imprinting" && opt$method == "combined") {
			scaletypes = c('imp','sym','pos','sym')
		} else {
			scaletypes = c('pos','sym','pos','sym')
		}
		for (i in 1:length(matrices)) {
			plotheatmap(matrices[[i]], newcluslist, myvclust, paste(opt$outprefix,'_heatmap_',kstr,'_factor',ff,'_',outstr[i],sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
			if (opt$customxorder == FALSE) plotheatmap(matrices[[i]], myhclust_f, myvclust, paste(opt$outprefix,'_heatmap_factor',ff,'_',outstr[i],'_reclus',sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
			if (opt$kmeans == FALSE && opt$customxorder == FALSE) plotheatmap(matrices[[i]], myhclust, myvclust, paste(opt$outprefix,'_heatmap_factor',ff,'_',outstr[i],sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
		}

		# repeat for sig subset, if exists
		if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
			matrices = list(sigM, sigZ, sigP, sigS)
			for (i in 1:length(matrices)) {
				plotheatmap(matrices[[i]], newsigcluslist, sigmyvclust, paste(opt$outprefix,'_heatmap_k',max(newsigcluslist),'_factor',ff,'_sigonly_',outstr[i],sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_sig, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
				plotheatmap(matrices[[i]], sigmyhclust_f, sigmyvclust, paste(opt$outprefix,'_heatmap_factor',ff,'_sigonly_',outstr[i],'_reclus',sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_sig, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
				if (opt$kmeans == FALSE) plotheatmap(matrices[[i]], sigmyhclust, sigmyvclust, paste(opt$outprefix,'_heatmap_factor',ff,'_sigonly_',outstr[i],sep=''), colorlist = colors[[i]], clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_sig, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
			}
		}	
	
		# repeat for the 'D' matrix, if exists
		if (mmode == "imprinting" && opt$method == "separate") {	
			plotheatmap(D, newcluslistD, myvclustD, paste(opt$outprefix,'_heatmap_',kstrD,'_factor',ff,'_difscore',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
			if (opt$customxorder == FALSE) plotheatmap(D, myhclustD_f, myvclustD, paste(opt$outprefix,'_heatmap_factor',ff,'_difscore_reclus',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
			if (opt$kmeans == FALSE && opt$customxorder == FALSE) plotheatmap(D, myhclustD, myvclustD, paste(opt$outprefix,'_heatmap_factor',ff,'_difscore',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) { 
				plotheatmap(sigD, newsigcluslistD, sigmyvclustD, paste(opt$outprefix,'_heatmap_k',max(newsigcluslistD),'_factor',ff,'_sigonly_difscore',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)	
				plotheatmap(sigD, sigmyhclustD_f, sigmyvclustD, paste(opt$outprefix,'_heatmap_factor',ff,'_sigonly_difscore_reclus',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)	
				if (opt$kmeans == FALSE) plotheatmap(sigD, sigmyhclustD, sigmyvclustD, paste(opt$outprefix,'_heatmap_factor',ff,'_sigonly_difscore',sep=''), colorlist = colorsD, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = upperD, scaletype = 'sym', showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)	
			}
			# also make version of plots showing individual mat and pat values, ordered according to D
			if (opt$clustersamples == TRUE) {
				newylbls = as.vector(rbind(paste(myvclustD$labels[myvclustD$order],"_m",sep=''),paste(myvclustD$labels[myvclustD$order],"_p",sep='')))
			} else if (class(myvclustD) != "logical") {
				newylbls = as.vector(rbind(paste(names(myvclustD),"_m",sep=''),paste(names(myvclustD),"_p",sep='')))
			} else {
				newylbls = as.vector(rbind(paste(colnames(D),"_m",sep=''),paste(colnames(D),"_p",sep='')))
			}
			if (! setequal(newylbls,colnames(M))) stop("internal error; column names of M != newylbls")
			newylblsclus = 1:length(newylbls)
			names(newylblsclus) = newylbls
			vbreaks = unname(newylblsclus[newylblsclus %% 2 == 0]); vbreaks = vbreaks[1:(length(vbreaks)-1)]		
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
				newylblssig = as.vector(rbind(paste(sigmyvclustD$labels[sigmyvclustD$order],"_m",sep=''),paste(sigmyvclustD$labels[sigmyvclustD$order],"_p",sep='')))
				if (! setequal(newylblssig,colnames(sigM))) stop("internal error; column names of sigM != newylblssig")
				newylblsclussig = 1:length(newylblssig)
				names(newylblsclussig) = newylblssig
				vbreakssig = unname(newylblsclussig[newylblsclussig %% 2 == 0]); vbreakssig = vbreakssig[1:(length(vbreakssig)-1)]
			}				

			matrices = list(M, Z, P, S)
			for (i in 1:length(matrices)) {
				plotheatmap(matrices[[i]], newcluslistD, newylblsclus, paste(opt$outprefix,'_heatmap_',kstrD,'_factor',ff,'_difscore_',outstr[i],sep=''), colorlist = colors[[i]], vbreaks = vbreaks, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
				if (opt$customxorder == FALSE) plotheatmap(matrices[[i]], myhclustD_f, newylblsclus, paste(opt$outprefix,'_heatmap_factor',ff,'_difscore_',outstr[i],'_reclus',sep=''), colorlist = colors[[i]], vbreaks = vbreaks, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
				if (opt$kmeans == FALSE && opt$customxorder == FALSE) plotheatmap(matrices[[i]], myhclustD, newylblsclus, paste(opt$outprefix,'_heatmap_factor',ff,'_difscore_',outstr[i],sep=''), colorlist = colors[[i]], vbreaks = vbreaks, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
			}
			if (length(siggenelist) > 2 && opt$customxorder == FALSE) {
				matrices = list(sigM, sigZ, sigP, sigS)
				for (i in 1:length(matrices)) {
					plotheatmap(matrices[[i]], newsigcluslistD, newylblsclussig, paste(opt$outprefix,'_heatmap_k',max(newsigcluslistD),'_factor',ff,'_sigonly_difscore_',outstr[i],sep=''), colorlist = colors[[i]], vbreaks = vbreakssig, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
					plotheatmap(matrices[[i]], sigmyhclustD_f, newylblsclussig, paste(opt$outprefix,'_heatmap_factor',ff,'_sigonly_difscore_',outstr[i],'_reclus',sep=''), colorlist = colors[[i]], vbreaks = vbreakssig, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
					if (opt$kmeans == FALSE) plotheatmap(matrices[[i]], sigmyhclustD, newylblsclussig, paste(opt$outprefix,'_heatmap_factor',ff,'_sigonly_difscore_',outstr[i],sep=''), colorlist = colors[[i]], vbreaks = vbreakssig, clustersep = opt$clustersep, orderByCustom = orderByCustom, width = opt$width*1.5, height = opt$height, rowannot = rowannot, paletteLength = 50, scaleupper = maxvals[i], scaletype = scaletypes[i], showRowNames = opt$showrownames, fontsize = fontsize_all, nacolor = opt$NAcolor, allowmissinggenes = opt$allowmissinggenes)
				}
			}
		}
	}
}




	
	
	
	
	
	
	
	
	
	
	
	













