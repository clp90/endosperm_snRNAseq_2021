#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(fpc))
suppressPackageStartupMessages(library(Rtsne))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(princurve))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(maptools))
suppressPackageStartupMessages(library(spatstat))
suppressPackageStartupMessages(library(RANN))
suppressPackageStartupMessages(require(pscl))
suppressPackageStartupMessages(require(MASS))
suppressPackageStartupMessages(require(boot))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(plyr))

# version 1.0 (05/25/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 05/25/2019
# -------------------------

# Description:
# -------------------------
# Broadly, this script is designed to try to fit a 'trajectory' (when the trajectory roughly
# corresponds to time, this can be called 'trajectory') curve to the single-cell
# data provided, using tSNE to get x-y values for each cell/nucleus and fitting a curve using
# the R package princurve. Broadly, tries to recapitulate part of the analysis shown in Fig. 2
# from this paper: https://science.sciencemag.org/content/364/6435/52. 

# The big parts of this analysis:
# (1) choosing tSNE parameters (perplexity and #PCs) - user must choose; run this script
# 	  with >1 value of --perplexity, --PCs or --k to inspect plots with all combinations of
#	  requested values; choose one set, and re-run script to do full analysis.
# (2) 

# Required inputs:
# -------------------------
# expr_matrix: This script accepts a matrix of raw counts (cols:cells x rows:genes) or pre-normalized 
# values (e.g. FPKM, etc.), filters based on user parameters, then performs requested 
# normalization, if any (--CPM or --FPKM; make sure not to request normalization if values 
# are already normalized)

# Notes:
# -------------------------




# Usage: single_cell_trajectory_analysis.R [options] expr_matrix.txt outprefix

parser = ArgumentParser()
parser$add_argument("expr_matrix", nargs=1, help = "tab-delimited text file matrix of raw expression counts (cols:cells x rows:genes)")
parser$add_argument("outprefix", nargs=1, help = "prefix for output files")
parser$add_argument("--seed", default = 1234567, help = "Set random number generator seed; runs with same seed and other inputs should always produce same results.", type="integer")
parser$add_argument("--minreads", default = 10, help = "Genes with fewer than this total value across all cells, are censored. Note that this filter is applied before any transformations (CPM, log2) requested.", type="integer")
parser$add_argument("--mincells", default = 5, help = "Genes detected in fewer than this many cells are censored.", type="integer")
parser$add_argument("--minexpr", default = 1000, help = "Cells detecting fewer than this many genes are censored.", type="integer")
parser$add_argument("--CPM", default=FALSE, action="store_true", help = "Normalize counts by CPM before building trajectory. Note that counts are used when testing for correlation with the trajectory, regardless of this setting.")
parser$add_argument("--log2", default=FALSE, action="store_true", help = "Take log2 of matrix values before building trajectory. Note that counts are used when testing for correlation with the trajectory, regardless of this setting.")
parser$add_argument("--pseudocount", default=1, help = "If taking log2 of values, this number is added first to avoid taking log2(0)", type="integer")
parser$add_argument("--genelist", help = "use only these genes to build trajectory; all genes are later analysed for differential expression across trajectory. Can contain additional columns with gene metadata (e.g. cell cycle phase the gene is expressed in, etc.) -> script will generate heatmap marking each of these gene groups")
parser$add_argument("--topN", default = 2000, help = "use the top --topN most variable genes for trajectory analysis; overridden by providing a --genelist")
parser$add_argument("--perplexity", default = "50", help = "Perplexity parameter for t-SNE; if comma-separated list provided (e.g. 10,20,50) will try all values", type="character")
parser$add_argument("--PCs", default = "3", help = "Must be ≤ 10. Use this many of the top PCs from regular PCA as input to tSNE; if comma-separated list provided (e.g. 2,3,5,10) will try all values", type="character")
parser$add_argument("--k", help = "Number of clusters k to look for in tSNE plot (default uses elbow method to find best number)", type="character")
parser$add_argument("--kmax", default = 10, help = "Upper limit on number of clusters k tested", type="integer")
parser$add_argument("--smoother", default = "smooth_spline", help = "'Smoother' for princurve to use, must be either 'smooth_spline', 'lowess', or 'periodic_lowess' (default smooth_spline). Or, use 'dijkstra' and I'll use a custom function - try this if the princurve options aren't working well.", type="character")
parser$add_argument("--neighbors", default = 5, help = "When building trajectories, consider this many nearest neighbors when finding path thru points", type="integer")
parser$add_argument("--fracdraw", default = 0.5, help = "When smoothing trajectories by bootstrapping, sample this fraction from each cluster before repeating k-means clustering + finding trajectory", type="double")
parser$add_argument("--numiter", default = 100, help = "Number of times to repeat analysis by bootstrapping; larger numbers will take longer but produce smoother curves", type="double")
parser$add_argument("--exact", default=FALSE, action="store_true", help = "Use a guaranteed-correct (but very slow for k > 5) method for finding shortest path through nodes (default uses slightly less accurate, much faster method)")
parser$add_argument("--sigthreshold", default = 0.01, help = "Significance threshold used to identify genes that vary significantly wrt the trajectory (header of this script has more info on models, etc.)", type="double")
parser$add_argument("--cell_info", help = "tab-delimited text file containing one row per cell + at least one descriptor (e.g. age, tissue, etc.)", type = "character")
parser$add_argument("--kgenes", default = 5, help = "tab-delimited text file containing one row per cell + at least one descriptor (e.g. age, tissue, etc.)", type = "integer")
parser$add_argument("--testing", default=FALSE, action="store_true", help = "A test run that just builds the consensus trajectory and reports gene expression of the genelist used (--genelist if provided) along the trajectory. Use to fine-tune your settings without running the long part.")

opt <- parser$parse_args()

set.seed(opt$seed)

fexpr_matrix = opt$expr_matrix
outprefix = opt$outprefix
minreads = opt$minreads
mincells = opt$mincells
minexpr = opt$minexpr
topN = opt$topN
takelog2 = opt$log2
pseudocount = opt$pseudocount
perplexity = as.numeric(strsplit(opt$perplexity,',')[[1]])
PCs = as.numeric(strsplit(opt$PCs,',')[[1]])
kmax = opt$kmax
CPM = opt$CPM
smoother = opt$smoother
fracdraw = opt$fracdraw
numiter = opt$numiter
neighbors = opt$neighbors
exact = opt$exact
sigthreshold = opt$sigthreshold
cell_info = opt$cell_info
kgenes = opt$kgenes
testing = opt$testing

if (length(opt$k) != 0) {
	karray = as.numeric(strsplit(opt$k,',')[[1]])
} else {
	karray = c()
}

if (smoother != 'dijkstra' && smoother != 'smooth_spline' && smoother != 'lowess' && smoother != 'periodic_lowess') {
	stop(smoother,' is not a valid value for the --smoother option (see usage)')
}

if (max(PCs) > 10) {
	stop("Error: the value of --PCs must be <= 10 (provided value was",PCs,")")
}

args <- commandArgs(trailingOnly = TRUE)	
if (length(opt) <= 1) {
	cat("Usage: single_cell_trajectory_analysis.R [options] expr_matrix.txt outprefix\n")
	cat("----------------------\n")
	cat("expr_matrix.txt : tab-delimited text file matrix of raw expression counts (cols:cells x rows:genes), or normalized values (e.g. FPKM)\n")
	cat("outprefix : prefix for output files\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("-TODO\n")
	cat("----------------------\n")
}


cat("\nRunning single_cell_trajectory_analysis.R v.1.0 (05/25/2019)\n")
cat("-----------------------\n")
cat("Current directory:",getwd(),"\n")
cat("Expression matrix:",fexpr_matrix,"\n")
cat("Prefix for output files:",outprefix,"\n")
cat("-----------------------\n")
cat("Transformation settings:\n")
if (CPM == TRUE) {
	cat("Transforming read counts into CPM\n")
}
if (takelog2 == TRUE) {
	cat("Taking log2 of values, using pseudocount of",pseudocount,"\n")
}
cat("-----------------------\n")
cat("Quality filtering settings:\n")
cat("Keeping only cells expressing at least this many genes:",minexpr,"\n")
cat("Keeping only genes expressed in at least this many cells:",mincells,"\n")
cat("Keeping only genes with at least this many reads across cells:",minreads,"\n")
cat("-----------------------\n")
cat("t-SNE settings:\n")
cat("Perplexity value(s):",opt$perplexity,"\n")
cat("# of PCA PCs to use:",opt$PCs,"\n")
cat("# of clusters k to look for:",opt$k,"\n")
cat("-----------------------\n")
if (length(opt$genelist) != 0) {
	cat("Use genes in this list to build trajectory:",opt$genelist,"\n")
	cat("-----------------------\n")
}
if (length(cell_info) != 0) {
	cat("Make plots of tSNE over factors in this file:",cell_info,"\n")
	cat("-----------------------\n")
}
cat("\n")

# FUNCTIONS
# ---------------------

# getOrder
# ---------------------
# get order of points resulting from a particular path through the graph 
# tsnematrix = df with 3 columns: V1 and V2 (tSNE coordinates) and cluster (cluster from kmeans above); rownames = cellnames
# path = df with 3 columns: V1 and V2 (coordinates of points along path, in tSNE space) and order (order of points along path); rownames not used
getOrder <- function(tsnematrix, path) {
	tomapx = cbind(tsnematrix$V1,tsnematrix$V2)
	rownames(tomapx) = rownames(tsnematrix)
	
	# for each point in the path, find it's cluster (which is the cluster of
	# the nearest point in the actual data)
	pathx = cbind(path$V1,path$V2)
	maptoclus = as.data.frame(nn2(tomapx, query = pathx, k = 1, searchtype = "standard", eps = 0)$nn.idx)	
	maptoclus$cluster = tsnematrix[maptoclus$V1,]$cluster
	path$cluster = maptoclus$cluster
	numclus = max(as.integer(maptoclus$cluster))
	
	# for each cluster, get all points in cluster and map to closest part of path (also in cluster)
	res = c()
	for (i in 1:numclus) {
		cellsubset = subset(tsnematrix, cluster == i)
		clusidx = which(path$cluster == i)

		# add 1 adjacent node from path outside of each end of cluster (some points may map to segment between two pathpoints in dif. clusters)
		startidx = max(1,(min(clusidx)-1))
		endidx = min(length(path$cluster),(max(clusidx)+1))
		pathsubset = path[startidx:endidx,]
		
		# get closest two points in path for every point in cellsubset
		cellxx = cbind(cellsubset$V1,cellsubset$V2)
		pathxx = cbind(pathsubset$V1,pathsubset$V2)
	
		restmp = as.data.frame(nn2(pathxx, query = cellxx, k = 2, searchtype = "standard", eps = 0)$nn.idx)
		rownames(restmp) = rownames(cellsubset)
		
		# map indexes back to original 'order' value in path
		restmp$newV1 = pathsubset$order[restmp$V1]
		restmp$newV2 = pathsubset$order[restmp$V2]
		restmp$V1 = NULL; restmp$V2 = NULL
		colnames(restmp) = c('V1','V2')

		# append to result
		res = rbind(res,restmp)
	}
	
	# res[rownames(res) == 'P18_5D',]
	
	# match back order in original tsnematrix
#COM	ord = data.frame(order = 1:length(rownames(tsnematrix)))
#COM	rownames(ord) = rownames(tsnematrix)
#COM	res = merge(res, ord, by=0, allx=TRUE)
#COM	res = res[order(res$order),]
#COM	res$order = NULL
#COM	rownames(res) = res$Row.names
#COM	res$Row.names = NULL

	# project all points onto their nearest path segments
	projections = c()
	segmentidx = c()

	for (i in 1:nrow(res)) {
		point = res[i,]
		lomatch = min(point$V1,point$V2)
		himatch = max(point$V1,point$V2)
	
		if (himatch - lomatch > 1) {
			# query's closest two points on the path are not adjacent; find segment with shortest
			# distance to this point by testing all segments between these two
#			cat(rownames(res)[i],', with i=',i,'has multiple segments that need to be resolved\n')
			mindist = Inf
			for (j in lomatch:(himatch-1)) {
				prj = nearestPointOnSegment(rbind(pathx[j,],pathx[j+1,]),tomapx[i,])
				if (prj[3] < mindist) {
					mindist = prj[3]
					projection = c(prj[1],prj[2])
					finalj = j
				}
			}	
			projections = rbind(projections,projection)
			rownames(projections)[i] = rownames(res)[i]
			segmentidx = c(segmentidx,finalj)
			names(segmentidx)[i] = rownames(res)[i]
		} else if (himatch - lomatch == 1) {
			# only one segment was closest
			prj = nearestPointOnSegment(rbind(pathx[point$V1,],pathx[point$V2,]),tomapx[i,])
			projection = c(prj[1],prj[2])
			projections = rbind(projections,projection)
			rownames(projections)[i] = rownames(res)[i]
			segmentidx = c(segmentidx,lomatch)
			names(segmentidx)[i] = rownames(res)[i]
		} else {
			stop("Internal error: could not map point to path. This should never happen! Something's gone really wrong.")
		}
	}

	# for every segment on the path in order, get all points that map to that
	# segment and get their internal order on the segment
	finalorder = c()
	for (i in 1:(nrow(pathx)-1)) {
#		cat('Projecting all cells nearest segment',i,'\n')
		cellist = names(segmentidx[segmentidx == i])
		prj = subset(projections, rownames(projections) %in% cellist)
		prj = rbind(pathx[i:(i+1),],prj)
		rownames(prj)[1:2] = c('start','end')
		sorted = prj[order(prj[,1], prj[,2]),]
		if (rownames(sorted)[1] == "end") {
			sorted = prj[order(prj[,1], prj[,2],decreasing=TRUE),]
		}
		sorder = rownames(subset(sorted, !(rownames(sorted) %in% c('start','end'))))
		finalorder = c(finalorder,sorder)
	}
	return(finalorder)
}


# isomap.incidence.matrix
# ---------------------
# function copied from this stackoverflow answer by user probaPerception:
# https://stackoverflow.com/questions/41598666/r-find-shortest-geodesic-path-between-2-points-of-a-2-d-point-cloud
isomap.incidence.matrix <- function (d, eps=NA, k=NA) {
  stopifnot(xor( is.na(eps), is.na(k) ))
  d <- as.matrix(d)
  if(!is.na(eps)) {
    im <- d <= eps
  } else {
    im <- apply(d,1,rank) <= k+1
    diag(im) <- FALSE
  }
  im | t(im)
}

# getPathLength
# ---------------------
# get length of path through a series of nodes, given a distance matrix between all nodes
# distm = NxN distance matrix where dist[n,m] = distance from n to m
# path = vector of length <= N with values <= M, representing a path through some or all of the N nodes in the dist matrix
getPathLength <- function(distm, path) {
	distance = 0
	for (i in 1:(length(path)-1)) {
		distance = distance + distm[path[i],path[i+1]]
	}
	return(distance)
}

# euclideanDist
# ---------------------
# basic euclidean distance between two points
euclideanDist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


# getDijkstra
# ---------------------
# note, a lot of this code is copied or adapted from this stackoverflow answer by user probaPerception:
# https://stackoverflow.com/questions/41598666/r-find-shortest-geodesic-path-between-2-points-of-a-2-d-point-cloud
# Inputs:
# tsnematrix = df with 3 columns: V1 and V2 (tSNE coordinates) and cluster (cluster from kmeans above); rownames = cellnames
# Output: 
# An ordered list of nodes corresponding to the shortest path through the centroid of each cluster
getDijkstra <- function(tsnematrix, kneighbors, exact) {

	# get number of clusters kk
	kk = max(as.integer(tsnematrix$cluster))
	if (kk == 1) {
		stop("Only one cluster found; this function currently doesn't support analysis when only one cluster exists")
	}

	# get centroid of each cluster (used later)
	centroidsx = vector(mode="double", length=kk)
	centroidsy = vector(mode="double", length=kk)
	for (i in 1:kk) {
		c1 <- subset(tsnematrix, cluster==i)
		centroidsx[i] = mean(c1$V1)
		centroidsy[i] = mean(c1$V2)
	}
	centroids = as.data.frame(cbind(centroidsx,centroidsy))

	# (1) get basic order of the clusters in the graph to try to traverse through
	# - basically, do mini djikstra's traversal through the clusters, and try to
	# put together physically adjacent clusters along the plot
	dists = matrix(, nrow = kk, ncol = kk)

	for (i in 1:kk) {
		for (j in i:kk) {
			c1 <- subset(tsnematrix, cluster==i)
			c2 <- subset(tsnematrix, cluster==j)
			dd = crossdist(c1$V1,c1$V2,c2$V1,c2$V2)
			minval = min(apply(dd,1,min))
			dists[i,j] = minval
			dists[j,i] = minval
		}
	}
	
	# find shortest path that goes through all clusters
	# note - this is the traveling salesman problem, only the brute-force approach
	# has been found. The answer is guaranteed correct, but super slow to get. To speed up, 
	# have the option of trying to identify a cluster that's likely on the endpoint, and 
	# only examine paths starting from that cluster. MUCH faster, but may fail.
	# Choose wisely. Note that if kk = 2 then there's no need (both clusters are 'ends').
	dists_graph = graph_from_adjacency_matrix(dists, mode = c("undirected"), weight = TRUE)
	maxpath = c()
	mindist = Inf

	if (exact == TRUE || kk <= 2) {	
		for (i in 1:kk) {
			# get all paths starting with node i
			sp = all_simple_paths(dists_graph,i)
			# reduce to those that traverse all kk nodes
			idx = which(lengths(sp)==kk)
		
			for (pp in 1:length(idx)) {
				dd = getPathLength(dists,sp[[idx[pp]]])
				if (dd < mindist) {
					maxpath = sp[[idx[pp]]]
					mindist = dd
				}
			}
		}
	} else {
		biggestdif = 0
		startfrom = NA
		for (i in 1:kk) {
			rr = sort(dists[,i])
			if (rr[3] - rr[2] > biggestdif) {
				biggestdif = rr[3] - rr[2]
				startfrom = i
			}
			if (is.na(startfrom)) {
				stop("Internal error: something's gone terribly wrong in getDijkstra(), couldn't find starting point. You may want to disable this part (which aims to greatly speed up the function) with the --exact option.")	
			}
		}
		
		sp = all_simple_paths(dists_graph,startfrom)
		# reduce to those that traverse all kk nodes
		idx = which(lengths(sp)==kk)
	
		for (pp in 1:length(idx)) {
			dd = getPathLength(dists,sp[[idx[pp]]])
			if (dd < mindist) {
				maxpath = sp[[idx[pp]]]
				mindist = dd
			}
		}
	}
	
	if (length(maxpath) == Inf) {
		stop("Internal error: could not find path through graph that included all clusters. Consider changing --smoother.")
	}

	# get endpoints of this path
	startclus = as.list(maxpath[1])[[1]]; endclus = as.list(maxpath[kk])[[1]]
	second = as.list(maxpath[2])[[1]]; secondlast = as.list(maxpath[kk-1])[[1]]

	# get path's starting node
	c1 <- subset(tsnematrix, cluster==startclus)
	c2 <- subset(tsnematrix, cluster==second)
	dd = crossdist(c1$V1,c1$V2,c2$V1,c2$V2)
	startidx = which(apply(dd,1,mean) == max(apply(dd,1,mean)))
	pathstart = rownames(c1[startidx,])
#	cat('Starting node is',pathstart,'with coordinates V1 =',c1[startidx,]$V1,'and V2=',c1[startidx,]$V2,'\n')
	startcoord = c(c1[startidx,]$V1, c1[startidx,]$V2)

	# get path's ending node
	c1 <- subset(tsnematrix, cluster==endclus)
	c2 <- subset(tsnematrix, cluster==secondlast)
	dd = crossdist(c1$V1,c1$V2,c2$V1,c2$V2)
	endidx = which(apply(dd,1,mean) == max(apply(dd,1,mean)))
	pathend = rownames(c1[endidx,])
#	cat('Ending node is',pathend,'with coordinates V1 =',c1[endidx,]$V1,'and V2=',c1[endidx,]$V2,'\n')
	endcoord = c(c1[endidx,]$V1, c1[endidx,]$V2)

	# start by getting path through centroids (result will be linked to start and end node above later)
	fullpath = c()

	for (i in 1:(kk-1)) {
		startnode = maxpath[i]
		endnode = maxpath[i+1]
#		cat('Using Dijkstra\'s to find shortest path between cluster',startnode,'and',endnode,'centroids\n')

		c1 = subset(tsnematrix, cluster %in% c(startnode,endnode))
		x = c(centroids$centroidsx[startnode],centroids$centroidsy[startnode])
		x = rbind(x,c(centroids$centroidsx[endnode],centroids$centroidsy[endnode]))
		rownames(x) = c('start','end')
		y = cbind(c1$V1,c1$V2)
		rownames(y) = rownames(c1)
		x = rbind(x,y)

		dist_matrix = as.matrix(dist(x))
		neighbors = isomap.incidence.matrix(dist_matrix, k=kneighbors)
		accessible = dist_matrix * neighbors
		access_graph <- graph_from_adjacency_matrix(accessible, mode = c("undirected"), weight = TRUE)
	
		sp = get.shortest.paths(access_graph, 1, to=2, mode = "all")$vpath[[1]]
		path = names(unlist(as.list(sp)))
		path = path[2:(length(path)-1)]
		
		# merge into fullpath
		if (any(path %in% fullpath)) {
			# at least one element in the new path already exists in the old; sometimes
			# backtracking from the centroid is a thing. As a cheap solution: if all the elements
			# shared are at the beginning of -path-, just drop them before appending
			shared = intersect(path,fullpath)
			idx = which(path %in% shared)
			maxidx = max(idx)
			if (identical(idx,1:maxidx)) {
				if (maxidx < length(path)) {
					fullpath = c(fullpath,path[(maxidx+1):length(path)])
				} else {
					cat('No path found between cluster ',startnode,' and ',endnode,' centers, skipping...\n')
				}
			} else {
				newend = match(path[max(idx)],fullpath)
				fullpath = c(fullpath[1:newend],path[(maxidx+1):length(path)])
			}
		
		} else {
			fullpath = c(fullpath,path)
		}
	}

	# finish path by connecting the start and end nodes from above
	
	# connect start node
	startnode = pathstart
	endnode = maxpath[1]
	c1 = subset(tsnematrix, cluster == endnode)
	ss = c1[rownames(c1) == startnode,]
	x = c(ss$V1,ss$V2)
	x = rbind(x,c(centroids$centroidsx[endnode],centroids$centroidsy[endnode]))
	rownames(x) = c(startnode,'end')
	y = cbind(c1[rownames(c1) != startnode,]$V1,c1[rownames(c1) != startnode,]$V2)
	rownames(y) = rownames(c1[rownames(c1) != startnode,])
	x = rbind(x,y)

	dist_matrix = as.matrix(dist(x))
	neighbors = isomap.incidence.matrix(dist_matrix, k=10)
	accessible = dist_matrix * neighbors
	access_graph <- graph_from_adjacency_matrix(accessible, mode = c("undirected"), weight = TRUE)

	sp = get.shortest.paths(access_graph, 1, to=2, mode = "all")$vpath[[1]]
	path = names(unlist(as.list(sp)))
	fullpath = c(path[1:(length(path)-1)],fullpath)

	# connect to last node
	startnode = maxpath[kk]
	endnode = pathend
	c1 = subset(tsnematrix, cluster == startnode)
	ss = c1[rownames(c1) == endnode,]
	x = c(centroids$centroidsx[startnode],centroids$centroidsy[startnode])
	x = rbind(x,c(ss$V1,ss$V2))
	rownames(x) = c('start',endnode)
	y = cbind(c1[rownames(c1) != endnode,]$V1,c1[rownames(c1) != endnode,]$V2)
	rownames(y) = rownames(c1[rownames(c1) != endnode,])
	x = rbind(x,y)

	dist_matrix = as.matrix(dist(x))
	neighbors = isomap.incidence.matrix(dist_matrix, k=10)
	accessible = dist_matrix * neighbors
	access_graph <- graph_from_adjacency_matrix(accessible, mode = c("undirected"), weight = TRUE)

	sp = get.shortest.paths(access_graph, 1, to=2, mode = "all")$vpath[[1]]
	path = names(unlist(as.list(sp)))
	fullpath = c(fullpath,path[2:length(path)])

	return(fullpath)
}

# runDijkstras
# ---------------------
# This is a wrapper for getDijkstra that does the error handling properly
runDijkstra <- function(tsnematrix, neighbors,exact) {
	out <- tryCatch(
		{
			fullpath = getDijkstra(tsnematrix,neighbors,exact)
			return(fullpath)
		},
		error=function(cond) {
			cat("Dijkstra's failed with the following inputs:\n")
			cat('tsnematrix:\n')
			print(tsnematrix)
			cat('neighbors:\n')
			print(neighbors)
			cat('exact:\n')
			print(exact)
			stop(cond)
		},
		warning=function(cond) {
			if (grepl('Couldn\'t reach some vertices', cond, fixed=TRUE)) {
				cat('Warning: could not find path through the points this run (couldn\'t reach all points).\n')
			} else {			
				cat("Dijkstra's produced the following warnings, proceeding anyway:\n")
			}
			return(NA)
		}
	)    
}


# sampleByCluster
# ---------------------
# draw a random sample, such that the same fraction of each cluster is sampled 
# tsnematrix = df with 3 columns: V1 and V2 (tSNE coordinates) and cluster (cluster from kmeans above); rownames = cellnames
# frac = fraction of points per cluster to sample
sampleByCluster <- function(tsnematrix, frac) {	
	kk = max(as.integer(tsnematrix$cluster))	
	ss = c()
	for (i in 1:kk) {
		c1 = subset(tsnematrix, cluster == i)
		nn = floor(nrow(c1)*frac)
		ss = rbind(ss, c1[sample(nrow(c1), nn), ])
	}
	return(ss)
}


# randomSample
# ---------------------
# draw a random sample from the tsnematrix
# tsnematrix = df with 2 columns: V1 and V2 (tSNE coordinates); rownames = cellnames
# frac = fraction of points per cluster to sample
randomSample <- function(tsnematrix, frac) {	
	nn = floor(nrow(tsnematrix)*frac)
	return(tsnematrix[sample(nrow(tsnematrix), nn), ])
}


# getPeakIdxs
# ---------------------
# finds all local maxima in the vector y, including at edges
getPeakIdxs <- function(y) {
	maxima = which(diff(sign(diff(y)))==-2)
	if (diff(head(y, 2)) < 0) {
		maxima = c(maxima,1)
	}
	if (diff(tail(y, 2)) > 0) {
		maxima = c(maxima,length(y))
	}
	return(sort(maxima))
}


# MAIN
# ---------------------

orig_matrix = read.table(fexpr_matrix, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)

if (all(sapply(orig_matrix, is.integer)) == FALSE) {
	stop("Error: input matrix contains non-integer values. Input matrix must have raw count data (not fpkm or other normalized data).")
}

if (length(cell_info) != 0) {
	cell_info = read.table(cell_info, row.names=1, sep='\t', header=TRUE, stringsAsFactors=FALSE)
}

# Filter matrix based on user inputs (orig_matrix -> expr_matrix)
cat("Filtering input matrix for low quality cells and very lowly expressed genes\n")
totgenes = nrow(orig_matrix)
totcells = ncol(orig_matrix)
cat("Original matrix contains",totgenes,"genes and",totcells,"cells\n")
filt_matrix = orig_matrix[, colSums(orig_matrix>0)>=minexpr]
if (ncol(filt_matrix) == 0) {
	stop("Oops, no cells remaining after filtering out cells expressing fewer than ",minexpr," genes")
}
totcells_filt = ncol(filt_matrix)
filt_matrix = filt_matrix[rowSums(filt_matrix)>=minreads, ]
if (ncol(filt_matrix) == 0) {
	stop("Oops, no genes remaining after filtering out all genes with fewer than ",minreads," reads detected across all cells")
}
filt_matrix = filt_matrix[rowSums(filt_matrix>0)>=mincells, ]
if (ncol(filt_matrix) == 0) {
	stop("Oops, no genes remaining after filtering out all genes expressed in fewer than ",mincells," cells")
}
totgenes_filt = nrow(filt_matrix)
cat("QC filtered matrix contains",totgenes_filt,"genes and",totcells_filt,"cells\n")


# Normalize reads if requested
if (CPM == TRUE) {
	totreads = colSums(filt_matrix)
	expr_matrix = t(t(filt_matrix) / totreads*1000000)
} else {
	expr_matrix = filt_matrix
}


# Take log2 if requested
if (takelog2 == TRUE) {
	expr_matrix = log2(expr_matrix + pseudocount)
	zeroval = log2(pseudocount)
} else {
	zeroval = 0
}


# Are we optimizing --perplexity and/or --PCs?
if (length(perplexity) > 1 || length(PCs) > 1 || length(karray) > 1) {
	optmode = TRUE
} else {
	optmode = FALSE
}


# Get set of genes to build trajectory over (expr_matrix -> expr_matrix_touse)
if (length(opt$genelist) != 0) {
	genelist = read.table(opt$genelist, header=TRUE, stringsAsFactors=FALSE,sep='\t')	
	expr_matrix_touse = subset(expr_matrix, rownames(expr_matrix) %in% genelist[,1])
	if (ncol(genelist) > 1) {
		genelist_factors = colnames(genelist)[2:length(colnames(genelist))]
	} else {
		genelist_factors = c()
	}
	cat("Subsetting matrix to only use genes listed in --genelist;",nrow(expr_matrix_touse),"rows remain.\n")
	if (nrow(expr_matrix_touse) != nrow(genelist)) {
		cat("WARNING: only", nrow(expr_matrix_touse), "out of", nrow(genelist),"genes in --genelist were found in the input matrix and passed QC filters - only those", nrow(expr_matrix_touse), "genes will be used in analysis.\n")
	}
} else {
	avgexpr = rowMeans(expr_matrix)
	variability = apply(expr_matrix,1,var)
	coefvar2 = variability/avgexpr^2
	genelist = names(sort(coefvar2,decreasing = TRUE)[1:topN])
	expr_matrix_touse = subset(expr_matrix, rownames(expr_matrix) %in% genelist)
	cat("Subsetting matrix to only use top",topN,"most variable genes\n")
	genelist_factors = c()
}

# Extra filtering after reducing dataset, to get rid of any columns (cells) with 
# all zeros
colsums = colSums(expr_matrix_touse)
expr_matrix_touse = expr_matrix_touse[, colSums(expr_matrix_touse)>zeroval]
cat("Subsetting matrix to drop any cells with no expression in selected genes;",ncol(expr_matrix_touse),"cells remain.\n")

if (length(cell_info) != 0) {
	tokeep = colnames(expr_matrix_touse)
	cell_info = cell_info[tokeep, ]
	cell_factors = c('cluster',colnames(cell_info))
} else {
	cell_factors = c('cluster')
}

# Run tSNE. If in optimization mode, run all combos here and exit.
for (perp in perplexity) {
	for (PC in PCs) {
		tsne = Rtsne(t(expr_matrix_touse), normalize = FALSE, pca_scale = TRUE, initial_dims = PC, dims = 2, perplexity=perp, verbose=TRUE, max_iter = 500)
		tsnematrix = tsne$Y
		rownames(tsnematrix) = colnames(expr_matrix_touse)

		# Use user-provided number if requested
		if (length(opt$k) != 0) {
			cat("Testing these k:",opt$k,"\n")
		} else {
			# learn optimal number of clusters k
			wss <- sapply(1:kmax, function(k){kmeans(tsnematrix, k, nstart=50,iter.max = 15 )$tot.withinss})
			wssd = diff(wss)

			# Use crude elbow method to ID # of clusters
			mwss = 5*max(wssd)
			wssd = wssd[wssd < mwss]
			k_opt = length(wssd)-1
			cat("Optimal k:",k_opt,"\n")
			if (k_opt > 10) {
				cat("Since k > 10, setting k = 10\n")
			}
			karray = c(k_opt)
		}
		
		for (kk in karray) {
			kkres = kmeans(tsnematrix, kk, nstart=50,iter.max = 150 )
			dftsnematrix = as.data.frame(tsnematrix)
			dftsnematrix$cluster = as.factor(kkres$cluster)
					
			# fit princurve or dijkstra's smoother != 'dijkstra' && smoother != 'smooth_spline' && smoother != 'lowess' && smoother != 'periodic_lowess'
			if (smoother == "dijkstra") {
				cat('Finding path through tSNE projection over',PC,'PCs with perplexity =',perp,', using k =',kk,', using the dijkstra approach\n')
				fullpath = runDijkstra(dftsnematrix,neighbors,exact)	# this is just a demo run, this is bootstrapped below
				if (! all(is.na(fullpath))) {
					# make plot with coloring by cluster + fitted line
					fullpathorder = data.frame("order" = 1:length(fullpath))
					rownames(fullpathorder) = fullpath
					pathpoints = merge(fullpathorder, dftsnematrix, by=0, allx=TRUE)
					pathpoints = pathpoints[order(pathpoints$order),]
					rownames(pathpoints) = pathpoints$Row.names
				}				
				
				if (length(cell_info) != 0) {
					dftsnematrix = merge(dftsnematrix,cell_info,by=0,merged=TRUE)
					rownames(dftsnematrix) = dftsnematrix$Row.names
					dftsnematrix$Row.names = NULL
				}

				# make plots	
				V1 = "V1"; V2 = "V2"
				for (cellfactor in cell_factors) {
					if (all(is.na(fullpath))) {
						a = ggplot(dftsnematrix)+ geom_point(aes_string(x=V1,y=V2,colour=cellfactor))+ theme_bw()
					} else {							
						a = ggplot(dftsnematrix)+ geom_point(aes_string(x=V1,y=V2,colour=cellfactor))+ theme_bw()
						a = a + geom_path(data = pathpoints, mapping = aes(x = V1, y = V2))
					}
					ff=paste(outprefix,"_tSNE_perp",perp,"_",PC,"PCs_k",kk,"_over_",cellfactor,".pdf",sep='')
					pdf(ff, width = 8.5, height = 8)
					print(a)
					graphics.off()				
					if (cellfactor != "cluster") {
						dftsnematrix[[cellfactor]] = NULL
					}
				}
								
			} else {
				cat('Finding path through tSNE projection over',PC,'PCs with perplexity =',perp,', using k =',kk,', by fitting a principal curve\n')
				# fit using princurve (default)
				# TODO - plotting over cell_factors not yet implemented for princurve
				x = cbind(dftsnematrix$V1,dftsnematrix$V2)
				lfit = principal_curve(x,smoother=smoother)
				lfitted = as.data.frame(lfit$s)

				# make plot with coloring by cluster + princurve
				a = ggplot(dftsnematrix)+ geom_point(aes(x=V1,y=V2,colour=cluster))+ theme_bw()
				a = a + geom_point(data = lfitted, mapping = aes(x = V1, y = V2))

				if (optmode == TRUE) {
					ff=paste(outprefix,"_tSNE_perp",perp,"_",PC,"PCs_k",kk,".pdf",sep='')
				} else {
					ff=paste(outprefix,"_tSNE_perp",perp,"_",PC,"PCs_k",kk,"_unsmoothed.pdf",sep='')
				}					
				pdf(ff, width = 8.5, height = 8)
				print(a)
				graphics.off()
			}	
		
			# output tSNE coordinates + clusters
			write.table(dftsnematrix, file=paste(outprefix,"_tSNE_perp",perp,"_",PC,"PCs_clusters.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)		
		}
	}
} 

if (optmode == TRUE) {
	# stop here if just performing parameter optimization
	cat("Done performing tSNE analysis with provided --perplexity and --PCs parameters. Please inspect the results and re-run with a single value per parameter to do full analysis.\n")
	quit(status=0, save='no')
}

# if not in optmode, continue. 
if (smoother == 'dijkstra') {

	pathlist = c()

	# smooth the fitted curve by running many times on a subset of points
	numiter = 200
	n = 1
	while (n <= numiter) {
#		cat('Iteration',n,'\n')
		tsnesubset = sampleByCluster(dftsnematrix, fracdraw)
		tsnesubset$cluster = NULL
		
		kkres = kmeans(tsnesubset, kk, nstart=50,iter.max = 150)
		tsnesubset$cluster = as.factor(kkres$cluster)
	
		fullpath = runDijkstra(tsnesubset,neighbors,exact)
		if (length(fullpath) == 1 && is.na(fullpath[1])) {
			cat('Things are fine if this error only appears sporadically, but if it\'s happening a lot, try increasing --fracdraw and/or --neighbors\n')
			cat('Repeating this iteration (#',n,')\n')
		} else {
			fullpathorder = data.frame("order" = 1:length(fullpath))
			rownames(fullpathorder) = fullpath
			pathpoints = merge(fullpathorder, tsnesubset, by=0, allx=TRUE)
			pathpoints = pathpoints[order(pathpoints$order),]
			rownames(pathpoints) = pathpoints$Row.names
			pathpoints$Row.names = NULL
		
			# save trajectory for this run
			pathpoints$cluster = NULL
			pathlist[[paste('list',n,sep='')]] = pathpoints

			n = n + 1
		}
	}
	
	# get consensus 'smooth' trajectory based on all the numiter individual trajectories
	# by using bipartite matching
	listA = pathlist$list1
	rownames(listA) = paste('A',1:length(listA$V1),sep='')
	consensus = listA		# build consensus starting with the very first list

	for (i in 2:numiter) {
#		cat('Iteration',i,'\n')
		listA = consensus
		listB = pathlist[[paste('list',i,sep='')]]	
		rownames(listB) = paste('B',1:length(listB$V1),sep='')
		
		dd = crossdist(listA$V1,listA$V2,listB$V1,listB$V2)
		ddinv = -1*dd + max(dd)		# maximum.bipartite.matching only maximizes, whereas we want to min, so transform so lowest dist becomes highest value
	
		bpg <- graph.full.bipartite(length(rownames(listA)),length(rownames(listB)))
		V(bpg)$name <- c(rownames(listA), rownames(listB))
		E(bpg)$weight = t(ddinv)
	
		matched = maximum.bipartite.matching(bpg)$matching
	
		consensus = c()
		for (i in 1:length(matched)) {
			if (!is.na(names(matched[i])) && !is.na(matched[i]) && grepl("B", matched[i], fixed=TRUE)) {
				a = names(matched[i])
				b = matched[i]		
				pp = cbind(listA[rownames(listA)==a,]$V1,listA[rownames(listA)==a,]$V2)
				pp = rbind(pp,cbind(listB[rownames(listB)==b,]$V1,listB[rownames(listB)==b,]$V2))		
				consensus = rbind(consensus,colMeans(pp))
			}
		}	
	
		consensus = as.data.frame(consensus)
		rownames(consensus) = paste('A',1:length(consensus$V1),sep='')	
	}
	
	# consensus may be slightly out of order due to small changes during iteration; fix it
	res = nn2(consensus, k = nrow(consensus), searchtype = "standard", eps = 0)$nn.idx

	# for each point, find it's next closest -that has not yet been incorporated into the chain-
	finalconsensus = consensus[1,]
	included = c(1)
	i = 1
	while (i < length(res[,1])) {
		rr = res[i,]
#		cat('Looking for neighbor to node',i,'\n')
		found = FALSE
		for (j in 1:nrow(consensus)) {
			if (!(rr[j] %in% included)) {
#				cat('Adding edge',i,'->',rr[j],'\n')
				finalconsensus = rbind(finalconsensus, consensus[rr[j],])
				included = c(included, rr[j])
				i = rr[j]
				found = TRUE
				break
			}
		}
		if (found == FALSE) {
			stop("Internal error while smoothing trajectory - couldn't find smooth path")
		}
	}
	
	# Now, get consensus ordering of points along the final consensus path.
	finalconsensus$order = c(1:length(rownames(finalconsensus)))
	res = getOrder(dftsnematrix, finalconsensus)

	ordering = data.frame(trajectory = c(1:length(res)))
	rownames(ordering) = res
	dftsnematrix = merge(dftsnematrix, ordering, by=0, all=TRUE)
	rownames(dftsnematrix) = dftsnematrix$Row.names
	dftsnematrix$Row.names = NULL
	
	a = ggplot(dftsnematrix)+ geom_point(aes(x=V1,y=V2,colour=trajectory))+ theme_bw() + scale_color_distiller(palette = "Spectral")
	a = a + geom_path(data = finalconsensus, mapping = aes(x = V1, y = V2), size = 1.5)		
	ff=paste(outprefix,"_tSNE_perp",perp,"_",PC,"PCs_k",kk,"_final_trajectory.pdf",sep='')
	pdf(ff, width = 8.5, height = 8)
	print(a)
	graphics.off()

	a = ggplot(dftsnematrix)+ geom_point(aes(x=V1,y=V2,colour=trajectory))+ theme_bw() + scale_color_distiller(palette = "Spectral")
	ff=paste(outprefix,"_tSNE_perp",perp,"_",PC,"PCs_k",kk,"_final_trajectory_noline.pdf",sep='')
	pdf(ff, width = 8.5, height = 8)
	print(a)
	graphics.off()
		
	# output final consensus trajectory
	write.table(dftsnematrix, file=paste(outprefix,"_trajectory.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)		
	write.table(finalconsensus, file=paste(outprefix,"_trajectory_coordinates.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)		

} else {
	stop('Additional analyses using princurve not yet implemented, please use --smoother dijkstra for full analysis with bootstrapping.')
}

# OK, now that we have the the cell ordering along the trajectory, look at gene expression along
# the trajectory

# Heat map + clustering of genes in the genelist subset to get a sense of results
#hmcolorscale <- colorRampPalette(c("blue","yellow"), space = "rgb")(30)		# color palette for all heatmaps
hmcolorscale = viridis(30)

tmp = t(expr_matrix_touse)
expr_time = merge(dftsnematrix, tmp, by=0, all=TRUE)
rownames(expr_time) = expr_time$Row.names
expr_time$Row.names = NULL
expr_time <- t(expr_time[order(expr_time$trajectory),5:length(colnames(expr_time))])

colorder = 1:ncol(expr_time)
colorschemey = colorRampPalette(brewer.pal(11, "Spectral"))(length(colorder))

if (length(genelist_factors) > 0) {
	genelist$order = 1:length(rownames(genelist))
	tmp = merge(data.frame(genelist[,1],genelist$order),expr_time, by.x=1,by.y=0,merged=TRUE)
	rownames(tmp) = tmp[,1]
	tmp[,1] = NULL
	tmp <- tmp[order(tmp[,1]),] 
	tmp[,1] = NULL

	for (fact in genelist_factors) {
		forder = unique(genelist[[fact]])
		ff = factor(genelist[[fact]],levels = forder)
		
		colorschemex = colorRampPalette(brewer.pal(11, "Spectral"))(nlevels(ff))
		names(colorschemex) <- levels(ff)
			
		pdf(paste(outprefix,"_heatmap_over_gene_subset_w_",fact,".pdf",sep=''), width = 8.5, height = 8)
		heatmap.2(as.matrix(tmp), Rowv=FALSE, Colv=FALSE, col = hmcolorscale, density.info="none", trace="none",dendrogram='none', RowSideColors=colorschemex[ff], ColSideColors=colorschemey[colorder]) 
		legend(x="bottomleft", legend=forder, col=colorschemex[forder], pch=15)		
		graphics.off()
	}
} 

hr <- hclust(as.dist(1-cor(t(expr_time), method="pearson")), method="complete")

pdf(paste(outprefix,"_heatmap_over_gene_subset.pdf",sep=''), width = 8.5, height = 8)
heatmap.2(as.matrix(expr_time), Rowv=as.dendrogram(hr), Colv=FALSE, col = hmcolorscale, density.info="none", trace="none",dendrogram='row', ColSideColors=colorschemey[colorder]) 
graphics.off()

if (testing == TRUE) {
	# stop here if just testing trajectory params
	cat("Done performing tSNE analysis and building trajectory. Please inspect the results and re-run without the --testing option to do the full analysis.\n")
	quit(status=0, save='no')
}


# Find genes whose expression changes along the trajectory axis
# Regression approach

# need to go back to raw counts for this part (-zeroinfl- models count data, not normalized)
# first, get rid of any columns (cells) in original filtered matrix not in this final
# matrix (occurs if --genelist specified and some cells had no data in that genelist)
cat('Identifying genes whose expression varies significantly along the trajectory\n')
filt_matrix = filt_matrix[, colnames(filt_matrix) %in% rownames(dftsnematrix)]

# use edgeR to convert raw counts to normalized 'pseudocounts' (not to be confused with
# the --pseudocount parameter to this script, which is a number added before log-transforming
# to prevent taking log of 0) -> using procedure outlined here: https://rdrr.io/bioc/tweeDEseq/src/R/normalizeCounts.R
d = DGEList(counts=filt_matrix)
d = calcNormFactors(d)
d = estimateCommonDisp(d)
lib.size = d$samples$lib.size * d$samples$norm.factors
d2 = edgeR::DGEList(counts=filt_matrix, lib.size=lib.size)
dispersion = rep(d$common.dispersion, nrow(filt_matrix))
dispersion <- pmax(dispersion, 1e-06)
pcs = ceiling(edgeR::equalizeLibSizes(d2, disp=dispersion)$pseudo.counts-0.5)

for (i in 1:10) {
	for (j in 1:10) {
		if (filt_matrix[j,i] == 0) {
			pcs[j,i] = 0
		}
	}
}
	
# format data for analysis
expr_time = t(pcs)
totreads = rowSums(expr_time)
expr_time = merge(dftsnematrix, expr_time, by=0, all=TRUE)
rownames(expr_time) = expr_time$Row.names
expr_time$Row.names = NULL
psdt <- expr_time$trajectory

genes = colnames(expr_time)[5:length(colnames(expr_time))]
#genes = genes[1:5000]
#genes = genelist[,1]

minobs = 10					# hurdle doesn't do well with very few obs.
zeroinflcutoff = 0.01		# hurdle doesn't do well with very few zeros, either
noeval = 0
linear = 0
poly2 = 0
regresults = c()

n = 1
for (gene in genes) {
	if (n %% 100 == 0) {
		cat('Testing ',n,'th gene for correlation with trajectory...\n',sep='')
	}

	if (sum(expr_time[[gene]]>0) >= minobs) {
		if (sum(expr_time[[gene]]==0) / length(expr_time[[gene]]==0) <= zeroinflcutoff) {
			# corner case: very few zeros in cell. Run an ordinary negbin regression instead since data aren't zero-inflated.
			reg1 <- glm.nb(expr_time[[gene]] ~ poly(psdt,2))
			reg2 <- glm.nb(expr_time[[gene]] ~ psdt)

			if (summary(reg1)$aic < summary(reg2)$aic) {
				regres <- c(coef(summary(reg1))[2,4],coef(summary(reg1))[3,4])
				linear = linear + 1
				testedmodel = 'negbin_linear'
			} else {
				regres <- c(coef(summary(reg2))[2,4])
				poly2 = poly2 + 1
				testedmodel = 'negbin_poly2'
			}
		} else {
			tryCatch(
				{
					reg1 <- hurdle(expr_time[[gene]] ~ psdt, dist = 'negbin')
				},
				error = function(e) { 
					reg1 <<- c()
				}
			)	
			tryCatch(
				{
					reg2 <- hurdle(expr_time[[gene]] ~ poly(psdt,2), dist = 'negbin')
				},
				error = function(e) { 
					reg2 <<- c()
				}
			)	
			
			if (length(reg1) == 0 && length(reg2) == 0) {
				# both regressions failed; gene could not be modeled
				noeval = noeval + 1
				testedmodel = 'not_tested'
				pval = NA
			} else if (length(reg1) == 0) {
				# only the poly model could be tested, use it
				regres <- c(coef(summary(reg2))$count[2,4],coef(summary(reg2))$count[3,4],coef(summary(reg2))$zero[2,4],coef(summary(reg2))$zero[3,4])
				poly2 = poly2 + 1
				testedmodel = paste('hurdle_poly2',sep='')
			} else if (length(reg2) == 0) {
				# only the linear model could be tested, use it
				regres <- c(coef(summary(reg1))$count[2,4],coef(summary(reg1))$zero[2,4])
				linear = linear + 1
				testedmodel = paste('hurdle_linear',sep='')
			} else {
				# both models worked; get pvals from the one with best loglikelihood			
				if (summary(reg1)$loglik > summary(reg2)$loglik) {
					# linear model has better loglikelihood despite fewer parameters/dfs, use it
					regres <- c(coef(summary(reg1))$count[2,4],coef(summary(reg1))$zero[2,4])
					linear = linear + 1
					testedmodel = paste('hurdle_linear',sep='')
				} else {
					regres <- c(coef(summary(reg2))$count[2,4],coef(summary(reg2))$count[3,4],coef(summary(reg2))$zero[2,4],coef(summary(reg2))$zero[3,4])
					poly2 = poly2 + 1
					testedmodel = paste('hurdle_poly2',sep='')
				}
			}
		}
		
		if (length(regres) == 0 || all(is.na(regres))) {
			stop("Internal error: hurdle regression for",gene,"failed; please try alternate models")
		}
	
		pval = min(regres,na.rm=TRUE)
	} else {
		noeval = noeval + 1
		testedmodel = 'not_tested'
		pval = NA
	}	
	
	regresults = rbind(regresults,c(gene,testedmodel,pval))
	n = n+1
}

regresultsdf = data.frame(row.names=regresults[,1],model=regresults[,2],pval=as.double(regresults[,3]), stringsAsFactors = FALSE)
regresultsdf$padj = p.adjust(regresultsdf$pval, method='bonferroni')
sig_genes = rownames(regresultsdf[Vectorize(isTRUE)(regresultsdf$padj < sigthreshold),])

# output final consensus trajectory
write.table(regresultsdf, file=paste(outprefix,"_stat_tests.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)		
write.table(sig_genes, file=paste(outprefix,"_significant_genes.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)		
a = table(regresultsdf$model)
cat('Out of',nrow(regresultsdf)-a[names(a) == 'not_tested'],'genes tested,',length(sig_genes),'varied significantly along the trajectory\n')

# useful debugging plots for looking at expression patterns of single gene:
# totreads = colSums(expr_time)
# plot(factor(expr_time[[gene]] == 0) ~ expr_time$trajectory, main = "Zero")
# plot(log(expr_time[[gene]]/totreads*1000000) ~ expr_time$trajectory, subset = expr_time[[gene]] > 0, main = "Count")
# plot(expr_time[[gene]] ~ expr_time$trajectory, subset = expr_time[[gene]] > 0, main = "Count")
# hist(log(expr_time[[gene]]/totreads*1000000))
# hist(expr_time[expr_time[gene]>0,gene])

# get the significant genes
cat('Making plots of expression of significant genes along trajectory\n')
toplot = sig_genes

# make heatmaps of significant genes
tmp = t(expr_matrix)
tmp = tmp[rownames(tmp) %in% rownames(dftsnematrix),]
expr_matrix_plt = merge(dftsnematrix, tmp, by=0, all=TRUE)
rownames(expr_matrix_plt) = expr_matrix_plt$Row.names
expr_matrix_plt$Row.names = NULL
expr_matrix_plt <- t(expr_matrix_plt[order(expr_matrix_plt$trajectory),5:length(colnames(expr_matrix_plt))])

# histogram with genes clustered by hierarchical clustering
expr_matrix_subset = expr_matrix_plt[rownames(expr_matrix_plt) %in% toplot,]
hr <- hclust(as.dist(1-cor(t(expr_matrix_subset), method="spearman")), method="single")

colorder = 1:ncol(expr_matrix_subset)
colorschemey = colorRampPalette(brewer.pal(11, "Spectral"))(length(colorder))

pdf(paste(outprefix,"_heatmap_sig_genes_hierarchical.pdf",sep=''), width = 8.5, height = 8)
heatmap.2(as.matrix(expr_matrix_subset), Rowv=as.dendrogram(hr), Colv=FALSE, col = hmcolorscale, density.info="none", trace="none",dendrogram='row', ColSideColors=colorschemey[colorder]) 
graphics.off()
             

# kmeans clustering of the genes instead of hierarchical
kkres = kmeans(expr_matrix_subset, kgenes, nstart=100,iter.max = 100)

roworder = as.data.frame(kkres$cluster)
tmp = merge(roworder,expr_matrix_subset,by=0,merged=TRUE)
rownames(tmp) = tmp[,1]
tmp[,1] = NULL
tmp <- tmp[order(tmp[,1]),] 

forder = unique(tmp[,1])
ff = factor(tmp[,1],levels = forder)
tmp[,1] = NULL
colorschemex = colorRampPalette(brewer.pal(11, "Spectral"))(nlevels(ff))
names(colorschemex) <- levels(ff)

#png(paste(outprefix,"_heatmap_sig_genes.png",sep=''), width = 8.5, height = 8, units = 'in', res = 300)
#heatmap.2(as.matrix(tmp), Rowv=FALSE, Colv=FALSE, col = hmcolorscale, density.info="none", trace="none",dendrogram='none', RowSideColors=colorschemex[ff], ColSideColors=colorschemey[colorder]) 
#legend(x="bottomleft", legend=forder, col=colorschemex[forder], pch=15)		
#graphics.off()

# for each k-means cluster, get average signal in every cell and build average expression
# profile of the cluster
tmp2 = merge(roworder,expr_matrix_subset,by=0,merged=TRUE)
rownames(tmp2) = tmp2[,1]
tmp2[,1] = NULL

modes = c()
for (k in 1:kgenes) {
	subsetk = tmp2[tmp2[,1] == k,2:ncol(tmp2)]
	avgpercell = colMeans(subsetk)
	profilek = data.frame(cell=1:ncol(subsetk),expr=avgpercell,stringsAsFactors=FALSE)
	
	spl = smooth.spline(profilek$cell, profilek$expr,df=floor(length(avgpercell)*4/1000)+1)
	maxima = getPeakIdxs(spl$y)
	maxorder = order(-spl$y[maxima])
	
	if (length(maxima) == 1) {
		modes = rbind(modes,c(k,maxima[1],0,mean(profilek$expr)))	
	} else {
		modes = rbind(modes,c(k,maxima[maxorder[1]],maxima[maxorder[2]],mean(profilek$expr)))	
	}		
#	plot(profilek$cell, profilek$expr)
#	lines(spl$x,spl$y)
}

# further flatten into quintiles
modes[,2:3] = floor(modes[,2:3]/length(avgpercell)*5)
newkorder = modes[order(modes[,2], modes[,3], modes[ ,4],decreasing = FALSE),1]
newkorder = data.frame(neworder = 1:kgenes, oldorder = newkorder)

# plot the k-means clusters in this new order
colnames(roworder)[1] = 'oldorder'
newkorder <- join(roworder, newkorder, by='oldorder', type='left', match='all')
rownames(newkorder) = rownames(roworder)
newkorder$oldorder = NULL

tmp = merge(newkorder,expr_matrix_subset,by=0,merged=TRUE)
rownames(tmp) = tmp[,1]
tmp[,1] = NULL
tmp <- tmp[order(tmp[,1]),] 

forder = unique(tmp[,1])
ff = factor(tmp[,1],levels = forder)
tmp[,1] = NULL
colorschemex = colorRampPalette(brewer.pal(11, "Spectral"))(nlevels(ff))
names(colorschemex) <- levels(ff)

pdf(paste(outprefix,"_heatmap_sig_genes.pdf",sep=''), width = 8.5, height = 8)
heatmap.2(as.matrix(tmp), Rowv=FALSE, Colv=FALSE, col = hmcolorscale, density.info="none", trace="none",dendrogram='none', RowSideColors=colorschemex[ff], ColSideColors=colorschemey[colorder]) 
legend(x="bottomleft", legend=forder, col=colorschemex[forder], pch=15)		
graphics.off()


# finally, make plots of average expression profiles along trajectory for each
# cluster, plus their spline fit
tmp2 = merge(newkorder,tmp,by=0,merged=TRUE)
rownames(tmp2) = tmp2[,1]
tmp2[,1] = NULL

for (k in 1:kgenes) {
	subsetk = tmp2[tmp2[,1] == k,2:ncol(tmp2)]
	avgpercell = colMeans(subsetk)
	profilek = data.frame(cell=1:ncol(subsetk),expr=avgpercell,stringsAsFactors=FALSE)
	
	spl = smooth.spline(profilek$cell, profilek$expr, df=floor(length(avgpercell)*4/1000)+1)
	spldf = data.frame(V1=spl$x,V2=spl$y)

	n = floor(length(avgpercell)*4/1000)+1
	runsum = aggregate(profilek,list(rep(1:(nrow(profilek)%/%n+1),each=n,len=nrow(profilek))),mean)[-1]

	pdf(paste(outprefix,"_genes_k",k,".pdf",sep=''), width = 8.5, height = 4)
		a = ggplot(runsum)+ geom_point(aes(x=cell,y=expr,colour=cell),shape=1)+ theme_bw() + scale_color_distiller(palette = "Spectral")
		a = a + geom_path(data = spldf, mapping = aes(x = V1, y = V2), size = 1.5)	
		a = a + xlab('trajectory')
		a = a + ylab('average log2(CPM)')
		print(a)	
	graphics.off()
}

colnames(newkorder)[1] = 'cluster'
write.table(newkorder, file=paste(outprefix,"_significant_genes_wcluster.txt",sep=''), sep='\t', quote=FALSE, col.names=NA)		


# other useful plots:
# cc <- viridis(10)[as.numeric(cut(log2(expr_time[[gene]]+1),breaks = 10))]
# plot(cbind(expr_time$V1, expr_time$V2), col=cc)











