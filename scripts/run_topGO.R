#!/usr/bin/env Rscript

# version 1.0 (09/07/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 05/27/2020
# -------------------------

# Description:
# -------------------------
# This is a simple wrapping script for running GO-enrichement analyses using the R package 'topGO',
# following the guide here:
# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf
# and this BioStars post with more plant-specific advice:
# https://www.biostars.org/p/250927/

# List of required libraries:
# -------------------------
liblist = c('argparse','topGO','biomaRt','Rgraphviz','grid','ggplot2','pheatmap')
	
# Notes:
# -------------------------
# (1) If you want to use a different database than the default one (which is set for A. thaliana), here's
# how to find the correct settings for --GOmart, --GOdataset, and --GOhost:
# Note that hosts, marts and datasets are nested, so a host will host multiple marts, each of which will
# contain multiple datasets.

# NOTE - if your organism is an animal, you probably don't need to do any of this - just see note (2)!

# Open R and load the package 'biomaRt':
# > library("biomaRt")

# To list the marts available for a particular host, do (here using plants.ensembl.org as an example):
# > listMarts(host="plants.ensembl.org")

# To list all the datasets available for a specific mart, do (here using plants_mart within plants.ensembl.org as an example):
# > m <- useMart("plants_mart", host="plants.ensembl.org")
# > listDatasets(m)

# Your dataset must have attributes 'ensembl_gene_id' and 'go_id', which you can check for a specific dataset (here using athaliana_eg_gene within plants_mart as example):
# > m <- useMart("plants_mart", dataset="athaliana_eg_gene", host="plants.ensembl.org")
# > listAttributes(ensembl)[listAttributes(ensembl)$name == 'ensembl_gene_id',]
# > listAttributes(ensembl)[listAttributes(ensembl)$name == 'go_id',]

# If the output of the last two lines is -not- an empty dataframe, you should be all set.

# (2) If your organism is an animal, it's probably in the general ensembl mart,
# which you can check in R by doing:
# > library("biomaRt")
# > listDatasets(ensembl)

# If your dataset is in that list, call this script with that dataset for --GOdataset, and set the --ensembl flag option.
# This will ignore --GOhost and --GOmart and use the ensembl mart. So for example, for M. musculus, do:
# run_topGO.R mygenelist.txt outprefix --GOdataset mmusculus_gene_ensembl --ensembl


# -------------------------
# Usage: run_topGO.R [options] genelist.txt outprefix
# -------------------------
# Check that all required libraries are installed and load them
if(length(setdiff(liblist, rownames(installed.packages()))) > 0) {
	missinglibs = setdiff(liblist, rownames(installed.packages()))
	missinglibs = paste(as.character(missinglibs), sep="' '", collapse=", ")
	stop("the following required libraries are not installed: ",missinglibs,". Please install them, then re-run this script.")
}
cat("Loading all required libraries...")
for (lib in liblist) {
	suppressMessages(library(lib, character.only = TRUE))
}
cat("DONE\n")
	
# Read in user-supplied arguments
parser = ArgumentParser()

# Required arguments
parser$add_argument("genelist", nargs=1, help = "NO HEADER! A file with a list of genes for GO analysis; if this file has two columns, second column assumed to provide groups of genes for different samples/conditions, with the value of column 2 indicating the name of the sample/condition.", type="character")
parser$add_argument("outprefix", nargs=1, help = "prefix for output file(s)")

# Other arguments
parser$add_argument("--background", help = "NO HEADER! List of genes to use as background dataset (else uses all genes downloaded with biomart as ref)", type="character")
parser$add_argument("--samplename", default = "sample", help = "name of single sample (if not provided as second column in genelist), if --expr_data provided should match the column name, or if --samp_map provided should match second column of that file", type="character")
parser$add_argument("--sampleorder", help = "(only used if >1 sample provided) in heatmaps, plot samples in this order (default orders them using hierarchical clustering) - note that only the samples in --sampleorder will be used, even if more samples are present in the data matrix/samp_map!", type="character")
parser$add_argument("--GOmart", default = "plants_mart", help = "BioMart database to use, default plant database", type="character")
parser$add_argument("--GOdataset", default = "athaliana_eg_gene", help = "BioMart dataset to use, default A. thaliana dataset", type="character")
parser$add_argument("--GOhost", default = "plants.ensembl.org", help = "Host to use for downloading GO annotations with bioMart", type="character")
parser$add_argument("--pval", default = 0.01, help = "Host to use for downloading GO annotations with bioMart", type="double")
parser$add_argument("--algorithm", default = "elim", help = "topGO algorithm to use (see topGO manual)", type="character")
parser$add_argument("--statistic", default = "fisher", help = "topGO test statistic to use (see topGO manual)", type="character")
parser$add_argument("--expr_data", help = "YES HEADER! matrix of expression data, genes (rows) x samples (columns), with column 1 = gene IDs and row 1 = sample IDs", type="character")
parser$add_argument("--samp_map", help = "NO HEADER! mapping of sample IDs in genelist to sample IDs in expr_data (if not same); can have multiple sample IDs in expr_data assigned to same sample_id in genelist (e.g. replicates); these will be averaged for final plots", type="character")
parser$add_argument("--pvalmax", help = "Maximum value for color scale in -log10(p-values) heatmap or axis in barchart", type="double")
parser$add_argument("--fillupper", help = "Upper limit of value for color scale for heatmap and/or dotplot (default automatically set to highest value in plot)", type="double")
parser$add_argument("--filllower", help = "Lower limit of value for color scale for heatmap and/or dotplot (default automatically set to lowest value in plot)", type="double")
parser$add_argument("--sizeupper", help = "Upper limit of value for size scale for dotplot (default automatically set to highest value in plot)", type="double")
parser$add_argument("--dotsize", default = 15, help = "Upper limit for size of points in dot plot (default 10)", type="double")
parser$add_argument("--topN", help = "When a large number of terms is identified across all samples, an additional plot of the 5 top terms per sample only will also be created", type="integer")
parser$add_argument("--ensembl", default=FALSE, action="store_true", help = "use one of the ensembl datasets for downloading GO-terms (see notes above)")
parser$add_argument("--single_cell", default=FALSE, action="store_true", help = "(only used if --expr_data provided); activates single-cell mode for the expression plots, so that dot plots (color = expr, size = fraction nonzero) is used instead of heatmap to plot avg. expression over significant GO-terms")
parser$add_argument("--GEOmean", default=FALSE, action="store_true", help = "(only used if --expr_data provided); take log2 of expression values before calculating mean (result is similar to a geometric mean instead of arithmetic mean); a pseudocount of 1 is added to each value")
parser$add_argument("--includezeros", default=FALSE, action="store_true", help = "(only relevant if --single_cell is on, otherwise this is always considered TRUE) - include zero values when calculating mean expression (normally zeros are excluded from mean calculation, represented instead by frac nonmissing = dot size))")
parser$add_argument("--useallGOtermgenes", default=FALSE, action="store_true", help = "(only used for expression plot) - when calculating avg. expression over each GO term, use all genes that fall under GO term, including those that were not provided as part of the list of significant terms (genelist) - default only uses the GO-term genes that are also in genelist")

opt <- parser$parse_args()

zz <- file(paste(opt$outprefix,"_ext_logfile.txt",sep=''), open = "wt")

# FUNCTIONS
# ---------------------
runtopGO <- function(genes, golist, algorithm = 'elim', statistic = 'fisher') {
	sink(zz, type = "message")
	GOdata = new("topGOdata", ontology='BP', allGenes = genelbl, annot = annFUN.gene2GO, gene2GO = golist)
	results = runTest(GOdata, algorithm = algorithm, statistic = statistic)
	sink(type = "message")
	return(list(GOdata,results))
}

wide_to_long <- function(mat) {
	for (cc in 1:ncol(mat)) {
		tmp = mat[,cc]
		tmp = data.frame(samp=rep(colnames(mat)[cc],length(tmp)), val=tmp)
		tmp$GO.ID = rownames(mat)
		if (cc != 1) {
			res = rbind(res,tmp)
		} else {
			res = tmp
		}
	}
	return(res)
}

# MAIN
# ---------------------

# Check all inputs ok
if (is.null(opt$genelist)) stop("Must provide a list of genes as first argument")
if (is.null(opt$outprefix)) stop("Must provide a prefix for output files as second argument")
if (opt$single_cell == FALSE) opt$includezeros = TRUE

# Read in genelist
cat("Reading in input files...")
genelist = read.table(opt$genelist, header=FALSE, sep="\t", stringsAsFactors = FALSE)
samplesprovided = TRUE
if (ncol(genelist) == 1) {
	# only one gene group provided, check all ok
	genelist$V2 = opt$samplename
	samplesprovided = FALSE
} else if (ncol(genelist) != 2) {
	stop("Input 'genelist' file must have either one or two columns")
} else if (ncol(genelist) == 2) {
	# second column exists in file, but has only one value
	if (length(unique(genelist$V2)) == 1) {
		samplesprovided = FALSE
		cat("Warning: genelist has two columns, but all values of second column are the same, so treating this as a single sample\n")
	}
}

# Read in sample order, if provided
if (! is.null(opt$sampleorder)) {
	colorder = unlist(strsplit(opt$sampleorder,","))
	if (length(unique(colorder)) != length(colorder)) stop("\nCannot list same sample more than once in --sampleorder")
	if (samplesprovided == TRUE) {
		for (ss in colorder) {
			if (! (ss %in% unique(genelist$V2))) stop("\nSample '",ss,"' (from --sampleorder) was not found in second column of genelist",sep='')
		}
	}
	genelist = genelist[genelist$V2 %in% colorder,]
}

# Read in expression data, if provided, check all ok
if (! is.null(opt$expr_data)) {
	expr_data = data.matrix(read.table(opt$expr_data, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
	if (! is.null(opt$samp_map)) {
		samp_map = read.table(opt$samp_map, header=FALSE, sep="\t", row.names=1, stringsAsFactors = FALSE)		
		if (samplesprovided == TRUE) {
			for (ss in unique(genelist$V2)) {
				if (! (ss %in% samp_map$V2)) stop("\nSample '",ss,"' (from genelist) was not found in second column of provided --samp_map",sep='')
			}
			samp_map = samp_map[samp_map$V2 %in% unique(genelist$V2),,drop=FALSE]	
		} else {
			if (!(opt$samplename %in% samp_map$V2)) stop("\nSample '",opt$samplename,"' not found in second column of --samp_map",sep='')
			samp_map = samp_map[samp_map$V2 == opt$samplename,,drop=FALSE]
		}
		for (ss in rownames(samp_map)) {
			if (! (ss %in% colnames(expr_data))) stop("\nColumn '",ss,"' from first column of --samp_map not found in --expr_matrix",sep='')
		}
	} else {
		if (samplesprovided == TRUE) {
			for (ss in unique(genelist$V2)) {
				if (! (ss %in% colnames(expr_data))) stop("\nNo matching column for sample '",ss,"' in --expr_data",sep='')
			}
		} else {
			if (! is.null(opt$sampleorder)) {
				for (ss in colorder) {
					if (! (ss %in% colnames(expr_data))) stop("\nNo matching column for sample '",ss,"' (from --sampleorder) was found in --expr_data",sep='')
				}
			}
			if (! (opt$samplename %in% colnames(expr_data))) stop("\nNo matching column for sample '",opt$samplename,"' in --expr_data",sep='')
		}
	}
	if (opt$GEOmean == TRUE) {
		expr_data = log2(expr_data + 1)		# add pseudocount of 1
	}
}

cat("DONE\n")

# get GO data from provided database using biomart
if (opt$ensembl == FALSE) {
	cat("Downloading GO dataset",opt$GOdataset,"(",opt$GOmart,") from",opt$GOhost,"using biomart...\n")
	m = useMart(opt$GOmart, dataset=opt$GOdataset, host=opt$GOhost)
} else {
	cat("Downloading GO dataset",opt$GOdataset," from ensembl using biomart...\n")
	m = useMart("ensembl",dataset=opt$GOdataset)
}	

godata = getBM(attributes = c("ensembl_gene_id","go_id"), mart = m)

# filter out genes with no GO terms
godata <- godata[godata$go_id != '',]

# reformat
golist <- by(godata$go_id, godata$ensembl_gene_id, function(x) as.character(x))

# get the list of background genes
if (is.null(opt$background)) {
	background = sort(unique(as.character(godata$ensembl_gene_id)))
} else {
	backgrdf = read.table(opt$background, header=FALSE, sep="\t", stringsAsFactors = FALSE)
	background = backgrdf[,1]
}

# for each group of genes (column 2 in genelist), run GO analysis
genegroups = unique(genelist[,2])
for (gg in genegroups) {
	gsubset = genelist[genelist$V2 == gg,]
	genell = gsubset[,1]
	ngenes = nrow(gsubset)
	cat("Running topGO analysis on the",ngenes,gg,"genes...")
	genelbl = factor(as.integer(background %in% genell))
	names(genelbl) = background
	
	# run GO analysis
	tryCatch({
		res = runtopGO(genelbl, golist, algorithm = opt$algorithm, statistic = opt$statistic)
		GOdata = res[[1]]; results = res[[2]]
	}, error = function(err) {
		cat("Error occurred while running topGO, see log file",paste(opt$outprefix,"_ext_logfile.txt",sep=''),"\n")
	})
	cat("DONE\n")
	
	# save p-values to big matrix
	pvals_list = data.frame(row.names = names(attributes(results)$score), pvals = attributes(results)$score, stringsAsFactors = FALSE)
	colnames(pvals_list)[1] = gg
	if (exists("pvals_all")) {
		pvals_all = merge(pvals_all,pvals_list,by=0, all = TRUE)
		rownames(pvals_all) = pvals_all$Row.names
		pvals_all$Row.names = NULL
	} else {
		pvals_all = pvals_list
	}

	# get out significant GO-terms, passing thresholds
	passfail = attributes(results)$score < opt$pval
	numpass = length(passfail[passfail == TRUE])
	if (numpass > 0) {
		results_table = GenTable(GOdata, elimFisher = results, topNodes = numpass, numChar = 100)
		colnames(results_table)[6] = "pval"
		results_table = cbind(samp = gg, results_table)
				
		if (exists("all_sig_res")) {
			all_sig_res = rbind(all_sig_res, results_table)
		} else {
			all_sig_res = results_table
		}
		
		# make flowchart-type plot of sig. GO terms
		cat("\t- Outputting summary graph of results for",gg,"to: ")
		printGraph(GOdata, results, firstSigNodes = numpass, fn.prefix = paste(opt$outprefix,'_graph_',gg,sep=''), useInfo = "all", pdfSW = TRUE)
	}	
	cat("\t- ",gg," summary: ",numpass," significant GO-terms with p-value < ",opt$pval," detected\n",sep='')
}

all_sig_res$pval[all_sig_res$pval == "< 1e-30"] = "0"
all_sig_res$pval = as.numeric(all_sig_res$pval)
all_sig_res$logpval = -1*log10(all_sig_res$pval)
annot = unique(all_sig_res[,c(2,3)])

# do I also make plot of top N significant terms?
maketopN = FALSE
if (! is.null(opt$topN)) {
	maketopN = TRUE
	tmp = all_sig_res[order(all_sig_res$samp, -all_sig_res$logpval),]
	topNlist = c()
	for (samp in unique(all_sig_res$samp)) {
		ss = tmp[tmp$samp == samp,]
		topNlist = c(topNlist, ss[1:opt$topN,]$GO.ID)
	}
	topNlist = unique(topNlist)

	# merge back to full dataset
	all_sig_res = as.data.frame(cbind(all_sig_res,ifelse(all_sig_res$GO.ID %in% topNlist, 1, 0)))
	colnames(all_sig_res)[length(colnames(all_sig_res))] = paste('top',opt$topN,sep='')
}

# make p-value plots
cat("Making plots of GO-terms and p-values...")
if (samplesprovided == FALSE) {
	# make barchart
	mydata = all_sig_res
	mydata$logpval[mydata$logpval > opt$pvalmax] = opt$pvalmax
	mydata = mydata[order(mydata$logpval),]
	mydata$xlbls = paste(mydata$GO.ID,mydata$Term,sep=': ')
	pdf(paste(opt$outprefix,"_pval_barchart.pdf",sep=''), width = 8, height = max(4,nrow(mydata)/50 + 2))
	a = ggplot(data=mydata, aes(x=xlbls, y=logpval))
	a = a + geom_bar(stat="identity", color="black", position = position_dodge()) + coord_flip() 
	a = a + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("p-value") + xlab("GO term") + theme(text = element_text(size=6))
	a = a + scale_x_discrete(limits=mydata$xlbls)
	print(a)
	graphics.off()	
	pvalorder = mydata$GO.ID
	
	if (maketopN == TRUE) {
		mydata = mydata[mydata$GO.ID %in% topNlist,]
		pdf(paste(opt$outprefix,"_pval_barchart_top",opt$topN,".pdf",sep=''), width = 8, height = max(4,nrow(mydata)/50 + 2))
		a = ggplot(data=mydata, aes(x=xlbls, y=logpval))
		a = a + geom_bar(stat="identity", color="black", position = position_dodge()) + coord_flip() 
		a = a + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("p-value") + xlab("GO term") + theme(text = element_text(size=6))
		a = a + scale_x_discrete(limits=mydata$xlbls)
		print(a)
		graphics.off()	
		pvalorder_topN = mydata$GO.ID
	}
} else {
	# make heatmap of the significant terms
	sig_terms = unique(all_sig_res$GO.ID)
	pvals_sig = pvals_all[rownames(pvals_all) %in% sig_terms,]
	pvals_sig_log10 = -1 * log10(pvals_sig)
	to_output = pvals_sig_log10
	to_output = merge(annot, to_output, by.x=1, by.y = 0)
	write.table(to_output,file=paste(opt$outprefix,"_log10_pval_matrix.txt",sep=''),col.names = TRUE,row.names = FALSE,quote = F,sep='\t')
	
	if (is.null(opt$pvalmax)) opt$pvalmax = max(pvals_sig_log10)
	breaksList = seq(0, min(opt$pvalmax,max(pvals_sig_log10)), by = min(opt$pvalmax,max(pvals_sig_log10)) / 50)
	paletteLength = length(breaksList)
	hmcolorscale = colorRampPalette(c("white","#FFFEC6","#F9E19B","#F2B26E","#E3754F","#D1382C"))(paletteLength)
	
	row_annot = annot[match(rownames(pvals_sig_log10),annot$GO.ID),]
	cluster_cols = TRUE
	if (! is.null(opt$sampleorder)) {
		pvals_sig_log10 = pvals_sig_log10[,colorder]
		cluster_cols = FALSE
		pval_colorder = colorder
		pval_colorder_topN = colorder
		pval_colorder_bool = colorder
		pval_colorder_bool_topN = colorder
	}
	
	fontsize_col = min(14,3 + 50/ncol(pvals_sig_log10))
	fontsize_row = min(3 + 50/nrow(pvals_sig_log10))
	pdf(paste(opt$outprefix,"_log10_pval_heatmap.pdf",sep=''), width = 8, height = 8)
	p1 = pheatmap(pvals_sig_log10, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", cluster_cols=cluster_cols, fontsize_row = fontsize_row, fontsize_col = fontsize_col, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
	show(p1)
	graphics.off()
		
	# write p-values matrix to file
	write.table(pvals_sig_log10,file=paste(opt$outprefix,"_log10_pval_matrix.txt",sep=''),col.names = TRUE,row.names = TRUE,quote = F,sep='\t')
	
	# also, get order of rows and columns here so we can replicate with the expression data
	pval_roworder = rownames(pvals_sig_log10)[p1$tree_row$order]
	if (is.null(opt$sampleorder)) pval_colorder = colnames(pvals_sig_log10)[p1$tree_col$order]
		
	# make subset plot if requested
	if (maketopN == TRUE) {
		pvals_sig_log10_ss = pvals_sig_log10[rownames(pvals_sig_log10) %in% topNlist,]
		row_annot_topN = row_annot[row_annot$GO.ID %in% topNlist,]
		row_annot_topN = row_annot_topN[match(rownames(pvals_sig_log10_ss),row_annot_topN$GO.ID),]
		pdf(paste(opt$outprefix,"_log10_pval_heatmap_top",opt$topN,".pdf",sep=''), width = 8, height = 8)
		p1 = pheatmap(pvals_sig_log10_ss, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", cluster_cols = cluster_cols, fontsize_row = fontsize_row, fontsize_col = fontsize_col, labels_row = paste(row_annot_topN$GO.ID,row_annot_topN$Term,sep=': '))
		show(p1)
		graphics.off()
		pval_roworder_topN = rownames(pvals_sig_log10_ss)[p1$tree_row$order]
		if (is.null(opt$sampleorder)) pval_colorder_topN = colnames(pvals_sig_log10_ss)[p1$tree_col$order]
	}

	# make boolean version
	cutoff = -1 * log10(opt$pval)
	breaksList = seq(0, 1, by = 1 / 50)
	paletteLength = length(breaksList)
	hmcolorscale = colorRampPalette(c("white","#FFFEC6","#F9E19B","#F2B26E","#E3754F","#D1382C"))(paletteLength)

	pvals_sig_bool = ifelse(pvals_sig_log10 > cutoff,1,0)
	pdf(paste(opt$outprefix,"_log10_pval_heatmap_bool.pdf",sep=''), width = 8, height = 8)
	p1 = pheatmap(pvals_sig_bool, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", cluster_cols = cluster_cols, fontsize_row = fontsize_row, fontsize_col = fontsize_col, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
	show(p1)
	graphics.off()
		
	# also, get order of rows and columns here so we can replicate with the expression data
	pval_roworder_bool = rownames(pvals_sig_bool)[p1$tree_row$order]
	if (is.null(opt$sampleorder)) pval_colorder_bool = colnames(pvals_sig_bool)[p1$tree_col$order]
	
	# also make a version where rows are just binary sorted, instead of clustered
	pvals_sig_bool_tmp = pvals_sig_bool
	pvals_sig_bool_tmp = data.frame(pvals_sig_bool_tmp[,pval_colorder_bool])
	pvals_sig_bool_tmp = data.frame(pvals_sig_bool_tmp[do.call(order,c(as.data.frame(pvals_sig_bool_tmp), list(decreasing=TRUE))),])
	row_annot_tmp = row_annot[match(rownames(pvals_sig_bool_tmp),row_annot$GO.ID),]
	
	pdf(paste(opt$outprefix,"_log10_pval_heatmap_bool_sortrows.pdf",sep=''), width = 8, height = 8)
	p1 = pheatmap(pvals_sig_bool_tmp, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", cluster_cols = FALSE, cluster_rows = FALSE, fontsize_row = fontsize_row, fontsize_col = fontsize_col, labels_row = paste(row_annot_tmp$GO.ID,row_annot_tmp$Term,sep=': '))
	show(p1)
	graphics.off()
	
	# make subset plot if requested
	if (maketopN == TRUE) {
		pvals_sig_bool = pvals_sig_bool[rownames(pvals_sig_bool) %in% topNlist,]
		row_annot_topN = row_annot_topN[match(rownames(pvals_sig_bool),row_annot_topN$GO.ID),]
		pdf(paste(opt$outprefix,"_log10_pval_heatmap_bool_top",opt$topN,".pdf",sep=''), width = 8, height = 8)
		p1 = pheatmap(pvals_sig_bool, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", cluster_cols = cluster_cols, fontsize_row = fontsize_row, fontsize_col = fontsize_col, labels_row = paste(row_annot_topN$GO.ID,row_annot_topN$Term,sep=': '))
		show(p1)
		graphics.off()
		pval_roworder_bool_topN = rownames(pvals_sig_bool)[p1$tree_row$order]
		if (is.null(opt$sampleorder)) pval_colorder_bool_topN = colnames(pvals_sig_bool)[p1$tree_col$order]
		
		# also make a version where the GO terms are sorted based on left-to-right order of 
		# samples, rather than clustering them
		pvals_sig_bool_tmp = pvals_sig_bool
		pvals_sig_bool_tmp = data.frame(pvals_sig_bool_tmp[,pval_colorder_bool_topN])
		pvals_sig_bool_tmp = data.frame(pvals_sig_bool_tmp[do.call(order,c(as.data.frame(pvals_sig_bool_tmp), list(decreasing=TRUE))),])
		row_annot_topN_tmp = row_annot_topN[match(rownames(pvals_sig_bool_tmp),row_annot_topN$GO.ID),]
		
		pdf(paste(opt$outprefix,"_log10_pval_heatmap_bool_top",opt$topN,"_sortrows.pdf",sep=''), width = 8, height = 8)
		p1 = pheatmap(pvals_sig_bool_tmp, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", cluster_cols = FALSE, cluster_rows = FALSE, fontsize_row = fontsize_row, fontsize_col = fontsize_col, labels_row = paste(row_annot_topN_tmp$GO.ID,row_annot_topN_tmp$Term,sep=': '))
		show(p1)
		graphics.off()
	}	
}
cat("DONE\n")

# for each significant GO term in the list, output list of all genes in each input list that were
# in that GO term
all_sig_res$genes_in_term = rep("",length(all_sig_res$samp))
for (i in 1:nrow(all_sig_res)) {
	gogenes = genesInTerm(GOdata, all_sig_res[i,]$GO.ID)[[1]]
	ggenes = genelist[genelist$V2 == all_sig_res[i,]$samp,]$V1
	gogenes = paste(gogenes[gogenes %in% ggenes], collapse = ',')
	all_sig_res$genes_in_term[i] = gogenes
}

# if expression data was provided, make those plots too
if (! is.null(opt$expr_data)) {
	cat("Since expression data also provided, making plots of expression over significant GO-terms...")

	if (is.null(opt$samp_map)) {
		mydata_mean = expr_data[,colnames(expr_data) %in% genegroups,drop=FALSE]
		if (opt$single_cell == TRUE) mydata_frac = ifelse(mydata_mean > 0,1,0)	
	} else {	
		mydata = merge(samp_map,t(expr_data),by=0)
		rownames(mydata) = mydata$Row.names
		mydata$Row.names = NULL
		mydata_mean = aggregate(mydata[, 2:ncol(mydata)], list(mydata[,1]), FUN=mean)
		rownames(mydata_mean) = mydata_mean$Group.1
		mydata_mean$Group.1 = NULL
		mydata_mean = t(mydata_mean)
		
		if (opt$single_cell == TRUE) {
			if (opt$includezeros == FALSE) {
				mydata_mean_nozero = aggregate(mydata[, 2:ncol(mydata)], list(mydata[,1]), FUN=(function(x){ifelse(sum(x==0)>0 & sum(x !=0) >0, mean(x[x>0]), mean(x))}))
				rownames(mydata_mean_nozero) = mydata_mean_nozero$Group.1
				mydata_mean_nozero$Group.1 = NULL
				mydata_mean_nozero = t(mydata_mean_nozero)		
			}
			mydata_frac = aggregate(mydata[, 2:ncol(mydata)] > 0, list(mydata[,1]), mean)
			rownames(mydata_frac) = mydata_frac$Group.1
			mydata_frac$Group.1 = NULL
			mydata_frac = t(mydata_frac)
		}
	}
		
	# now get averages across each GO-term; since multiple go terms can have the same genes,
	# will do it one-by-one in loop
	sig_terms = unique(all_sig_res$GO.ID)
	allgotermgenes = genesInTerm(GOdata, sig_terms)
	for (goterm in sig_terms) {
		genes = allgotermgenes[[goterm]]
		# keep only the significant ones
		if (opt$useallGOtermgenes == FALSE) genes = genes[genes %in% genelist$V1]
				
		res = colMeans(mydata_mean[rownames(mydata_mean) %in% genes,,drop=FALSE],na.rm=TRUE)	
		res = as.data.frame(t(res))
		rownames(res) = goterm
		if (exists("mydata_mean_final")) {
			mydata_mean_final = rbind(mydata_mean_final,res)
		} else {
			mydata_mean_final = res
		}
	
		if (opt$single_cell == TRUE) {
			res = colMeans(mydata_frac[rownames(mydata_frac) %in% genes,,drop=FALSE],na.rm=TRUE)
			res = as.data.frame(t(res))
			rownames(res) = goterm
			if (exists("mydata_frac_final")) {
				mydata_frac_final = rbind(mydata_frac_final,res)
			} else {
				mydata_frac_final = res
			}
			if (opt$includezeros == FALSE) {
				res = colMeans(mydata_mean_nozero[rownames(mydata_mean_nozero) %in% genes,,drop=FALSE],na.rm=TRUE)
				res = as.data.frame(t(res))
				rownames(res) = goterm
				if (exists("mydata_mean_nozero_final")) {
					mydata_mean_nozero_final = rbind(mydata_mean_nozero_final,res)
				} else {
					mydata_mean_nozero_final = res
				}
			}
		}	
	}	
	
	# prepare final datasets, with expression data, to output as table	
	res = wide_to_long(mydata_mean_final)
	if (opt$GEOmean == TRUE) { colnames(res) = c('samp','avg. log2(expr)','GO.ID') } else { colnames(res) = c('samp','avg. expr','GO.ID') } 
	all_sig_res = merge(all_sig_res, res, by=c("GO.ID","samp"))

	# write mean expression matrix to file
	write.table(mydata_mean_final,file=paste(opt$outprefix,"_avg_expr_matrix.txt",sep=''),col.names = TRUE,row.names = TRUE,quote = F,sep='\t')

	if (opt$single_cell == TRUE) {
		res = wide_to_long(mydata_frac_final)
		colnames(res) = c('samp','avg. frac. nonzero','GO.ID') 
		all_sig_res = merge(all_sig_res, res, by=c("GO.ID","samp"))
		
		# write mean expression matrix to file
		write.table(mydata_frac_final,file=paste(opt$outprefix,"_avg_frac_nonzero_matrix.txt",sep=''),col.names = TRUE,row.names = TRUE,quote = F,sep='\t')

		if (opt$includezeros == FALSE) {
			res = wide_to_long(mydata_mean_nozero_final)
			if (opt$GEOmean == TRUE) { colnames(res) = c('samp','avg. log2(expr) - excl. zeros','GO.ID') } else { colnames(res) = c('samp','avg. expr - excl. zeros','GO.ID') } 
			all_sig_res = merge(all_sig_res, res, by=c("GO.ID","samp"))
			
			# write mean expression matrix to file
			write.table(mydata_mean_nozero_final,file=paste(opt$outprefix,"_avg_expr_exclzero_matrix.txt",sep=''),col.names = TRUE,row.names = TRUE,quote = F,sep='\t')
		}
	}
	
	row_annot = annot[match(rownames(mydata_mean_final),annot$GO.ID),]
	if (opt$single_cell == FALSE && samplesprovided == TRUE) {
		if (is.null(opt$fillupper)) opt$fillupper = max(mydata_mean_final,na.rm=TRUE)
		if (is.null(opt$filllower)) opt$filllower = min(mydata_mean_final,na.rm=TRUE)
		chunksize = (opt$fillupper-opt$filllower) / 50
		theobreaks = seq(-1*max(opt$fillupper, abs(opt$filllower)),max(opt$fillupper, abs(opt$filllower)),chunksize)
		hmcolorscale = colorRampPalette(c("#3D42B6","#4F74B0","#80ABCD","#B4D8E7","#E4F2F7","white","#FFFEC6","#F9E19B","#F2B26E","#E3754F","#D1382C"))(length(theobreaks))		
		
		if (max(opt$fillupper, abs(opt$filllower)) == opt$fillupper) {
			breaksList = tail(theobreaks,51)
			hmcolorscale = tail(hmcolorscale,51)
		} else {
			breaksList = head(theobreaks,51)
			hmcolorscale = head(hmcolorscale,51)
		}
				
		pdf(paste(opt$outprefix,"_expr_heatmap.pdf",sep=''), width = 8, height = 8)
		p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
		show(p1)
		graphics.off()
		
		# make version sorted according to p-values heatmap
		mydata_mean_final = mydata_mean_final[pval_roworder, pval_colorder]
		row_annot = annot[match(rownames(mydata_mean_final),annot$GO.ID),]
		pdf(paste(opt$outprefix,"_expr_heatmap_pvalsort.pdf",sep=''), width = 8, height = 8)
		p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, cluster_cols=FALSE, cluster_rows=FALSE, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
		show(p1)
		graphics.off()
		
		# make version sorted according to boolean p-values heatmap
		mydata_mean_final = mydata_mean_final[pval_roworder_bool, pval_colorder_bool]
		row_annot = annot[match(rownames(mydata_mean_final),annot$GO.ID),]
		pdf(paste(opt$outprefix,"_expr_heatmap_pvalsortbool.pdf",sep=''), width = 8, height = 8)
		p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, cluster_cols=FALSE, cluster_rows=FALSE, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
		show(p1)
		graphics.off()	
		
		# make subset plots if requested
		if (maketopN == TRUE) {
			mydata_mean_final = mydata_mean_final[pval_roworder_topN,pval_colorder_topN]
			row_annot = annot[match(rownames(mydata_mean_final),annot$GO.ID),]
			pdf(paste(opt$outprefix,"_expr_heatmap_top",opt$topN,".pdf",sep=''), width = 8, height = 8)
			p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
			show(p1)
			graphics.off()

			pdf(paste(opt$outprefix,"_expr_heatmap_pvalsort_top",opt$topN,".pdf",sep=''), width = 8, height = 8)
			p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, cluster_cols=FALSE, cluster_rows=FALSE, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
			show(p1)
			graphics.off()

			mydata_mean_final = mydata_mean_final[pval_roworder_bool_topN, pval_colorder_bool_topN]
			row_annot = annot[match(rownames(mydata_mean_final),annot$GO.ID),]
			pdf(paste(opt$outprefix,"_expr_heatmap_pvalsortbool_top",opt$topN,".pdf",sep=''), width = 8, height = 8)
			p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, cluster_cols=FALSE, cluster_rows=FALSE, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
			show(p1)
			graphics.off()	
		}				
	} else if (opt$single_cell == TRUE) {
		
		# perform clustering to get optimal order of x and y values for plot (comparable to the heatmap)
		if (length(sig_terms) > 40) {
			cat("\n--single_cell mode requested, but over 100 significant terms (",length(sig_terms),") were identified - dot plot will be replaced by heatmap of avg. values (zeros included)\n")
			if (is.null(opt$fillupper)) opt$fillupper = max(mydata_mean_final,na.rm=TRUE)
			if (is.null(opt$filllower)) opt$filllower = min(mydata_mean_final,na.rm=TRUE)
			chunksize = (opt$fillupper-opt$filllower) / 50
			theobreaks = seq(-1*max(opt$fillupper, abs(opt$filllower)),max(opt$fillupper, abs(opt$filllower)),chunksize)
			hmcolorscale = colorRampPalette(c("#3D42B6","#4F74B0","#80ABCD","#B4D8E7","#E4F2F7","white","#FFFEC6","#F9E19B","#F2B26E","#E3754F","#D1382C"))(length(theobreaks))		
		
			if (max(opt$fillupper, abs(opt$filllower)) == opt$fillupper) {
				breaksList = tail(theobreaks,51)
				hmcolorscale = tail(hmcolorscale,51)
			} else {
				breaksList = head(theobreaks,51)
				hmcolorscale = head(hmcolorscale,51)
			}
		
			pdf(paste(opt$outprefix,"_expr_heatmap.pdf",sep=''), width = 8, height = 8)
			p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
			show(p1)
			graphics.off()
		
			# make version sorted according to p-values heatmap
			mydata_mean_final = mydata_mean_final[pval_roworder, pval_colorder]
			row_annot = annot[match(rownames(mydata_mean_final),annot$GO.ID),]
			pdf(paste(opt$outprefix,"_expr_heatmap_pvalsort.pdf",sep=''), width = 8, height = 8)
			p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, cluster_cols=FALSE, cluster_rows=FALSE, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
			show(p1)
			graphics.off()
		
			# make version sorted according to boolean p-values heatmap
			mydata_mean_final = mydata_mean_final[pval_roworder_bool, pval_colorder_bool]
			row_annot = annot[match(rownames(mydata_mean_final),annot$GO.ID),]
			pdf(paste(opt$outprefix,"_expr_heatmap_pvalsortbool.pdf",sep=''), width = 8, height = 8)
			p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, cluster_cols=FALSE, cluster_rows=FALSE, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
			show(p1)
			graphics.off()	
		} else {
			# make dotplot - reshape to long df
			if (opt$includezeros == FALSE) { expr_touse = mydata_mean_nozero_final } else { expr_touse = mydata_mean_final }
			mydata_final = data.frame(samp = c(), goterm = c(), mmean = c(), frac = c())
			for (samp in colnames(mydata_mean_final)) {
				mydata_tmp = data.frame(samp = rep(samp,nrow(expr_touse)), goterm = rownames(expr_touse), mmean = unname(expr_touse[,samp]), frac = unname(mydata_frac_final[,samp]))
				mydata_final = rbind(mydata_final,mydata_tmp)
			}

			# cluster samples by overall expression; we'll use this order to make dotplot
			mydata_mean_forclus = mydata_mean_final
			dd = dist(mydata_mean_forclus, method = 'euclidean')
			myhclust = hclust(dd, method = 'complete')
			yorder = myhclust$labels[myhclust$order]
		
			dd = dist(t(mydata_mean_forclus), method = 'euclidean')
			myhclust = hclust(dd, method = 'complete')
			xorder = myhclust$labels[myhclust$order]
			
			xorder = data.frame(samp = xorder, x = 1:length(xorder))
			yorder = data.frame(goterm = yorder, y = 1:length(yorder))
		
			mydata_final = merge(mydata_final,xorder,by="samp")
			mydata_final = merge(mydata_final,yorder,by="goterm")
		
			# set upper limits on mmean and frac, if requested
			if (is.null(opt$fillupper)) {
				fmax = max(mydata_final$mmean,na.rm=TRUE)
			} else {
				fmax = opt$fillupper
				mydata_final$mmean[mydata_final$mmean > fmax] = fmax
			}
			if (is.null(opt$filllower)) {
				fmin = min(mydata_final$mmean,na.rm=TRUE)
			} else {
				fmin = opt$filllower
				mydata_final$mmean[mydata_final$mmean < fmin] = fmin
			}
			if (is.null(opt$sizeupper)) {
				smax = max(mydata_final$frac,na.rm=TRUE)
			} else {
				smax = opt$sizeupper
				mydata_final$frac[mydata_final$frac > smax] = smax
			}
		
			row_annot = annot[match(unique(mydata_final$goterm),annot$GO.ID),]
			row_annot$order = unique(mydata_final$y)
			row_annot = row_annot[order(row_annot$order),]

			pdf(paste(opt$outprefix,'_expr_dotplot.pdf',sep=''), width = 20, height = 15)
			p1 = ggplot(mydata_final, aes(x=x, y=y)) + geom_point(aes(size=frac, fill=mmean),pch=21) + theme_bw()
			p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax)) + scale_fill_gradientn(colours=c("royalblue", "khaki", "firebrick"), limits=c(fmin,fmax))
			p1 = p1 + scale_y_continuous(breaks=1:nrow(yorder), labels=paste(row_annot$GO.ID,row_annot$Term,sep=': '), limits = c(0.5,nrow(yorder)+0.5)) + theme(panel.grid.minor.y = element_blank())
			p1 = p1 + scale_x_continuous(breaks=1:nrow(xorder), labels=xorder$samp, limits = c(0.5,nrow(xorder)+0.5)) + theme(panel.grid.minor.x = element_blank()) + xlab("Samples") + ylab("GO-terms")
			p1 = p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
			if (opt$GEOmean == FALSE) p1 = p1 + labs(fill = "avg. expr", size = "frac. non-zero")
			if (opt$GEOmean == TRUE) p1 = p1 + labs(fill = "avg. log2(expr)", size = "frac. non-zero")
			show(p1)
			graphics.off()
				
			# make versions in same order as p-value heatmaps
			newyorder = data.frame(term = rev(pval_roworder), py = 1:length(pval_roworder))
			newxorder = data.frame(samp = pval_colorder, px = 1:length(pval_colorder))
			mydata_final = merge(mydata_final, newyorder, by = 1)
			mydata_final = merge(mydata_final, newxorder, by.x = 2, by.y = 1)
		
			row_annot = row_annot[match(newyorder$term,row_annot$GO.ID),]

			pdf(paste(opt$outprefix,'_expr_dotplot_pvalsort.pdf',sep=''), width = 20, height = 15)
			p1 = ggplot(mydata_final, aes(x=px, y=py)) + geom_point(aes(size=frac, fill=mmean),pch=21) + theme_bw()
			p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax)) + scale_fill_gradientn(colours=c("royalblue", "khaki", "firebrick"), limits=c(fmin,fmax))
			p1 = p1 + scale_y_continuous(breaks=1:nrow(newyorder), labels=paste(row_annot$GO.ID,row_annot$Term,sep=': '), limits = c(0.5,nrow(yorder)+0.5)) + theme(panel.grid.minor.y = element_blank())
			p1 = p1 + scale_x_continuous(breaks=1:nrow(newxorder), labels=newxorder$samp, limits = c(0.5,nrow(xorder)+0.5)) + theme(panel.grid.minor.x = element_blank()) + xlab("Samples") + ylab("GO-terms")
			p1 = p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
			if (opt$GEOmean == FALSE) p1 = p1 + labs(fill = "avg. expr", size = "frac. non-zero")
			if (opt$GEOmean == TRUE) p1 = p1 + labs(fill = "avg. log2(expr)", size = "frac. non-zero")
			show(p1)
			graphics.off()

			newyorder = data.frame(term = rev(pval_roworder_bool), bpy = 1:length(pval_roworder_bool))
			newxorder = data.frame(samp = pval_colorder_bool, bpx = 1:length(pval_colorder_bool))
			mydata_final = merge(mydata_final, newyorder, by.x = 2, by.y = 1)
			mydata_final = merge(mydata_final, newxorder, by.x = 2, by.y = 1)
		
			row_annot = row_annot[match(newyorder$term,row_annot$GO.ID),]

			pdf(paste(opt$outprefix,'_expr_dotplot_pvalsortbool.pdf',sep=''), width = 20, height = 15)
			p1 = ggplot(mydata_final, aes(x=x, y=y)) + geom_point(aes(size=frac, fill=mmean),pch=21) + theme_bw()
			p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax)) + scale_fill_gradientn(colours=c("royalblue", "khaki", "firebrick"), limits=c(fmin,fmax))
			p1 = p1 + scale_y_continuous(breaks=1:nrow(newyorder), labels=paste(row_annot$GO.ID,row_annot$Term,sep=': '), limits = c(0.5,nrow(yorder)+0.5)) + theme(panel.grid.minor.y = element_blank())
			p1 = p1 + scale_x_continuous(breaks=1:nrow(newxorder), labels=newxorder$samp, limits = c(0.5,nrow(xorder)+0.5)) + theme(panel.grid.minor.x = element_blank()) + xlab("Samples") + ylab("GO-terms")
			p1 = p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
			if (opt$GEOmean == FALSE) p1 = p1 + labs(fill = "avg. expr", size = "frac. non-zero")
			if (opt$GEOmean == TRUE) p1 = p1 + labs(fill = "avg. log2(expr)", size = "frac. non-zero")
			show(p1)
			graphics.off()		
		}
		
		if (maketopN == TRUE) {
			mydata_mean_final = mydata_mean_final[pval_roworder_topN,pval_colorder_topN]
			row_annot = annot[match(rownames(mydata_mean_final),annot$GO.ID),]

			if (length(pval_roworder_topN) > 40) {
				pdf(paste(opt$outprefix,"_expr_heatmap_top",opt$topN,".pdf",sep=''), width = 8, height = 8)
				p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
				show(p1)
				graphics.off()

				pdf(paste(opt$outprefix,"_expr_heatmap_pvalsort_top",opt$topN,".pdf",sep=''), width = 8, height = 8)
				p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, cluster_cols=FALSE, cluster_rows=FALSE, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
				show(p1)
				graphics.off()

				mydata_mean_final = mydata_mean_final[pval_roworder_bool_topN, pval_colorder_bool_topN]
				row_annot = annot[match(rownames(mydata_mean_final),annot$GO.ID),]
				pdf(paste(opt$outprefix,"_expr_heatmap_pvalsortbool_top",opt$topN,".pdf",sep=''), width = 8, height = 8)
				p1 = pheatmap(mydata_mean_final, color = hmcolorscale, breaks = breaksList, show_rownames=TRUE, na_col="grey50", fontsize_row = 3, cluster_cols=FALSE, cluster_rows=FALSE, labels_row = paste(row_annot$GO.ID,row_annot$Term,sep=': '))
				show(p1)
				graphics.off()	
			} else {
				mydata_frac_final = mydata_frac_final[pval_roworder_topN,pval_colorder_topN]
				if (opt$includezeros == FALSE) { 
					mydata_mean_nozero_final = mydata_mean_nozero_final[pval_roworder_topN,pval_colorder_topN]
					expr_touse = mydata_mean_nozero_final 
				} else { 
					expr_touse = mydata_mean_final 
				}
				mydata_final = data.frame(samp = c(), goterm = c(), mmean = c(), frac = c())
				for (samp in colnames(mydata_mean_final)) {
					mydata_tmp = data.frame(samp = rep(samp,nrow(expr_touse)), goterm = rownames(expr_touse), mmean = unname(expr_touse[,samp]), frac = unname(mydata_frac_final[,samp]))
					mydata_final = rbind(mydata_final,mydata_tmp)
				}

				mydata_mean_forclus = mydata_mean_final
				dd = dist(mydata_mean_forclus, method = 'euclidean')
				myhclust = hclust(dd, method = 'complete')
				yorder = myhclust$labels[myhclust$order]
		
				dd = dist(t(mydata_mean_forclus), method = 'euclidean')
				myhclust = hclust(dd, method = 'complete')
				xorder = myhclust$labels[myhclust$order]
			
				xorder = data.frame(samp = xorder, x = 1:length(xorder))
				yorder = data.frame(goterm = yorder, y = 1:length(yorder))
		
				mydata_final = merge(mydata_final,xorder,by="samp")
				mydata_final = merge(mydata_final,yorder,by="goterm")
			
				# set upper limits on mmean and frac, if requested
				if (is.null(opt$fillupper)) {
					fmax = max(mydata_final$mmean,na.rm=TRUE)
				} else {
					fmax = opt$fillupper
					mydata_final$mmean[mydata_final$mmean > fmax] = fmax
				}
				if (is.null(opt$filllower)) {
					fmin = min(mydata_final$mmean,na.rm=TRUE)
				} else {
					fmin = opt$filllower
					mydata_final$mmean[mydata_final$mmean < fmin] = fmin
				}
				if (is.null(opt$sizeupper)) {
					smax = max(mydata_final$frac,na.rm=TRUE)
				} else {
					smax = opt$sizeupper
					mydata_final$frac[mydata_final$frac > smax] = smax
				}

				row_annot = annot[match(unique(mydata_final$goterm),annot$GO.ID),]
				row_annot$order = unique(mydata_final$y)
				row_annot = row_annot[order(row_annot$order),]
			
				pdf(paste(opt$outprefix,'_expr_dotplot_top',opt$topN,'.pdf',sep=''), width = 20, height = 15)
				p1 = ggplot(mydata_final, aes(x=x, y=y)) + geom_point(aes(size=frac, fill=mmean),pch=21) + theme_bw()
				p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax)) + scale_fill_gradientn(colours=c("royalblue", "khaki", "firebrick"), limits=c(fmin,fmax))
				p1 = p1 + scale_y_continuous(breaks=1:nrow(yorder), labels=paste(row_annot$GO.ID,row_annot$Term,sep=': '), limits = c(0.5,nrow(yorder)+0.5)) + theme(panel.grid.minor.y = element_blank())
				p1 = p1 + scale_x_continuous(breaks=1:nrow(xorder), labels=xorder$samp, limits = c(0.5,nrow(xorder)+0.5)) + theme(panel.grid.minor.x = element_blank()) + xlab("Samples") + ylab("GO-terms")
				p1 = p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
				if (opt$GEOmean == FALSE) p1 = p1 + labs(fill = "avg. expr", size = "frac. non-zero")
				if (opt$GEOmean == TRUE) p1 = p1 + labs(fill = "avg. log2(expr)", size = "frac. non-zero")
				show(p1)
				graphics.off()
			
				# make versions in same order as p-value heatmaps
				newyorder = data.frame(term = rev(pval_roworder_topN), py = 1:length(pval_roworder_topN))
				newxorder = data.frame(samp = pval_colorder_topN, px = 1:length(pval_colorder_topN))
				mydata_final = merge(mydata_final, newyorder, by = 1)
				mydata_final = merge(mydata_final, newxorder, by.x = 2, by.y = 1)
		
				row_annot = row_annot[match(newyorder$term,row_annot$GO.ID),]

				pdf(paste(opt$outprefix,'_expr_dotplot_pvalsort_top',opt$topN,'.pdf',sep=''), width = 20, height = 15)
				p1 = ggplot(mydata_final, aes(x=px, y=py)) + geom_point(aes(size=frac, fill=mmean),pch=21) + theme_bw()
				p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax)) + scale_fill_gradientn(colours=c("royalblue", "khaki", "firebrick"), limits=c(fmin,fmax))
				p1 = p1 + scale_y_continuous(breaks=1:nrow(newyorder), labels=paste(row_annot$GO.ID,row_annot$Term,sep=': '), limits = c(0.5,nrow(yorder)+0.5)) + theme(panel.grid.minor.y = element_blank())
				p1 = p1 + scale_x_continuous(breaks=1:nrow(newxorder), labels=newxorder$samp, limits = c(0.5,nrow(xorder)+0.5)) + theme(panel.grid.minor.x = element_blank()) + xlab("Samples") + ylab("GO-terms")
				p1 = p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
				if (opt$GEOmean == FALSE) p1 = p1 + labs(fill = "avg. expr", size = "frac. non-zero")
				if (opt$GEOmean == TRUE) p1 = p1 + labs(fill = "avg. log2(expr)", size = "frac. non-zero")
				show(p1)
				graphics.off()

				newyorder = data.frame(term = rev(pval_roworder_bool_topN), bpy = 1:length(pval_roworder_bool_topN))
				newxorder = data.frame(samp = pval_colorder_bool_topN, bpx = 1:length(pval_colorder_bool_topN))
				mydata_final = merge(mydata_final, newyorder, by.x = 2, by.y = 1)
				mydata_final = merge(mydata_final, newxorder, by.x = 2, by.y = 1)
		
				row_annot = row_annot[match(newyorder$term,row_annot$GO.ID),]

				pdf(paste(opt$outprefix,'_expr_dotplot_pvalsortbool_top',opt$topN,'.pdf',sep=''), width = 20, height = 15)
				p1 = ggplot(mydata_final, aes(x=x, y=y)) + geom_point(aes(size=frac, fill=mmean),pch=21) + theme_bw()
				p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax)) + scale_fill_gradientn(colours=c("royalblue", "khaki", "firebrick"), limits=c(fmin,fmax))
				p1 = p1 + scale_y_continuous(breaks=1:nrow(newyorder), labels=paste(row_annot$GO.ID,row_annot$Term,sep=': '), limits = c(0.5,nrow(yorder)+0.5)) + theme(panel.grid.minor.y = element_blank())
				p1 = p1 + scale_x_continuous(breaks=1:nrow(newxorder), labels=newxorder$samp, limits = c(0.5,nrow(xorder)+0.5)) + theme(panel.grid.minor.x = element_blank()) + xlab("Samples") + ylab("GO-terms")
				p1 = p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
				if (opt$GEOmean == FALSE) p1 = p1 + labs(fill = "avg. expr", size = "frac. non-zero")
				if (opt$GEOmean == TRUE) p1 = p1 + labs(fill = "avg. log2(expr)", size = "frac. non-zero")
				show(p1)
				graphics.off()			
			}
		}
		
	} else {
		mydata = mydata_mean_final
		colnames(mydata)[1] = "expr"
		mydata = mydata[order(mydata$expr),,drop=FALSE]
		mydata$x = 1:nrow(mydata)
		row_annot = annot[match(unique(rownames(mydata)),annot$GO.ID),]
		mydata$Term = row_annot$Term
		
		# make plots
		if (opt$GEOmean == TRUE) { ylab = "avg. log2(expr)" } else { ylab = "avg. expr" }
		pdf(paste(opt$outprefix,"_expr_barchart.pdf",sep=''), width = 8, height = max(4,nrow(mydata)/50 + 2))
		a = ggplot(data=mydata, aes(x=x, y=expr))
		a = a + geom_bar(stat="identity", color="black", position = position_dodge()) + coord_flip() 
		a = a + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(ylab) + xlab("GO term") + theme(text = element_text(size=6))
		a = a + scale_x_continuous(breaks=1:nrow(mydata), labels=paste(rownames(mydata),mydata$Term,sep=': '))
		print(a)
		graphics.off()
		
		# repeat, but sorted to match order in p-value plot
		mydata = mydata[match(pvalorder,rownames(mydata)),]
		mydata$x = 1:nrow(mydata)
		row_annot = annot[match(unique(rownames(mydata)),annot$GO.ID),]
		pdf(paste(opt$outprefix,"_expr_barchart_pvalsort.pdf",sep=''), width = 8, height = max(4,nrow(mydata)/50 + 2))
		a = ggplot(data=mydata, aes(x=x, y=expr))
		a = a + geom_bar(stat="identity", color="black", position = position_dodge()) + coord_flip() 
		a = a + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(ylab) + xlab("GO term") + theme(text = element_text(size=6))
		a = a + scale_x_continuous(breaks=1:nrow(mydata), labels=paste(rownames(mydata),mydata$Term,sep=': '))
		print(a)
		graphics.off()	
		
		# repeat, showing only the topN values
		if (maketopN == TRUE) {
			mydata = mydata[rownames(mydata) %in% pvalorder_topN, ]
			mydata = mydata[order(mydata$expr),,drop=FALSE]
			mydata$x = 1:nrow(mydata)
			row_annot = annot[match(unique(rownames(mydata)),annot$GO.ID),]
			pdf(paste(opt$outprefix,"_expr_barchart_top",opt$topN,".pdf",sep=''), width = 8, height = max(4,nrow(mydata)/50 + 2))
			a = ggplot(data=mydata, aes(x=x, y=expr))
			a = a + geom_bar(stat="identity", color="black", position = position_dodge()) + coord_flip() 
			a = a + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(ylab) + xlab("GO term") + theme(text = element_text(size=6))
			a = a + scale_x_continuous(breaks=1:nrow(mydata), labels=paste(rownames(mydata),mydata$Term,sep=': '))
			print(a)
			graphics.off()		
			
			mydata = mydata[match(pvalorder_topN,rownames(mydata)),]
			row_annot = annot[match(unique(rownames(mydata)),annot$GO.ID),]
			mydata$x = 1:nrow(mydata)
			pdf(paste(opt$outprefix,"_expr_barchart_pvalorder_top",opt$topN,".pdf",sep=''), width = 8, height = max(4,nrow(mydata)/50 + 2))
			a = ggplot(data=mydata, aes(x=x, y=expr))
			a = a + geom_bar(stat="identity", color="black", position = position_dodge()) + coord_flip() 
			a = a + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(ylab) + xlab("GO term") + theme(text = element_text(size=6))
			a = a + scale_x_continuous(breaks=1:nrow(mydata), labels=paste(rownames(mydata),mydata$Term,sep=': '))
			print(a)
			graphics.off()		
		}		
	}	
}

# output results
write.table(all_sig_res,file=paste(opt$outprefix,"_results_table.txt",sep=''),col.names = TRUE,row.names = FALSE,quote = F,sep='\t')
cat("DONE, have a nice day!\n")


