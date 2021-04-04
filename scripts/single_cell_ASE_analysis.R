#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(gmodels))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(maxLik))
suppressPackageStartupMessages(library(metap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(vcd))
suppressPackageStartupMessages(library(dplyr))

# version 1.0 (07/31/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 07/31/2019
# -------------------------

# Description:
# -------------------------
# This script is designed to perform imprinted expression analyses given single-cell
# allele-specific count data. 

# The big parts of this analysis (done separately for each gene):
# (1) normalize raw counts according to DEsingle's method
# (2) for each gene, fit the maternal and paternal counts to the optimal distribution (ZINB or NB)
# (3) use outcomes of step 2 to estimate coefficients - relationship between mat. and pat. parameter estimates
# (4) use LRT to test whether we can reject H0: mat and pat counts come from same distribution
#		- also tests 'adjusted' hypothesis based on the coefficients in (3)
# (5) adjust p-values for multiple hypothesis testing and classify each gene into:
#		- no data : too few allelic reads to evaluate
#		- no bias : could not reject H0 (of course, doesn't necessarily mean truly unbiased, just not enough evidence to reject)
# 		- weak mat/pat bias : passes p-value threshold, and log2(m/p) value is more than --log2fc_weak but less than --log2fc_med from median
#		- mat/pat bias : passes p-value threshold, and log2(m/p) value is more than --log2fc_med but less than --log2fc_strong from median
#		- strong mat/pat bias : passes p-value threshold, and log2(m/p) value is more than --log2fc_strong from median

# Required inputs:
# -------------------------
# acounts : total counts for each gene (regardless of allele), in tab-delimited format with rows == genes and columns == cells
# mcounts : maternal counts only, in tab-delimited format with rows == genes and columns == cells
# pcounts : paternal counts only, in tab-delimited format with rows == genes and columns == cells

# Notes:
# -------------------------
# TODO

# Usage: single_cell_ASE_analysis.R [options] acounts mcounts pcounts outprefix

parser = ArgumentParser()
parser$add_argument("acounts", nargs=1, help = "tab-delimited text file matrix of raw, total expression counts for axb cross direction (cols:cells x rows:genes); first col = gene IDs, first row = nuclei IDs")
parser$add_argument("mcounts", nargs=1, help = "tab-delimited text file matrix of raw maternal expression counts for axb cross direction (cols:cells x rows:genes); first col = gene IDs, first row = nuclei IDs")
parser$add_argument("pcounts", nargs=1, help = "tab-delimited text file matrix of raw paternal expression counts for axb cross direction (cols:cells x rows:genes); first col = gene IDs, first row = nuclei IDs")
parser$add_argument("outprefix", nargs=1, help = "prefix for output files")
parser$add_argument("--seed", default = 123456, help = "Set random number generator seed; runs with same seed and other inputs should always produce same results.", type="integer")
parser$add_argument("--minreads", default = 10, help = "Genes with fewer than this total number of allelic reads (maternal + paternal) are censored.", type="integer")
parser$add_argument("--mincellsfrac", default = 0.01, help = "Genes with allelic reads detected in fewer than fraction of cells (maternal + paternal) are censored (for total analysis); for analysis within cluster, at least this fraction of nuclei in cluster must have at least 1 allelic read.", type="double")
parser$add_argument("--minreadspercell", default = 10, help = "Cells with fewer than this many allelic reads across all genes are censored.", type="integer")
parser$add_argument("--nameA", default = "a", help = "Name of mother in cross (default 'a')", type="character")
parser$add_argument("--nameB", default = "b", help = "Name of father in cross (default 'b')", type="character")
parser$add_argument("--coef_nu", help = "Override calculated value of 'coef_nu' with this value", type="double")
parser$add_argument("--coef_mu_ZINB", help = "Override calculated value of 'coef_mu_ZINB' with this value", type="double")
parser$add_argument("--coef_mu_NB", help = "Override calculated value of 'coef_mu_ZINB' with this value", type="double")
parser$add_argument("--skip1", help = "Skip step 1 - fitting ZINB/NB models (must provide 'fits' output file with this argument", type="character")
parser$add_argument("--skip2", help = "Skip step 2 - evaluating ASE (must provide 'results' and 'paramslist' output file with this argument", type="character")
parser$add_argument("--pval", default = 0.05, help = "FDR-adjusted pvalue cutoff for significance", type="double")
parser$add_argument("--log2fc_weak", default = 0.5, help = "log2(m/p) difference from median to use as cutoff for bias - weak bias", type="double")
parser$add_argument("--log2fc_med", default = 2, help = "log2(m/p) difference from median to use as cutoff for bias - med. bias", type="double")
parser$add_argument("--log2fc_strong", default = 5, help = "log2(m/p) difference from median to use as cutoff for bias - strong bias", type="double")
parser$add_argument("--log2fc_median", help = "Script uses deviation from this median value to assess ASE; override actual median with this value (default calculated in-script)", type="double")
parser$add_argument("--mincountsnobias", default = 20, help = "Min counts to allow 'no bias' call; genes with less than this will have 'no bias' replaced with 'no data' in final output", type="double")
parser$add_argument("--nonorm", default=FALSE, action="store_true", help = "Input counts matrix is already normalized; do not normalize again")

opt <- parser$parse_args()

# Check and summarize inputs
# ----------------
set.seed(opt$seed)
mincellsfrac = opt$mincellsfrac
minreads = opt$minreads
minreadspercell = opt$minreadspercell
nameA = opt$nameA
nameB = opt$nameB

pvalcutoff = opt$pval
log2fc_weak = opt$log2fc_weak
log2fc_med = opt$log2fc_med
log2fc_strong = opt$log2fc_strong
mincountsnobias = opt$mincountsnobias

outprefix = opt$outprefix

paramslist = data.frame(parameter=character(),value=character(),stringsAsFactors=FALSE) 
paramslist[1,] = c('mincellsfrac',mincellsfrac)
paramslist[2,] = c('minreads',minreads)
paramslist[3,] = c('minreadspercell',minreadspercell)
paramslist[4,] = c('pval_cutoff',pvalcutoff)
paramslist[5,] = c('log2fc_weak',log2fc_weak)
paramslist[6,] = c('log2fc_med',log2fc_med)
paramslist[7,] = c('log2fc_strong',log2fc_strong)
paramslist[8,] = c('mincountsnobias',mincountsnobias)

args <- commandArgs(trailingOnly = TRUE)	
if (length(opt) <= 1) {
	cat("Usage: single_cell_ASE_analysis.R [options] acounts mcounts pcounts outprefix\n")
	cat("----------------------\n")
	cat("acounts.txt : tab-delimited text file matrix of raw, total expression counts (cols:cells x rows:genes)\n")
	cat("mcounts.txt : tab-delimited text file matrix of raw maternal expression counts (cols:cells x rows:genes)\n")
	cat("pcounts.txt : tab-delimited text file matrix of raw paternal expression counts (cols:cells x rows:genes)\n")
	cat("outprefix : prefix for output files\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("TODO\n")
	cat("----------------------\n")
}

cat("\nRunning single_cell_ASE_analysis.R v.1.0 (08/20/2019)\n")
cat("-----------------------\n")
cat("Current directory:",getwd(),"\n")
cat("total expression matrix:",opt$acounts,"\n")
cat("maternal expression matrix:",opt$mcounts,"\n")
cat("paternal expression matrix:",opt$pcounts,"\n")
cat("Prefix for output files:",outprefix,"\n")
cat("-----------------------\n")
cat("Additional options:\n")
cat("Min. total allelic reads across all cells in cross:",minreads,"\n")
cat("Min. fraction of cells with allelic reads detected for cross:",mincellsfrac,"\n")
cat("Min total allelic reads per cell:",minreadspercell,"\n")
cat("P-value cutoff for significance:",pvalcutoff,"\n")
cat("Min log2fc dif for weak bias:",log2fc_weak,"\n")
cat("Min log2fc dif for med. bias:",log2fc_med,"\n")
cat("Min log2fc dif for strong bias:",log2fc_strong,"\n")
if (! is.null(opt$coef_nu)) {
	cat("Overriding calculated coef_nu with this value:",opt$coef_nu,"\n")
}
if (! is.null(opt$coef_mu_ZINB)) {
	cat("Overriding calculated coef_mu_ZINB with this value:",opt$coef_mu_ZINB,"\n")
}
if (! is.null(opt$coef_mu_NB)) {
	cat("Overriding calculated coef_mu_NB with this value:",opt$coef_mu_NB,"\n")
}
cat("-----------------------\n")

if (! is.null(opt$skip1)) {
	skip1=TRUE
	results = read.table(opt$skip1, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
	cat("Skipping step 1 (fitting ZINB/NB), using this file in place of results:",opt$skip1,"\n")
} else { 
	skip1=FALSE
} 

if (! is.null(opt$skip2)) {
	skip2=TRUE
	filelist = unlist(strsplit(opt$skip2,","))
	if (length(filelist) != 2) { stop("Exactly two files are expected with option --skip2 (see usage info; provide them as a comma-separated list)") }
	results = read.table(filelist[1], header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
	paramslist = read.table(filelist[2], header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)
	cat("Skipping step 2 (evaluating ASE), using these files in place of results:",opt$skip2,"\n")
} else { 
	skip2=FALSE
} 


# FUNCTIONS
# ----------------
# All functions now live in single_cell_ASE_src.R
command_args = commandArgs(trailingOnly = FALSE)
scriptloc = dirname(sub("--file=", "", command_args[grep("--file=", command_args)]))
srcfile <- file.path(scriptloc, "single_cell_ASE_src.R")
source(srcfile)

# MAIN
# ----------------
# read in counts files
acounts = data.matrix(read.table(opt$acounts, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
mcounts = data.matrix(read.table(opt$mcounts, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
pcounts = data.matrix(read.table(opt$pcounts, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))

# check same genes (rownames) included in all matrices
if (setequal(rownames(acounts),rownames(mcounts)) == FALSE) {
	stop("Row names not equal in acounts and mcounts")
} else if (setequal(rownames(acounts),rownames(pcounts)) == FALSE) {
	stop("Row names not equal in acounts and pcounts")
} 

# check same nuclei (colnames) included in all axb, bxa matrices
if (setequal(colnames(acounts),colnames(mcounts)) == FALSE) {
	stop("Column names not equal in acounts and mcounts")
} else if (setequal(colnames(acounts),colnames(pcounts)) == FALSE) {
	stop("Column names not equal in acounts and pcounts")
}

# Order rows (genes) of all matrices by acounts order for consistency
mcounts = mcounts[order(rownames(acounts)), ]
pcounts = pcounts[order(rownames(acounts)), ]

# Check that m+p is always <= a counts
if (! all(mcounts + pcounts <= acounts)) {
	firstexample = which((mcounts + pcounts <= acounts) == FALSE,arr.ind = TRUE)[1,]
	cat("Warning: At least one value in axb 'all' matrix is greater than the maternal + paternal counts\n")
	cat("Example - row",firstexample[1],"and column",firstexample[2],"has:\n")
	cat(mcounts[firstexample[1],firstexample[2]],"maternal counts,",pcounts[firstexample[1],firstexample[2]],"paternal counts and",acounts[firstexample[1],firstexample[2]],"total counts.\n")

	if (length(which((mcounts + pcounts <= acounts) == FALSE,arr.ind = TRUE)) > 1) {
		secondexample = which((mcounts + pcounts <= acounts) == FALSE,arr.ind = TRUE)[2,]
		cat("Example - row",secondexample[1],"and column",secondexample[2],"has:\n")
		cat(mcounts[secondexample[1],secondexample[2]],"maternal counts,",pcounts[secondexample[1],secondexample[2]],"paternal counts and",acounts[secondexample[1],secondexample[2]],"total counts.\n")
	}
	
	tt = nrow(which((mcounts + pcounts <= acounts) == FALSE,arr.ind = TRUE))
	cat("Total occurrences:",tt,"\n")
	if (opt$nonorm == TRUE) {
		cat("This can occur through normalization since values are then converted to int using ceiling().")
	} else {
		stop("Likely error in input matrices.")
	}
}

# Filter input matrices
# ------------------
cat('\n')
cat("Evaluating allelic expression\n")
cat("------------------\n")
cat("Original matrices had",nrow(acounts),"genes and",ncol(acounts),"cells\n")
pass1 = (rowSums(mcounts)+rowSums(pcounts)>=minreads) & ((rowSums(mcounts>0)+rowSums(pcounts>0)) / length(colnames(mcounts)) > mincellsfrac)
pass2 = colSums(mcounts)+colSums(pcounts)>=minreadspercell
acounts_filt = acounts[pass1,pass2]
mcounts_filt = mcounts[pass1,pass2]
pcounts_filt = pcounts[pass1,pass2]
cat("After filtering out genes and cells with low allelic coverage,",nrow(acounts_filt),"genes and",ncol(acounts_filt),"cells remain\n")

# Making sure that input data are integers
acounts_filt <- round(as.matrix(acounts_filt)); mcounts_filt <- round(as.matrix(mcounts_filt)); pcounts_filt <- round(as.matrix(pcounts_filt)); 
storage.mode(acounts_filt) <- "integer"; storage.mode(mcounts_filt) <- "integer"; storage.mode(pcounts_filt) <- "integer"; 

# Normalize all three matrices based on total counts
# ------------------
if (opt$nonorm == FALSE) {
	res = normCounts(acounts_filt,mcounts_filt,pcounts_filt)
	acounts_norm = res[1][[1]]
	mcounts_norm = res[2][[1]]
	pcounts_norm = res[3][[1]]

	# output normalized counts
	write.table(acounts_norm,file=paste(outprefix,"_acounts_norm.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	write.table(mcounts_norm,file=paste(outprefix,"_mcounts_norm.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	write.table(pcounts_norm,file=paste(outprefix,"_pcounts_norm.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
} else {
	acounts_norm = acounts_filt
	mcounts_norm = mcounts_filt
	pcounts_norm = pcounts_filt
}

# Step (0): add pseudocount of 1 to largest maternal and largest paternal count values for each gene
# ------------------
# (ensures we can calculate log2, etc. later with minimal effect on overall distribution)
cat("Adding pseudocount to largest count value in each row (gene):\n")
for (i in 1:nrow(mcounts_norm)) {
	cat("\r",paste0("Processing gene ",i," of ",nrow(mcounts_norm)))
	idx = which.max(mcounts_norm[i,])
	mcounts_norm[i,idx] = mcounts_norm[i,idx] + 1
	idx = which.max(pcounts_norm[i,])
	pcounts_norm[i,idx] = pcounts_norm[i,idx] + 1
}

# Step (1) fit ZINB/NB to maternal and paternal counts separately and estimate params
# ------------------
if (skip1 == FALSE) {

	results <- matrix(data=NA, nrow = nrow(mcounts_norm), ncol = 17, dimnames = list(row.names(mcounts_norm), c("fracnonzero_m", "median_nonzero_m", "mean_m", "fracnonzero_p", "median_nonzero_p", "mean_p", "log2fc_m_over_p", "mu_m", "nu_m", "sigma_m", "fit_m", "logL_m", "mu_p", "nu_p", "sigma_p", "fit_p", "logL_p")))
	results <- as.data.frame(results, stringsAsFactors = FALSE)

	# calculate stats of count distributions (fraction nonzero, median, mean, log2fc)	
	results$fracnonzero_m = rowSums(mcounts_norm != 0) / ncol(mcounts_norm)
	results$median_nonzero_m = apply(mcounts_norm,1,function(x){median(x[x>0])})
	results$mean_m = rowMeans(mcounts_norm)
	results$fracnonzero_p = rowSums(pcounts_norm != 0) / ncol(pcounts_norm)
	results$median_nonzero_p = apply(pcounts_norm,1,function(x){median(x[x>0])})
	results$mean_p = rowMeans(pcounts_norm)
	results$log2fc_m_over_p = log2(results$mean_m / results$mean_p)

	# fit ZINB or NB distributions to each gene
	cat("Fitting best model to each gene:\n")
	for (i in 1:nrow(mcounts_norm)) {
		cat("\r",paste0("Fitting best model to data; processing gene ",i," of ",nrow(mcounts_norm)))

		# fit maternal counts
		res_m = fitdist(mcounts_norm[i,], "mother")
		results$nu_m[i] = res_m[[1]]
		results$mu_m[i] = res_m[[2]]
		results$sigma_m[i] = res_m[[3]]
		results$fit_m[i] = res_m[[4]]
		results$logL_m[i] = res_m[[5]]
		
		# fit paternal counts
		res_p = fitdist(pcounts_norm[i,], "father")
		results$nu_p[i] = res_p[[1]]
		results$mu_p[i] = res_p[[2]]
		results$sigma_p[i] = res_p[[3]]
		results$fit_p[i] = res_p[[4]]
		results$logL_p[i] = res_p[[5]]
	}
	cat("DONE\n")

	# done fitting model, output just the fits
	write.table(results,file=paste(outprefix,"_fits.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
} else {
	# check dimensions of 'results' file provided matches length of acounts_norm
	if (dim(results)[1] != dim(acounts_norm)[1]) { stop("File provided with --skip1 has",dim(results)[1],"rows while filtered count matrices have",dim(acounts_norm)[1],", these must be same (are you sure that you used the right 'results' file?)") }
}

# make plots of the nu, mu and sigma parameters in mat vs. pat (ZINB and NB fits separately)
# ------------------

# compare the maternal and paternal parameter estimates for genes fit to ZINB and NB separately
results_ZINB = results[results$fit_m == "ZINB" & results$fit_p == "ZINB" & results$mean_m > min(results$mean_m) & results$mean_p > min(results$mean_p),]
results_NB = results[results$fit_m == "NB" & results$fit_p == "NB" & results$mean_m > min(results$mean_m) & results$mean_p > min(results$mean_p),]

# nu parameter (actually do 1-nu since that's easier to interpret)
results_ZINB$invnu_m = 1 - results_ZINB$nu_m
results_ZINB$invnu_p = 1 - results_ZINB$nu_p
reg_nu = lm(results_ZINB$invnu_p ~ 0 + results_ZINB$invnu_m)
pdf(paste(outprefix,"_nu_fit_ZINB.pdf",sep=''), width = 8, height = 8)
smoothScatter(results_ZINB$invnu_m, results_ZINB$invnu_p, xlab='Estimated (1-nu) - maternal counts', ylab = 'Estimated (1-nu) - paternal counts', xlim = c(0,1), ylim=c(0,1))
abline(0,reg_nu$coefficients[1],col='black')
abline(0,1,col='gray60',lty=3)
eq <- paste0("(1-nu_p) = ",ifelse(sign(reg_nu$coefficients[1])==1, " + ", " - "), format(abs(reg_nu$coefficients[1]),digits=3), " (1-nu_m); R^2(adj) =",format(summary(reg_nu)$adj.r.squared,digits=3))
mtext(eq, 3, line=1)
graphics.off()

# mu parameter (ZINB fits only)
lim_mu = min(max(results_ZINB$mu_p),max(results_ZINB$mu_m))
reg_mu = lm(results_ZINB$mu_p ~ 0+results_ZINB$mu_m)
pdf(paste(outprefix,"_mu_fit_ZINB.pdf",sep=''), width = 8, height = 8)
smoothScatter(results_ZINB$mu_m, results_ZINB$mu_p, xlab='Estimated mu - maternal counts', ylab = 'Estimated mu - paternal counts', xlim=c(0,lim_mu),ylim=c(0,lim_mu))
abline(0,reg_mu$coefficients[1],col='black')
abline(0,1,col='gray60',lty=3)
eq <- paste0("mu_p = ",ifelse(sign(reg_mu$coefficients[1])==1, " + ", " - "), format(abs(reg_mu$coefficients[1]),digits=3), " mu_m; R^2(adj) =",format(summary(reg_mu)$adj.r.squared,digits=3))
mtext(eq, 3, line=1)
graphics.off()

# sigma parameter (ZINB fits only)
lim_sigma = min(max(results_ZINB$sigma_p),max(results_ZINB$sigma_m))
reg_sigma = lm(results_ZINB$sigma_p ~ 0+results_ZINB$sigma_m)
pdf(paste(outprefix,"_sigma_fit_ZINB.pdf",sep=''), width = 8, height = 8)
smoothScatter(results_ZINB$sigma_m, results_ZINB$sigma_p, xlab='Estimated sigma - maternal counts', ylab = 'Estimated sigma - paternal counts', xlim=c(0,lim_sigma),ylim=c(0,lim_sigma))
abline(0,reg_sigma$coefficients[1],col='black')
abline(0,1,col='gray60',lty=3)
eq <- paste0("sigma_p = ",ifelse(sign(reg_sigma$coefficients[1])==1, " + ", " - "), format(abs(reg_sigma$coefficients[1]),digits=3), " sigma_m; R^2(adj) =",format(summary(reg_sigma)$adj.r.squared,digits=3))
mtext(eq, 3, line=1)
graphics.off()

# mu parameter (NB fits only) - note, this is given to extreme high values so only fit lower 90th percentile of distn'
lim_mu = max(quantile(results_NB$mu_m, c(.9)),quantile(results_NB$mu_p, c(.9)))
reg_mu_NB = lm(results_NB[results_NB$mu_m < lim_mu & results_NB$mu_p < lim_mu,]$mu_p ~ 0+results_NB[results_NB$mu_m < lim_mu & results_NB$mu_p < lim_mu,]$mu_m)
pdf(paste(outprefix,"_mu_fit_NB.pdf",sep=''), width = 8, height = 8)
smoothScatter(results_NB$mu_m, results_NB$mu_p, xlab='Estimated mu - maternal counts', ylab = 'Estimated mu - paternal counts', xlim=c(0,lim_mu),ylim=c(0,lim_mu))
abline(0,reg_mu_NB$coefficients[1],col='black')
abline(0,1,col='gray60',lty=3)
eq <- paste0("mu_p = ",ifelse(sign(reg_mu_NB$coefficients[1])==1, " + ", " - "), format(abs(reg_mu_NB$coefficients[1]),digits=3), " mu_m; R^2(adj) =",format(summary(reg_mu_NB)$adj.r.squared,digits=3))
mtext(eq, 3, line=1)
graphics.off()

# sigma parameter (NB fits only)
lim_sigma = min(max(results_NB$sigma_p),max(results_NB$sigma_m))
reg_sigma = lm(results_NB$sigma_p ~ 0+results_NB$sigma_m)
pdf(paste(outprefix,"_sigma_fit_NB.pdf",sep=''), width = 8, height = 8)
smoothScatter(results_NB$sigma_m, results_NB$sigma_p, xlab='Estimated sigma - maternal counts', ylab = 'Estimated sigma - paternal counts', xlim=c(0,lim_sigma),ylim=c(0,lim_sigma))
abline(0,reg_sigma$coefficients[1],col='black')
abline(0,1,col='gray60',lty=3)
eq <- paste0("sigma_p = ",ifelse(sign(reg_sigma$coefficients[1])==1, " + ", " - "), format(abs(reg_sigma$coefficients[1]),digits=3), " sigma_m; R^2(adj) =",format(summary(reg_sigma)$adj.r.squared,digits=3))
mtext(eq, 3, line=1)
graphics.off()

# gather resulting predicted differences between mu and nu mat vs. pat
coef_nu = reg_nu$coefficients[1]; names(coef_nu) = NULL
coef_mu_ZINB = reg_mu$coefficients[1]; names(coef_mu_ZINB) = NULL
coef_mu_NB = reg_mu_NB$coefficients[1]; names(coef_mu_NB) = NULL

# if provided by user, override these
if (! is.null(opt$coef_nu)) {
	coef_nu = opt$coef_nu
}
if (! is.null(opt$coef_mu_ZINB)) {
	coef_mu_ZINB = opt$coef_mu_ZINB
}
if (! is.null(opt$coef_mu_NB)) {
	coef_mu_NB = opt$coef_mu_NB
}

# save coefs to paramslist
paramslist[9,] = c('coef_nu',coef_nu)
paramslist[10,] = c('coef_mu_ZINB',coef_mu_ZINB)
paramslist[11,] = c('coef_mu_NB',coef_mu_NB)


# Step (2) use estimated values from (1) to do LRT
# ------------------
if (skip2 == FALSE) {
	results$fit_H0 = rep(NA, nrow(mcounts_norm))
	results$nu_H0 = rep(NA, nrow(mcounts_norm))
	results$mu_H0 = rep(NA, nrow(mcounts_norm))
	results$sigma_H0 = rep(NA, nrow(mcounts_norm))
	results$logL_H1 = rep(NA, nrow(mcounts_norm))
	results$logL_H0 = rep(NA, nrow(mcounts_norm))
	results$pval = rep(NA, nrow(mcounts_norm))
	results$nu_H0_mod = rep(NA, nrow(mcounts_norm))
	results$mu_H0_mod = rep(NA, nrow(mcounts_norm))
	results$sigma_H0_mod = rep(NA, nrow(mcounts_norm))
	results$logL_H0_mod = rep(NA, nrow(mcounts_norm))
	results$pval_mod = rep(NA, nrow(mcounts_norm))
	results$logL_H0_2= rep(NA, nrow(mcounts_norm))
	results$pval_H0_2= rep(NA, nrow(mcounts_norm))
	results$logL_H0_2_mod= rep(NA, nrow(mcounts_norm))
	results$pval_H0_2_mod= rep(NA, nrow(mcounts_norm))
	results$logL_H0_3= rep(NA, nrow(mcounts_norm))
	results$pval_H0_3= rep(NA, nrow(mcounts_norm))
	results$logL_H0_3_mod= rep(NA, nrow(mcounts_norm))
	results$pval_H0_3_mod= rep(NA, nrow(mcounts_norm))

	for (i in 1:nrow(mcounts_norm)) {	
		cat("\r",paste0("Evaluating ASE; processing gene ",i," of ",nrow(mcounts_norm)))
		res = testH0(mcounts_norm[i,], pcounts_norm[i,], results[i,], coef_nu = coef_nu, coef_mu_ZINB = coef_mu_ZINB, coef_mu_NB = coef_mu_NB)
		results[i,] = res
	}

	cat("DONE\n")
}

# Step (3) classify genes into biased and unbiased based on outcome of LRT
# ------------------

# correct p-values using Benjamini-Hochberg method
results$padj = p.adjust(results$pval, method="fdr")
results$padj_mod = p.adjust(results$pval_mod, method="fdr")
results$padj_H0_2 = p.adjust(results$pval_H0_2, method="fdr")
results$padj_H0_2_mod = p.adjust(results$pval_H0_2_mod, method="fdr")
results$padj_H0_3 = p.adjust(results$pval_H0_3, method="fdr")
results$padj_H0_3_mod = p.adjust(results$pval_H0_3_mod, method="fdr")

# classify genes into categories
if (! is.null(opt$log2fc_median)) {
	log2fcmedval = opt$log2fc_median
} else {
	log2fcmedval = median(results[results$mean_m + results$mean_p > 0.1,]$log2fc_m_over_p)
}

fracnodatathreshold = mincountsnobias / ncol(mcounts_norm)

paramslist[12,] = c('log2fcmedval',log2fcmedval)
paramslist[13,] = c('fracnodatathreshold',fracnodatathreshold)

results$status = ifelse(results$padj_mod < pvalcutoff & results$log2fc_m_over_p > log2fcmedval + log2fc_weak, "weak mat bias", "TBD")
results$status = ifelse(results$padj_mod < pvalcutoff & results$log2fc_m_over_p < log2fcmedval - log2fc_weak, "weak pat bias", results$status)
results$status = ifelse(results$status == "TBD", "no bias", results$status)
results$status = ifelse(results$status == "no bias" & results$mean_m + results$mean_p < fracnodatathreshold, "no data", results$status)

results$status = ifelse(results$status == "weak mat bias" & results$padj_mod < pvalcutoff & results$log2fc_m_over_p > log2fcmedval + log2fc_med, "mat bias", results$status)
results$status = ifelse(results$status == "weak pat bias" & results$padj_mod < pvalcutoff & results$log2fc_m_over_p < log2fcmedval - log2fc_med, "pat bias", results$status)

results$status = ifelse(results$status == "mat bias" & results$padj_mod < pvalcutoff & results$log2fc_m_over_p > log2fcmedval + log2fc_strong, "strong mat bias", results$status)
results$status = ifelse(results$status == "pat bias" & results$padj_mod < pvalcutoff & results$log2fc_m_over_p < log2fcmedval - log2fc_strong, "strong pat bias", results$status)

results$status = ifelse(results$status == "no bias" & results$pval_mod < pvalcutoff & results$log2fc_m_over_p > log2fcmedval + log2fc_med, "pot. mat bias", results$status)
results$status = ifelse(results$status == "no bias" & results$pval_mod < pvalcutoff & results$log2fc_m_over_p < log2fcmedval - log2fc_med, "pot. pat bias", results$status)	

# also evaluate the secondary hypotheses H02 and H03
results$btype = ifelse(results$status == "no bias" | results$status == "no data" | results$status == "pot. mat bias" | results$status == "pot. pat bias", "N/A", "TBD")
results$btype = ifelse(results$btype == "TBD" & results$padj_H0_2_mod < pvalcutoff & results$padj_H0_3_mod >= pvalcutoff, "mag", results$btype)
results$btype = ifelse(results$btype == "TBD" & results$padj_H0_3_mod < pvalcutoff & results$padj_H0_2_mod >= pvalcutoff, "freq", results$btype)
results$btype = ifelse(results$btype == "TBD", "all", results$btype)

classified = results[,colnames(results) %in% c('status','btype')]

write.table(results,file=paste(outprefix,"_results.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
write.table(paramslist,file=paste(outprefix,"_paramslist.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
write.table(classified,file=paste(outprefix,"_classified.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')


