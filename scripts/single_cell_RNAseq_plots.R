#!/usr/bin/env Rscript

# version 1.0 (09/07/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 09/07/2019
# v.2.0: 12/03/2019
#	- added 'nuc' plot type
#	- fleshed out description section
# v.3.0: 12/09/2019
#	- added 'dot' plot type
# 	- changed input setup, able to accept more inputs (not just for ASE analyses but for other plot types
#	relevant to single cell analysis)
#	- added info that prints when no args provided
#	- renamed to single_cell_RNAseq_plots.R (from single_cell_ASE_make_plots.R, which is name of v1 and v2)
# -------------------------

# Description:
# -------------------------
# This simple script can make several useful plots for analyzing single-cell count data, including 
# allele-specific count data. The main arguments plottype and outfile are required for all plot types.

# Plot type 1: 'gof'
# Goodness-of-fit plot that compares the fit of H1 (separate mat and pat fits) and H0 (joint fit)
# distribution, for a single gene. 
# Required inputs:
#	--gene : which (single) gene to make plot for
#	--coef_nu : nu coefficient from main ASE analysis output (see single_cell_ASE_analysis.R)
#	--coef_mu_ZINB : mu coefficient for ZINB dist, from main ASE analysis output (see single_cell_ASE_analysis.R)
#	--coef_mu_NB : mu coefficient for NB dist, from main ASE analysis output (see single_cell_ASE_analysis.R)
# Optional params:
#	--xmax : maximum value to plot on x-axis

# Plot type 2: 'cmp'
# Goodness-of-fit plot that compares the ZINB and NB best fits, separately for both mat and pat counts
# Required inputs:
#	--gene : which (single) gene to make plot for
# Optional params:
#	--xmax : maximum value to plot on x-axis

# Plot type 3: 'bar'
# Plots total maternal and paternal counts per cell/nucleus, on the same plot as a bar chart, for
# a single gene, OR total counts (not allele-specific).
# Required inputs:
#	--gene : which (single) gene to make plot for
# Optional params:
#	--clusters : sorts nuclei (x-axis) by clustering, and adds vertical lines demarcating boundaries between nuclei clusters

# Plot type 4: 'nuc'
# Similar to 'bar' in terms of info content, but more easily shows multiple genes instead of just one.
# Each gene is represented by a row of dots, corresponding to cells/nuclei. Each dot's size is a function of 
# number of allelic reads are present, color is a function of % maternal. 
# Required inputs:
# 	--gene : list of genes (comma-delimited, no spaces) to make a plot for
# Optional params:
#	--clusters : splits nuclei according to clusters

# Plot type 5: 'dot'
# If --acounts or --CPM provided: dot plot over genes/gene groups vs. nuclei/sample groups, where dot color is expression 
# level and dot size is the fraction of informative cells/nuclei (averaged over sample & gene groups, where applicable). \
# If --mcounts and --pcounts provided, dot color is average % maternal, and dot size is the fraction of nuclei with allelic reads
# (again averaged over sample & gene groups, where applicable). Note that --CPM and --acounts are actually treated exactly
# the same (counts are assumed normalized e.g. by the method used in DESingle), the only difference is in the figure legend,
# which will correctly refer to CPM or normalized counts respectively.
# Required inputs:
# 	Either --acounts OR --CPM : a matrix (rows: genes x cols: nuclei) 
# Optional params:
#	--clusters : splits nuclei according to clusters
#	--gene OR --genegroups : combine multiple genes into these groups

# Plot type 6: 'lin'
# Plots average total (expr. mode) or mat and pat (imprinting mode) expression, as points, over samples. If clusters 
# provided, plots average over all nuclei in clusters instead, and also plots # informative nuclei per cluster. Points
# are connected by lines to show general trends.
# Required inputs:
# 	Either --acounts OR --CPM : a matrix (rows: genes x cols: nuclei) 
#	--gene : which (single) gene to make plot for
# Optional params:
#	--clusters : sorts nuclei (x-axis) by clustering, and adds vertical lines demarcating boundaries between nuclei clusters

# List of required libraries:
# -------------------------
liblist = c('ggplot2','gridExtra','vcd','maxLik','gamlss','reshape2')
	
# Notes:
# -------------------------
# (1) Output file will always be PDF

# (2) Input files should contain -normalized counts-. That is, counts normalized via the normCounts()
# 	function that currently lives in single_cell_ASE_src.R. Otherwise, your mileage may vary.

# Usage: single_cell_ASE_analysis.R [options] expr_matrix.txt outprefix

# -------------------------

# Check that all required libraries are installed
if(length(setdiff(liblist, rownames(installed.packages()))) > 0) {
	missinglibs = setdiff(liblist, rownames(installed.packages()))
	missinglibs = paste(as.character(missinglibs), sep="' '", collapse=", ")
	stop("the following required libraries are not installed: ",missinglibs,". Please install them, then re-run this script.")
}

# Load argument parser
suppressPackageStartupMessages(library(argparse))
	
# Read in user-supplied arguments
parser = ArgumentParser()

# Required arguments
parser$add_argument("plottype", nargs=1, help = "one of the plot types listed in file description (see above)")
parser$add_argument("outprefix", nargs=1, help = "prefix for output file(s); all plots will output as ${outprefix}_plot.pdf, data files as ${outprefix}_data.txt")

# Other arguments
parser$add_argument("--seed", default = 123456, nargs=1, help = "Seed from random number generator (ensures consistent output).")
parser$add_argument("--acounts", nargs=1, help = "Either this or --CPM can be used for 'dot' plot; tab-delimited text file matrix of filtered, normalized maternal expression counts (cols:cells x rows:genes); see output of single_cell_ASE_analysis.R")
parser$add_argument("--mcounts", nargs=1, help = "Required for 'gof', 'bar' and 'nuc' plots, can be used for 'dot' plot; tab-delimited text file matrix of filtered, normalized maternal expression counts (cols:cells x rows:genes); see output of single_cell_ASE_analysis.R")
parser$add_argument("--pcounts", nargs=1, help = "Required for 'gof', 'bar' and 'nuc' plots, can be used for 'dot' plot; tab-delimited text file matrix of filtered, normalized paternal expression counts (cols:cells x rows:genes); see output of single_cell_ASE_analysis.R")
parser$add_argument("--CPM", nargs=1, help = "Either this or --acounts can be used for 'dot' plot; tab-delimited text file matrix of CPM values (cols:cells x rows:genes)")
parser$add_argument("--genes", help = "Single gene or comma-delimited list of genes to use in plot (must be 1 gene for 'gof' or 'bar', can be multiple for 'nuc' or 'dot'", type="character")
parser$add_argument("--genefile", help = "NO HEADER! Can use in place of --gene to provide a file with a list of genes for 'nuc' or 'dot' plots; if this file has two columns, second column assumed to provide 'gene groups', and if plottype == 'dot', all genes in the same group will be averaged together to create one row in plot per group, instead of per gene. 2nd column ignored by all other plot types.", type="character")
parser$add_argument("--sampfile", help = "NO HEADER! (optional) - list of samples (columns in input matrices) to use in plot, all other samples in matrix omitted; if this file has two columns, second column assumed to provide 'sample groups' (e.g. clusters), and if plottype == 'dot', all genes in the same group will be averaged together to create one row in plot per group, instead of per gene.", type="character")
parser$add_argument("--ASE_params", help = "(required only if plottype = 'gof') - *_paramslist.txt file output from a run of single_cell_ASE_analysis.R, make sure to use mcounts/pcounts from same run", type="character")
parser$add_argument("--coef_nu", default=1, help = "(required only if plottype = 'gof' or 'lin'; ignored if --ASE_params supplied instead) - coef_nu from single_cell_ASE_analysis.R, make sure to use mcounts/pcounts from same run", type="double")
parser$add_argument("--coef_mu_ZINB", default=1, help = "(required only if plottype = 'gof' or 'lin'; ignored if --ASE_params supplied instead) - coef_mu_ZINB from single_cell_ASE_analysis.R, make sure to use mcounts/pcounts from same run", type="double")
parser$add_argument("--coef_mu_NB", default=1, help = "(required only if plottype = 'gof' or 'lin'; ignored if --ASE_params supplied instead) - coef_mu_NB from single_cell_ASE_analysis.R, make sure to use mcounts/pcounts from same run", type="double")
parser$add_argument("--xmax", help = "(optional, only used if plottype = 'gof') set max value shown on x-axis (use to truncate if tail is very long)", type="double")
parser$add_argument("--minreads", default=1, help = "(optional, only used if plottype = 'nuc', or 'dot' for imprinting data) - filter out nuclei with fewer than this many total allelic reads. For plottype 'lin', omits points for nuclei/clusters with fewer than this many total allelic reads.", type="integer")
parser$add_argument("--maxdims", default=100, help = "(optional, only used if plottype = 'dot') - creates an error if requested dot plot would have more than either 100 rows or columns (would be very hard to interpret this plot)", type="integer")
parser$add_argument("--ymin", help = "(optional, only used if plottype = 'lin') - min value to plot on y-axis", type="double")
parser$add_argument("--ymax", help = "(optional, only used if plottype = 'lin') - max value to plot on y-axis", type="double")
parser$add_argument("--xorder", help = "(optional, only used if plottype = 'dot' or 'lin') - order in which to plot samples/sample clusters (see --sampfile) on x-axis, from left-to-right", type="character")
parser$add_argument("--yorder", help = "(optional, only used if plottype = 'dot') - order in which to plot genes/gene groups (see --genefile) on y-axis, from top-to-bottom", type="character")
parser$add_argument("--dotsize", default = 10, help = "(optional, only used if plottype = 'dot') - upper limit for size of points in dot plot (default 10)", type="double")
parser$add_argument("--linewidth", default = 2, help = "(optional, only used if plottype = 'lin') - size of line in line plot", type="double")
parser$add_argument("--fillupper", help = "(optional, only used if plottype = 'dot' & --mcounts/pcounts provided) - upper limit of value for color scale (default automatically set to highest value in plot)", type="double")
parser$add_argument("--sizeupper", help = "(optional, only used if plottype = 'dot' & --mcounts/pcounts provided) - upper limit of value for size scale (default automatically set to highest value in plot)", type="double")
parser$add_argument("--samesize", default=FALSE, action="store_true", help = "Make all points in plot same size (only used in 'dot' plots).")
parser$add_argument("--outputdatalong", default=FALSE, action="store_true", help = "Output the final set of processed data used to output the plot (long format).")
parser$add_argument("--outputdatawide", default=FALSE, action="store_true", help = "Output the final set of processed data used to output the plot (wide format).")
parser$add_argument("--skipplot", default=FALSE, action="store_true", help = "Do not output plot (if used without --outputdata, script will output nothing).")
parser$add_argument("--allowmissinggenes", default=FALSE, action="store_true", help = "If using --genefile, allow some of those genes to be absent from input matrices")
parser$add_argument("--noprelog", default=FALSE, action="store_true", help = "(only used if plottype == 'dot' or 'lin') do NOT take log2 of count/CPM values BEFORE calculating averages")
parser$add_argument("--postlog", default=FALSE, action="store_true", help = "(only used if plottype == 'dot' or 'lin') DO take log2 of count/CPM values AFTER calculating averages (default no, only applies if --noprelog is on)")
parser$add_argument("--includezeros", default=FALSE, action="store_true", help = "(only used if plottype == 'dot' or 'lin') - include zero values when calculating mean expression (normally zeros are excluded from mean calculation, represented instead by frac nonmissing = dot size))")
parser$add_argument("--addsims", default=FALSE, action="store_true", help = "(only used if plottype == 'lin' and --mcounts/pcounts provided) - include in plot a line simulating paternal reads if drawn from maternal distribution")
parser$add_argument("--nreps", default=500, help = "(optional, only used if plottype = 'lin', --mcounts/pcounts provided, and --addsims used) - number of times to simulate random draws from NB/ZINB to simulate doubled paternal dosage", type="double")
parser$add_argument("--plotheight", help = "(optional) - height (in inches) of plot", type="double")
parser$add_argument("--plotwidth", help = "(optional) - width (in inches) of plot", type="double")

opt <- parser$parse_args()
set.seed(opt$seed)

# FUNCTIONS
# ---------------------

# All functions related to ASE live in single_cell_ASE_src.R
command_args = commandArgs(trailingOnly = FALSE)
scriptloc = dirname(sub("--file=", "", command_args[grep("--file=", command_args)]))
srcfile <- file.path(scriptloc, "single_cell_ASE_src.R")
source(srcfile)

# Plot comparing goodness-of-fit of H1 (separate fit for m and p counts) vs. H0 (joint fit)
make_gof_plot <- function(mcounts, pcounts, outprefix, coef_nu = 1, coef_mu_ZINB = 1, coef_mu_NB = 1, xmax = NULL) {

	# convert to vectors
	mcounts = as.numeric(mcounts)
	pcounts = as.numeric(pcounts)

	# get true count distributions
	m_obs <- table(mcounts)
	m_max <- max(as.numeric(names(m_obs)))

	p_obs <- table(pcounts)
	p_max <- max(as.numeric(names(p_obs)))

	maxv = max(m_max, p_max)
	rr = 0:maxv

	m_robs = rr[!(rr %in% names(m_obs))]
	m_obs2 = rep(0, length(m_robs))
	names(m_obs2) = m_robs
	m_obs = c(m_obs,m_obs2)
	m_obs = m_obs[order(as.numeric(names(m_obs)))]

	p_robs = rr[!(rr %in% names(p_obs))]
	p_obs2 = rep(0, length(p_robs))
	names(p_obs2) = p_robs
	p_obs = c(p_obs,p_obs2)
	p_obs = p_obs[order(as.numeric(names(p_obs)))]
	
	fracnonzero_m = 1 - (m_obs[1] / sum(m_obs))
	fracnonzero_p = 1 - (p_obs[1] / sum(p_obs))
	
	# get best fit for both
	m_fit = fitdist(mcounts,"mother")
	p_fit = fitdist(pcounts,"father")
	
	# get fitted distributions for both mat and pat, and predicted distributions based on other parents' fit
	if (m_fit[[4]] == "ZINB") {
		m_fit_actual = getProbZinbConv(as.numeric(names(m_obs)), m_fit[[1]], m_fit[[2]], m_fit[[3]])*sum(m_obs)
	} else {
		m_fit_actual = dNBI(as.numeric(names(m_obs)), mu = m_fit[[2]]*2, sigma = m_fit[[3]]/2)*sum(m_obs)
	}
	
	if (p_fit[[4]] == "ZINB") {
		p_fit_actual = dZINBI(as.numeric(names(p_obs)), nu = p_fit[[1]],mu = p_fit[[2]], sigma = p_fit[[3]])*sum(m_obs)
	} else {
		p_fit_actual = dNBI(as.numeric(names(m_obs)), mu = p_fit[[2]], sigma = p_fit[[3]])*sum(m_obs)
	}
	
	# get fits under joint (restricted) model
	if (m_fit[[4]] == "ZINB" && p_fit[[4]] == "ZINB") {
		res = fitZINBboth(mcounts,pcounts)
		fit_H0 = "ZINB"
	} else if (m_fit[[4]] == "NB" && p_fit[[4]] == "NB") {
		res = fitNBboth(mcounts,pcounts)
		fit_H0 = "NB"
	} else {
		res1 = fitZINBboth(mcounts,pcounts)
		res2 = fitNBboth(mcounts,pcounts)
		if (res1[1] > res2[1]) {
			res = res1
			fit_H0 = "ZINB"
		} else {
			res = res2
			fit_H0 = "NB"
		}
	}

	if (fit_H0 == "ZINB") {
		pred_m_fit = c(res[2], res[3], res[4])
		pred_p_fit = c(max(0.000001,min(1-0.00001,1-((1-res[2])*coef_nu))), max(0.000001,res[3]*coef_mu_ZINB), res[4])
		m_fit_pred = getProbZinbConv(as.numeric(names(m_obs)), pred_m_fit[1], pred_m_fit[2], pred_m_fit[3])*sum(m_obs)	
		p_fit_pred = dZINBI(as.numeric(names(p_obs)), nu = pred_p_fit[1],mu = pred_p_fit[2], sigma = pred_p_fit[3])*sum(m_obs)
		logL_H0 = res[1]
	} else {
		pred_m_fit = c(NA, res[2], res[3])
		pred_p_fit = c(NA, max(0.000001,res[2]*coef_mu_NB), res[3])
		m_fit_pred = dNBI(as.numeric(names(m_obs)), mu = pred_m_fit[2]*2, sigma = pred_m_fit[3]/2)*sum(m_obs)
		p_fit_pred = dNBI(as.numeric(names(m_obs)), mu = pred_p_fit[2], sigma = pred_p_fit[3])*sum(m_obs)
		logL_H0 = res[1]
	}		
			
	# build plot (top = ZINB, bot = NB)
	if (opt$skipplot == FALSE) {
		pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 10)
		pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 12)
		cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
		pdf(paste(outprefix,'_plot.pdf',sep=''), width = pwidth, height = pheight)
	}

	# top (true best fit)
	maxyLHS = max(fracnonzero_m,fracnonzero_p,1-(m_fit_pred[1]/sum(m_fit_pred)),1-(m_fit_actual[1]/sum(m_fit_actual)),1-(p_fit_pred[1]/sum(p_fit_pred)),1-(p_fit_actual[1]/sum(p_fit_actual)))
	maxyRHS = max(m_obs[2:length(m_obs)],p_obs[2:length(p_obs)],m_fit_actual[2:length(m_fit_actual)],m_fit_pred[2:length(m_fit_pred)],p_fit_actual[2:length(p_fit_actual)],p_fit_pred[2:length(p_fit_pred)])
	minyRHS = -1*maxyRHS

	mydf1 <- data.frame(vals = c(fracnonzero_m,fracnonzero_p), names=c("Mat","Pat"), expected = c(1-(m_fit_actual[1]/sum(m_fit_actual)),1-(p_fit_actual[1]/sum(p_fit_actual))))
	m_obs_nz = m_obs[names(m_obs) != '0']
	p_obs_nz = p_obs[names(p_obs) != '0']
	mydata1 = data.frame(group = c(rep("Mat",maxv),rep("Pat",maxv)), x = rep(as.numeric(names(m_obs_nz)), 2), y = c(m_obs_nz,p_obs_nz*-1), fit = c(m_fit_actual[2:length(m_fit_actual)],-1*p_fit_actual[2:length(p_fit_actual)]))
	if (opt$skipplot == FALSE) {
		p1 = ggplot(data=mydf1, aes(x=names, y=vals, fill=names)) + geom_bar(stat="identity") + ylab("Fraction of nonzero counts") + xlab("") + geom_point(aes(x=names, y=expected))+ theme_bw() + ylim(0,maxyLHS)
		p2 = ggplot(mydata1, aes(x=x, y=y, fill=group)) + geom_bar(stat="identity", position="identity") + xlab("Nonzero read counts") + ylab("Cells with indicated read count") + geom_line(aes(x=x, y=fit))+ theme_bw() + ylim(minyRHS,maxyRHS)

		tg1 = textGrob(paste("Mat: (1-nu)=",round(1-m_fit[[1]],4),'\n',"Pat: (1-nu)=",round(1-p_fit[[1]],4),sep=''), gp=gpar(fontsize=14), just=c("center", "center"))
		tg2 = textGrob(paste("Mat: mu=",round(m_fit[[2]],4),'; sigma=',round(m_fit[[3]],4),'\n',"Pat: mu=",round(p_fit[[2]],4),'; sigma=',round(p_fit[[3]],4),sep=''), gp=gpar(fontsize=14), just=c("center", "center"))
		tg3 = textGrob(paste("Mat fit: ",m_fit[[4]],"; Pat fit: ",p_fit[[4]],"\nH1 logL: ",round(m_fit[[5]]+p_fit[[5]],4),sep=''), gp=gpar(fontsize=14), just=c("center", "center"))
	}

	mydf2 <- data.frame(vals = c(fracnonzero_m,fracnonzero_p), names=c("Mat","Pat"), expected = c(1-(m_fit_pred[1]/sum(m_fit_pred)),1-(p_fit_pred[1]/sum(p_fit_pred))))
	mydata2 = data.frame(group = c(rep("Mat",maxv),rep("Pat",maxv)), x = rep(as.numeric(names(m_obs_nz)), 2), y = c(m_obs_nz,p_obs_nz*-1), fit = c(m_fit_pred[2:length(m_fit_pred)],-1*p_fit_pred[2:length(p_fit_pred)]))
	if (opt$skipplot == FALSE) {
		p3 = ggplot(data=mydf2, aes(x=names, y=vals, fill=names)) + geom_bar(stat="identity") + ylab("Fraction of nonzero counts") + xlab("") + geom_point(aes(x=names, y=expected))+ theme_bw() + ylim(0,maxyLHS)
		p4 = ggplot(mydata2, aes(x=x, y=y, fill=group)) + geom_bar(stat="identity", position="identity") + xlab("Nonzero read counts") + ylab("Cells with indicated read count") + geom_line(aes(x=x, y=fit))+ theme_bw() + ylim(minyRHS,maxyRHS)

		tg4 = textGrob(paste("Joint: (1-nu)=",round(1-pred_m_fit[1],4),'\n',sep=''), gp=gpar(fontsize=14), just=c("center", "top"))
		tg5 = textGrob(paste("Joint: mu=",round(pred_m_fit[2],4),'; sigma=',round(pred_m_fit[3],4),'\n',sep=''), gp=gpar(fontsize=14), just=c("center", "top"))
		tg6 = textGrob(paste("Joint fit: ",fit_H0,"\nH0 logL: ",round(logL_H0,4),sep=''), gp=gpar(fontsize=14), just=c("center", "center"))

		if (! is.null(xmax)) {
			p2 = p2 + xlim(0,xmax)
			p4 = p4 + xlim(0,xmax)
		}

		lay = as.matrix(rbind(c(1,2,3),c(4,5,5),c(6,7,8),c(9,10,10)))
		grid.arrange(tg1,tg2,tg3,p1,p2,tg4,tg5,tg6,p3,p4,layout_matrix = lay, heights=c(0.1, 0.85, 0.1, 0.85))
		graphics.off()
	}
	
	if (opt$outputdatalong == TRUE || opt$outputdatawide == TRUE) {
		write.table(mydf1,file=paste(outprefix,"_H1_zero.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
		write.table(mydata1,file=paste(outprefix,"_H1_nonzero.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
		write.table(mydf2,file=paste(outprefix,"_H0_zero.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
		write.table(mydata2,file=paste(outprefix,"_H0_nonzero.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	}
}

# Plot comparing goodness-of-fit of ZINB (top) vs. NB (bottom) for a single gene
make_cmp_plot <- function(mcounts, pcounts, outprefix, xmax = NULL) {

	# convert to vectors
	mcounts = as.numeric(mcounts)
	pcounts = as.numeric(pcounts)

	# get true count distributions
	m_obs <- table(mcounts)
	m_max <- max(as.numeric(names(m_obs)))

	p_obs <- table(pcounts)
	p_max <- max(as.numeric(names(p_obs)))

	maxv = max(m_max, p_max)
	rr = 0:maxv

	m_robs = rr[!(rr %in% names(m_obs))]
	m_obs2 = rep(0, length(m_robs))
	names(m_obs2) = m_robs
	m_obs = c(m_obs,m_obs2)
	m_obs = m_obs[order(as.numeric(names(m_obs)))]

	p_robs = rr[!(rr %in% names(p_obs))]
	p_obs2 = rep(0, length(p_robs))
	names(p_obs2) = p_robs
	p_obs = c(p_obs,p_obs2)
	p_obs = p_obs[order(as.numeric(names(p_obs)))]

	fracnonzero_m = 1 - (m_obs[1] / sum(m_obs))
	fracnonzero_p = 1 - (p_obs[1] / sum(p_obs))

	# fit both distributions to ZINB and NB
	zinb_fit_m = fitdist(mcounts, 'mother', disttouse = 'ZINB'); logL_zinb_m = zinb_fit_m[[5]]; aic_zinb_m = 6 - 2*logL_zinb_m
	nb_fit_m = fitdist(mcounts, 'mother', disttouse = 'NB'); logL_nb_m = nb_fit_m[[5]]; aic_nb_m = 4 - 2*logL_nb_m
	zinb_fit_p = fitdist(pcounts, 'father', disttouse = 'ZINB'); logL_zinb_p = zinb_fit_p[[5]]; aic_zinb_p = 6 - 2*logL_zinb_p
	nb_fit_p = fitdist(pcounts, 'father', disttouse = 'NB'); logL_nb_p = nb_fit_p[[5]]; aic_nb_p = 4 - 2*logL_nb_p
	
	# get fitted distributions to H1 for both ZINB and NB
	m_fit_zinb = getProbZinbConv(as.numeric(names(m_obs)), zinb_fit_m[[1]], zinb_fit_m[[2]], zinb_fit_m[[3]])*sum(m_obs)
	m_fit_nb = dNBI(as.numeric(names(m_obs)), mu = 2*nb_fit_m[[2]], sigma = 0.5*nb_fit_m[[3]])*sum(m_obs)
	p_fit_zinb = dZINBI(as.numeric(names(p_obs)),mu = zinb_fit_p[[2]], nu = zinb_fit_p[[1]], sigma = zinb_fit_p[[3]])*sum(p_obs)
	p_fit_nb = dNBI(as.numeric(names(p_obs)), mu = nb_fit_p[[2]], sigma = nb_fit_p[[3]])*sum(p_obs)
	
	# ZINB fit dfs
	mydf1 <- data.frame(vals = c(fracnonzero_m,fracnonzero_p), names=c("Mat","Pat"), expected = c(1-(m_fit_zinb[1]/sum(m_fit_zinb)),1-(p_fit_zinb[1]/sum(p_fit_zinb))))
	m_obs_nz = m_obs[names(m_obs) != '0']; p_obs_nz = p_obs[names(p_obs) != '0']
	mydata1 = data.frame(group = c(rep("Mat",maxv),rep("Pat",maxv)), x = rep(as.numeric(names(m_obs_nz)), 2), y = c(m_obs_nz,p_obs_nz*-1), fit = c(m_fit_zinb[2:length(m_fit_zinb)],-1*p_fit_zinb[2:length(p_fit_zinb)]))

	# NB fit dfs
	mydf2 <- data.frame(vals = c(fracnonzero_m,fracnonzero_p), names=c("Mat","Pat"), expected = c(1-(m_fit_nb[1]/sum(m_fit_nb)),1-(p_fit_nb[1]/sum(p_fit_nb))))
	mydata2 = data.frame(group = c(rep("Mat",maxv),rep("Pat",maxv)), x = rep(as.numeric(names(m_obs_nz)), 2), y = c(m_obs_nz,p_obs_nz*-1), fit = c(m_fit_nb[2:length(m_fit_nb)],-1*p_fit_nb[2:length(p_fit_nb)]))
	
	# Get limits for y-axis for consistency between plots
	maxyLHS = max(mydf1$vals, mydf1$expected, mydf2$vals, mydf2$expected)
	maxyRHS = max(mydata1$fit, mydata2$fit)
	minyRHS = -1*maxyRHS
	
	# build plot
	if (opt$skipplot == FALSE) {
		pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 10)
		pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 12)
		cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
		pdf(paste(outprefix,"_plot.pdf",sep=''), width = pwidth, height = pheight)

		# ZINB plot
		p1 = ggplot(data=mydf1, aes(x=names, y=vals, fill=names)) + geom_bar(stat="identity") + ylab("Fraction of nonzero counts") + xlab("") + geom_point(aes(x=names, y=expected))+ theme_bw() + ylim(0,maxyLHS)
		p2 = ggplot(mydata1, aes(x=x, y=y, fill=group)) + geom_bar(stat="identity", position="identity") + xlab("Nonzero read counts") + ylab("Cells with indicated read count") + geom_line(aes(x=x, y=fit))+ theme_bw() + ylim(minyRHS,maxyRHS)

		tg1 = textGrob(paste("Mat: (1-nu)=",round(1-zinb_fit_m[[1]],4),'\n',"Pat: (1-nu)=",round(1-zinb_fit_p[[1]],4),sep=''), gp=gpar(fontsize=14), just=c("center", "center"))
		tg2 = textGrob(paste("Mat: mu=",round(zinb_fit_m[[2]],4),'; sigma=',round(zinb_fit_m[[3]],4),'\n',"Pat: mu=",round(zinb_fit_p[[2]],4),'; sigma=',round(zinb_fit_p[[3]],4),sep=''), gp=gpar(fontsize=14), just=c("center", "center"))
		tg3 = textGrob(paste("Mat, H1: logL = ",round(logL_zinb_m,4),", AIC:",round(aic_zinb_m,4),'\n',"Pat, H1: logL = ",round(logL_zinb_p,4),", AIC:",round(aic_zinb_p,4),sep=''), gp=gpar(fontsize=14), just=c("center", "center"))

		# NB plot
		p3 = ggplot(data=mydf2, aes(x=names, y=vals, fill=names)) + geom_bar(stat="identity") + ylab("Fraction of nonzero counts") + xlab("") + geom_point(aes(x=names, y=expected))+ theme_bw() + ylim(0,maxyLHS)
		p4 = ggplot(mydata2, aes(x=x, y=y, fill=group)) + geom_bar(stat="identity", position="identity") + xlab("Nonzero read counts") + ylab("Cells with indicated read count") + geom_line(aes(x=x, y=fit))+ theme_bw() + ylim(minyRHS,maxyRHS)

		tg4 = textGrob(paste("Mat: mu=",round(nb_fit_m[[2]],4),'; sigma=',round(nb_fit_m[[3]],4),'\n',"Pat: mu=",round(nb_fit_p[[2]],4),'; sigma=',round(nb_fit_p[[3]],4),sep=''), gp=gpar(fontsize=14), just=c("center", "center"))
		tg5 = textGrob(paste("Mat, H1: logL = ",round(logL_nb_m,4),", AIC:",round(aic_nb_m,4),'\n',"Pat, H1: logL = ",round(logL_nb_p,4),", AIC:",round(aic_nb_p,4),sep=''), gp=gpar(fontsize=14), just=c("center", "center"))

		if (! is.null(xmax)) {
			p2 = p2 + xlim(0,xmax)
			p4 = p4 + xlim(0,xmax)
			p6 = p6 + xlim(0,xmax)
		}

		lay = as.matrix(rbind(c(1,2,3),c(4,5,5),c(NA,6,7),c(8,9,9)))
		grid.arrange(tg1,tg2,tg3,p1,p2,tg4,tg5,p3,p4,layout_matrix = lay, heights=c(0.1, 0.85, 0.1, 0.85))
		graphics.off()
	}

	if (opt$outputdatalong == TRUE || opt$outputdatawide == TRUE) {
		write.table(mydf1,file=paste(outprefix,"_ZINB_zero.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
		write.table(mydata1,file=paste(outprefix,"_ZINB_nonzero.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
		write.table(mydf2,file=paste(outprefix,"_NB_zero.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
		write.table(mydata2,file=paste(outprefix,"_NB_nonzero.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	}
}

# Bar chart showing allelic counts over clusters, for a single gene
# Inputs:
# 	- mcounts : integer vector of values corresponding to maternal allelic counts, with names() corresponding to nuclei
# 	- pcounts : integer vector of values corresponding to paternal allelic counts, with names() corresponding to nuclei
# 	- outprefix : string prefix for output files
# 	- clusters : optional dataframe with 1 column corresponding to cluster assignments, where each row corresponds to a nucleus and rowname is the nucleus ID
make_bar_plot_impr <- function(mcounts, pcounts, outprefix, clusters=NULL, xorder = NULL) {

	# get order of x-axis if provided
	if (! is.null(xorder)) {
		xorder = unlist(strsplit(xorder,",")[[1]])
		if (length(xorder) != length(unique(xorder))) stop("--xorder list cannot contain repeated values")
		if (! is.null(clusters)) {
			if(isFALSE(all(xorder %in% clusters$cluster))) stop("--xorder contains value(s) not found in second column of --sampfile\n")
		} else {
			cat("Warning: --xorder values provided, but will be ignored for 'bar' plot unless --sampfile also provided\n")
		}
	}
	
	# convert to vectors
	mcounts_sort = as.numeric(mcounts); names(mcounts_sort) = colnames(mcounts)
	pcounts_sort = as.numeric(pcounts); names(pcounts_sort) = colnames(pcounts)
		
	mydata = cbind(mcounts_sort,pcounts_sort)

	if (! is.null(clusters)) {
		mydata = merge(mydata,clusters,by=0)
		
		if (! is.null(xorder)) mydata = mydata[order(match(mydata$cluster,xorder)),]
		if (is.null(xorder)) mydata = mydata[order(mydata$cluster),]
		
		xlabels = unique(mydata$cluster)
		xlabpos = unlist(lapply(xlabels, function(x) mean(which(mydata$cluster == x))))
				
		mydata$ordertouse = 1:nrow(mydata)	
		mdata = mydata[,colnames(mydata) %in% c("Row.names","mcounts_sort","cluster","ordertouse")]
		pdata = mydata[,colnames(mydata) %in% c("Row.names","pcounts_sort","cluster","ordertouse")]

		pdata$pcounts_sort = -1*pdata$pcounts_sort
		colnames(mdata)[2] = "counts"
		colnames(pdata)[2] = "counts"
	} else {
		mydata = as.data.frame(mydata)

		mydata$ordertouse = 1:nrow(mydata)
		mdata = mydata[,colnames(mydata) %in% c("mcounts_sort","ordertouse")]
		pdata = mydata[,colnames(mydata) %in% c("pcounts_sort","ordertouse")]

		pdata$pcounts_sort = -1*pdata$pcounts_sort
		colnames(mdata)[1] = "counts"
		colnames(pdata)[1] = "counts"
	}
	
	mdata$group = rep("Mat",nrow(mdata))
	pdata$group = rep("Pat",nrow(pdata))
	
	mydata = rbind(mdata,pdata)
	
	maxy = max(abs(c(mydata$counts)))
	miny = -1*maxy
	
	if (opt$skipplot == FALSE) {
		pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 8)
		pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 18)
		cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
		pdf(paste(outprefix,'_plot.pdf',sep=''), width = pwidth, height = pheight)
		p1 = ggplot(mydata, aes(x=ordertouse, y=counts, fill=group)) + geom_bar(stat="identity", position="identity") + xlab("Individual Nuclei") + ylab("Read counts") + theme_bw() + ylim(miny,maxy)
	
		if (! is.null(clusters)) {
			pos = c(-0.5)
			for (i in unique(mydata$cluster)[1:length(unique(mydata$cluster))-1]) {
				pos = c(pos,which(mdata$cluster == i)[length(which(mdata$cluster == i))] + 0.5)
			}
			pos = c(pos,nrow(mdata)+0.5)
			cluslines <- data.frame(clus = c(0,unique(mydata$cluster)), pos = pos)
			p1 = p1 + geom_vline(aes(xintercept=pos), cluslines, linetype = "dashed")
			p1 = p1 + scale_x_continuous(breaks = xlabpos, labels = xlabels)
		}

		show(p1)
		graphics.off()
	}

	if (opt$outputdatalong == TRUE || opt$outputdatawide == TRUE) {
		write.table(mydata,file=paste(outprefix,"_counts.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	}
}

# Bar chart showing total expression over clusters, for a single gene
# Inputs:
# 	- acounts : integer vector of values corresponding to total counts, or CPM values, with names() corresponding to nuclei
# 	- outprefix : string prefix for output files
# 	- clusters : optional dataframe with 1 column corresponding to cluster assignments, where each row corresponds to a nucleus and rowname is the nucleus ID
make_bar_plot_expr <- function(acounts, outprefix, clusters=NULL, xorder = NULL) {

	# get order of x-axis if provided
	if (! is.null(xorder)) {
		xorder = unlist(strsplit(xorder,",")[[1]])
		if (length(xorder) != length(unique(xorder))) stop("--xorder list cannot contain repeated values")
		if (! is.null(clusters)) {
			if(isFALSE(all(xorder %in% clusters$cluster))) stop("--xorder contains value(s) not found in second column of --sampfile\n")
		} else {
			cat("Warning: --xorder values provided, but will be ignored for 'bar' plot unless --sampfile also provided\n")
		}
	}
	
	# convert to vectors
	mydata = as.numeric(acounts); names(mydata) = colnames(acounts)
	
	if (! is.null(clusters)) {
		mydata = merge(mydata,clusters,by=0)
				
		if (! is.null(xorder)) mydata = mydata[order(match(mydata$cluster,xorder)),]
		if (is.null(xorder)) mydata = mydata[order(mydata$cluster),]
		
		xlabels = unique(mydata$cluster)
		xlabpos = unlist(lapply(xlabels, function(x) mean(which(mydata$cluster == x))))
		
		mydata$ordertouse = 1:nrow(mydata)	
		colnames(mydata)[2] = "counts"
	} else {
		mydata = as.data.frame(mydata)

		mydata$ordertouse = 1:nrow(mydata)
		colnames(mydata)[1] = "counts"
	}
	
	maxy = max(abs(c(mydata$counts)))
	
	if (opt$skipplot == FALSE) {
		pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 8)
		pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 18)
		cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
		pdf(paste(outprefix,'_plot.pdf',sep=''), width = pwidth, height = pheight)
		p1 = ggplot(mydata, aes(x=ordertouse, y=counts)) + geom_bar(stat="identity", position="identity") + xlab("Individual Nuclei") + ylab("Read counts") + theme_bw() + ylim(0,maxy)
	
		if (! is.null(clusters)) {
			pos = c(-0.5)
			for (i in unique(mydata$cluster)[1:length(unique(mydata$cluster))-1]) {
				pos = c(pos,which(mydata$cluster == i)[length(which(mydata$cluster == i))] + 0.5)
			}
			pos = c(pos,nrow(mydata)+0.5)
			cluslines <- data.frame(clus = c(0,unique(mydata$cluster)), pos = pos)
			p1 = p1 + geom_vline(aes(xintercept=pos), cluslines, linetype = "dashed")
			p1 = p1 + scale_x_continuous(breaks = xlabpos, labels = xlabels)
		}

		show(p1)
		graphics.off()
	}

	if (opt$outputdatalong == TRUE || opt$outputdatawide == TRUE) {
		write.table(mydata,file=paste(outprefix,"_counts.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	}
}

# Dot plot of individual nuclei showing allelic bias over all informative nuclei for gene(s); can accept multiple genes as input
make_nuc_plot <- function(mcounts, pcounts, outprefix, clusters=NULL, minreads = 1) {

	totreads = mcounts + pcounts
	fracmat = mcounts / totreads

	numgenes = nrow(totreads)
	genenames = rownames(totreads)
	
	if (ngenes > 1) {
		mm = max(rowSums(totreads >= minreads))
	} else {
		mm = sum(totreads >= minreads)
	}
	
	if (is.null(clusters)) {
		
		maxnonzero = max(rowSums(totreads >= minreads))
		if (maxnonzero == 0) {
			stop("No data for any of these --genes in provided count matrices; aborting.")
		}
		y = 1
	
		for (gene in rownames(totreads)) {
			totreads_tmp = totreads[rownames(totreads) == gene,]
			fracmat_tmp = fracmat[rownames(fracmat) == gene,]
			i = unname(which(totreads_tmp >= minreads))
			if (length(i) == 0) {
				# no nuclei with this minimum coverage, plot will show no dots for this gene
				mydata_tmp = data.frame(tot = rep(NA,maxnonzero), frac = rep(NA,maxnonzero))
			} else {
				totreads_tmp = unname(totreads_tmp[i])
				fracmat_tmp = unname(fracmat_tmp[i])
		
				# sort so that the most paternal nuclei are left, maternal right
				totreads_tmp = totreads_tmp[order(fracmat_tmp)]
				fracmat_tmp = fracmat_tmp[order(fracmat_tmp)]
				mydata_tmp = data.frame(tot = totreads_tmp, frac = fracmat_tmp)
				mydata_tmp = mydata_tmp[order(mydata_tmp$frac),]
				
				# pad vector so that length is same as longest vector
				ll = length(totreads_tmp)
		
				mydata_pad = data.frame(tot = rep(NA,maxnonzero - ll), frac = rep(NA,maxnonzero - ll))
				mydata_tmp = rbind(mydata_tmp,mydata_pad)		
			}
			
			mydata_tmp$x = 1:maxnonzero
			mydata_tmp$y = rep(y,maxnonzero)
		
			# add to main data frame
			if (exists("mydata")) {
				mydata = rbind(mydata, mydata_tmp)
			} else {
				mydata = mydata_tmp
			}		
			y = y + 1 			
		}

	} else {
		
		totreads = merge(clusters,t(totreads),by=0)
		rownames(totreads) = totreads$Row.names
		totreads$Row.names = NULL
		
		fracmat = merge(clusters,t(fracmat),by=0)
		rownames(fracmat) = fracmat$Row.names
		fracmat$Row.names = NULL
		
		cluslist = sort(unique(totreads$cluster))
		xoffset = 0; splitters = c(); xlblvals = c(); xlbls = c()
		
		for (clus in cluslist) {
			tt = totreads[totreads$cluster == clus,]; tt$cluster = NULL; tt = t(tt)
			ff = fracmat[fracmat$cluster == clus,]; ff$cluster = NULL; ff = t(ff)

			maxnonzero = max(rowSums(tt >= minreads))
			if (maxnonzero == 0) maxnonzero = 1	
			y = 1; xlblvals = c(xlblvals, xoffset + (maxnonzero+1) / 2); xlbls = c(xlbls,as.character(clus))
			
			for (gene in rownames(tt)) {	
				totreads_tmp = tt[rownames(tt) == gene,]
				fracmat_tmp = ff[rownames(ff) == gene,]
				i = unname(which(totreads_tmp >= minreads))
				
				if (length(i) == 0) {
					# no nuclei with this minimum coverage in this cluster
					mydata_tmp = data.frame(tot = rep(NA,maxnonzero), frac = rep(NA,maxnonzero))		
				} else {
					totreads_tmp = unname(totreads_tmp[i])
					fracmat_tmp = unname(fracmat_tmp[i])
						
					# sort so that the most paternal nuclei are left, maternal right
					totreads_tmp = totreads_tmp[order(fracmat_tmp)]
					fracmat_tmp = fracmat_tmp[order(fracmat_tmp)]
					mydata_tmp = data.frame(tot = totreads_tmp, frac = fracmat_tmp)
					mydata_tmp = mydata_tmp[order(mydata_tmp$frac),]
				
					# pad vector so that length is same as longest vector
					ll = length(totreads_tmp)
		
					mydata_pad = data.frame(tot = rep(NA,maxnonzero - ll), frac = rep(NA,maxnonzero - ll))
					mydata_tmp = rbind(mydata_tmp,mydata_pad)		
				}
				
				mydata_tmp$x = xoffset + 1:maxnonzero
				mydata_tmp$y = rep(y,maxnonzero)
						
				# add to main data frame
				if (exists("mydata")) {
					mydata = rbind(mydata, mydata_tmp)
				} else {
					mydata = mydata_tmp
				}		
				
				y = y + 1 			
			}
			
			xoffset = xoffset + maxnonzero + 1
			splitters = c(splitters, xoffset)			
		}		
	}
	
	if (opt$skipplot == FALSE) {
		pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 3 + numgenes/4)
		pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 5 + (mm / 20))
		cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
		pdf(paste(outprefix,'_plot.pdf',sep=''), width = pwidth, height = pheight)
		p1 = ggplot(mydata, aes(x=x, y=y, color = frac)) + geom_point(aes(size=tot)) + ylab("Genes") + theme_bw()
		p1 = p1 + scale_size_continuous(range = c(0,4)) + scale_colour_gradientn(colours=c("blue","yellow","red"), values = c(0,0.667,1))
		p1 = p1 + scale_y_continuous(breaks=1:numgenes, labels=genenames, limits = c(0.5,numgenes+0.5)) + theme(panel.grid.minor.y = element_blank())
		if (! is.null(clusters)) {
			p1 = p1 + geom_vline(xintercept = splitters) + scale_x_continuous(name = "Clusters of individual nuclei", breaks = xlblvals, labels = xlbls)
		} else {
			p1 = p1 + xlab("Individual Nuclei")
		}
		show(p1)
		graphics.off()
	}

	if (opt$outputdatalong == TRUE || opt$outputdatawide == TRUE) {
		write.table(mydata,file=paste(outprefix,"_data.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	}
}

# Dot plot of average expression over individual nuclei or groups of nuclei, for individual or groups of genes
# If groups of nuclei or genes: size is the fraction of nuclei/genes expressed from that group. 
# Plots of single genes across single samples: all dots same size, dot color indicates expression
# Plots of single genes across sample clusters: dot size is proportional to fraction of samples in cluster with nonzero expression, color is average expression across all nonzero samples
make_dot_plot_expr <- function(expr_data, outprefix, clusters=NULL, genegroups=NULL, maxdims = 100, xorder = NULL, yorder = NULL, expr_type = "", noprelog = FALSE, includezeros = FALSE, fmax = NULL, smax = NULL) {

	# stop if xorder or yorder have repeats of same gene(s) or samples
	if (!is.null(xorder)) xorder = strsplit(xorder,",")[[1]]
	if (!is.null(yorder)) yorder = rev(strsplit(yorder,",")[[1]])
	if (! is.null(xorder) && length(xorder) != length(unique(xorder))) stop("--xorder list cannot contain repeated values")
	if (! is.null(yorder) && length(yorder) != length(unique(yorder))) stop("--yorder list cannot contain repeated values")

	# stop if more than maxdims genes/gene groups or maxdims samples/sample clusters requested, plot will be unintelligible
	if (is.null(genegroups)) {
		if (nrow(expr_data) > maxdims) stop("Error: dot plot would have ",nrow(expr_data)," rows (genes), which is greater than the value of --maxdims ",maxdims,". Please either reduce the number of genes or increase --maxdims.")
		if (! is.null(yorder) && setequal(rownames(expr_data),yorder) == FALSE) stop("if using --yorder to order genes, same genes must be in --yorder and in --genes/--genefile or count matrices")
	} else {
		ggroups = unique(genegroups[,1])
		if (length(ggroups) > maxdims) stop("Error: dot plot would have ",length(ggroups)," rows (gene groups), which is greater than the value of --maxdims ",maxdims,". Please either reduce the number of genes or increase --maxdims.")
		if (! is.null(yorder) && setequal(ggroups,yorder) == FALSE) stop("if using --yorder to order gene groups, same gene groups must be in --yorder and in --genes/--genefile or count matrices")
	}
	if (is.null(clusters)) {
		if (ncol(expr_data) > maxdims) stop("Error: dot plot would have ",ncol(expr_data)," columns (samples), which is greater than the value of --maxdims ",maxdims,". Please either reduce the number of samples or increase --maxdims.")
		if (! is.null(xorder) && setequal(colnames(expr_data),xorder) == FALSE) stop("if using --xorder to order samples, same samples must be in --sampfile or count matrices")
	} else {
		cclusters = unique(clusters[,1])
		if (length(cclusters) > maxdims) stop("Error: dot plot would have ",length(cclusters)," columns (sample clusters), which is greater than the value of --maxdims ",maxdims,". Please either reduce the number of sample clusters or increase --maxdims.")
		if (! is.null(xorder) && setequal(cclusters,xorder) == FALSE) stop("if using --xorder to order sample clusters, same sample clusters must be in --xorder and in --sampfile or count matrices")
	}
	
	# take log2 if requested
	if (noprelog == FALSE) expr_data = log2(expr_data + 1)		# add pseudocount of 1; since the pseudocount is always 1, log2(1) = 0 so code below for including/excluding zeros, etc. still works
	if (noprelog == FALSE) cat("Taking log2 of data before averaging\n")
				
	# get expression averages and frac_nonzero across all sample clusters (if provided)
	if (is.null(clusters)) {
		mydata_mean = expr_data
		mydata_frac = ifelse(expr_data > 0,1,0)
	} else {
		mydata = merge(clusters,t(expr_data),by=0)
		rownames(mydata) = mydata$Row.names
		mydata$Row.names = NULL
		
		mydata_frac = aggregate(mydata[, 2:ncol(mydata)] > 0, list(mydata[,1]), mean)
		rownames(mydata_frac) = mydata_frac$Group.1
		mydata_frac$Group.1 = NULL
		
		if (includezeros == FALSE) {
			mydata_mean = aggregate(mydata[, 2:ncol(mydata)], list(mydata[,1]), FUN=(function(x){ifelse(sum(x==0)>0 & sum(x !=0) >0, mean(x[x>0]), mean(x))}))
		} else {
			mydata_mean = aggregate(mydata[, 2:ncol(mydata)], list(mydata[,1]), FUN=mean)
		}
		rownames(mydata_mean) = mydata_mean$Group.1
		mydata_mean$Group.1 = NULL
	
		mydata_mean = t(mydata_mean)
		mydata_frac = t(mydata_frac)
	}
	
	# get expression averages and frac_nonzero across all gene groups (if provided)
	if (! is.null(genegroups)) {	
		mydata_frac = merge(genegroups,mydata_frac,by=0)
		rownames(mydata_frac) = mydata_frac$Row.names
		mydata_frac$Row.names = NULL

		mydata_mean = merge(genegroups,mydata_mean,by=0)
		rownames(mydata_mean) = mydata_mean$Row.names
		mydata_mean$Row.names = NULL
		
		if (! is.null(clusters)) {
			mydata_frac = aggregate(mydata_frac[, 2:ncol(mydata_frac)], list(mydata_frac[,1]), mean)
			mydata_mean = aggregate(mydata_mean[, 2:ncol(mydata_mean)], list(mydata_mean[,1]), mean)
		} else {
			mydata_frac = aggregate(mydata_frac[, 2:ncol(mydata_frac)] > 0, list(mydata_frac[,1]), mean)
			if (includezeros == FALSE) {
				mydata_mean = aggregate(mydata_mean[, 2:ncol(mydata_mean)], list(mydata_mean[,1]), FUN=(function(x){ifelse(sum(x==0)>0 & sum(x !=0) >0, mean(x[x>0]), mean(x))}))
			} else {
				mydata_mean = aggregate(mydata_mean[, 2:ncol(mydata_mean)], list(mydata_mean[,1]), FUN=mean)
			}
		}
		rownames(mydata_frac) = mydata_frac$Group.1
		mydata_frac$Group.1 = NULL
		rownames(mydata_mean) = mydata_mean$Group.1
		mydata_mean$Group.1 = NULL
	}
			
	# reshape matrices for plot to make final dataset
	mydata_final = data.frame(samp = c(), gene = c(), mmean = c(), frac = c())
	for (samp in colnames(mydata_mean)) {
		mydata_tmp = data.frame(samp = rep(samp,nrow(mydata_mean)), gene = rownames(mydata_mean), mmean = unname(mydata_mean[,samp]), frac = unname(mydata_frac[,samp]))
		mydata_final = rbind(mydata_final,mydata_tmp)
	}
	
	if (is.null(xorder)) xorder = unique(mydata_final$samp)
	if (is.null(yorder)) yorder = unique(mydata_final$gene)
	
	xorder = data.frame(samp = xorder, x = 1:length(xorder))
	yorder = data.frame(gene = yorder, y = 1:length(yorder))
		
	mydata_final = merge(mydata_final,xorder,by="samp")
	mydata_final = merge(mydata_final,yorder,by="gene")
	
	if (opt$samesize == TRUE) {
		mydata_final$frac = 1
	}
	
	# make plot
	if (opt$skipplot == FALSE) {
		pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 2 + (nrow(yorder) / 4))
		pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 4 + (nrow(xorder) / 2))
		if (is.null(fmax)) {
			fmax = max(mydata_final$mmean)
		} else {
			mydata_final$mmean[mydata_final$mmean > fmax] = fmax
		}
		if (is.null(smax)) {
			smax = max(mydata_final$frac)
		} else {
			mydata_final$frac[mydata_final$frac > smax] = smax
		}
		cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
		pdf(paste(outprefix,'_plot.pdf',sep=''), width = pwidth, height = pheight)
		p1 = ggplot(mydata_final, aes(x=x, y=y)) + geom_point(aes(size=frac, fill=mmean),pch=21) + theme_bw()
		p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax)) + scale_fill_gradientn(colours=c("royalblue", "khaki", "firebrick"), limits=c(0,fmax))
		p1 = p1 + scale_y_continuous(breaks=1:nrow(yorder), labels=yorder$gene, limits = c(0.5,nrow(yorder)+0.5)) + theme(panel.grid.minor.y = element_blank())
		p1 = p1 + scale_x_continuous(breaks=1:nrow(xorder), labels=xorder$samp, limits = c(0.5,nrow(xorder)+0.5)) + theme(panel.grid.minor.x = element_blank())
		if (is.null(clusters)) {
			p1 = p1 + xlab("Samples")
		} else {
			p1 = p1 + xlab("Sample clusters")
		}	
		if (is.null(genegroups)) {
			p1 = p1 + ylab("Genes")
		} else {
			p1 = p1 + ylab("Gene groups")
		}	
	
		if (expr_type == "CPM") {
			if (noprelog == FALSE) p1 = p1 + labs(fill = "avg. log2(CPM)", size = "frac. non-zero")
			if (noprelog == TRUE) p1 = p1 + labs(fill = "avg. CPM", size = "frac. non-zero")
		} else if (expr_type == "counts") {
			if (noprelog == FALSE) p1 = p1 + labs(fill = "avg. log2(norm. read counts)", size = "frac. non-zero")
			if (noprelog == TRUE) p1 = p1 + labs(fill = "avg. norm. read counts", size = "frac. non-zero")
		} else {
			if (noprelog == FALSE) p1 = p1 + labs(fill = "avg. log2(expression)", size = "frac. non-zero")
			if (noprelog == TRUE) p1 = p1 + labs(fill = "avg. expression", size = "frac. non-zero")
		} 

		show(p1)
		graphics.off()
	}
	
	if (opt$outputdatalong == TRUE) {
		write.table(mydata_final,file=paste(outprefix,"_data.txt",sep=''),col.names = TRUE,row.names = FALSE,quote = F,sep='\t')
	}
	if (opt$outputdatawide == TRUE) {
		mydata_final_mean = mydata_final[,1:3]
		mydata_final_mean = dcast(mydata_final_mean, gene ~ samp, value.var = "mmean")
		mydata_final_frac = mydata_final[,c(1,2,4)]
		mydata_final_frac = dcast(mydata_final_frac, gene ~ samp, value.var = "frac")
		
		write.table(mydata_final_mean,file=paste(outprefix,"_mean_data.txt",sep=''),col.names = TRUE,row.names = FALSE,quote = F,sep='\t')
		write.table(mydata_final_frac,file=paste(outprefix,"_frac_data.txt",sep=''),col.names = TRUE,row.names = FALSE,quote = F,sep='\t')
	}
}


# Dot plot of % maternal over individual nuclei or groups of nuclei, for individual or groups of genes
# If groups of nuclei or genes: size is the fraction of nuclei/genes expressed from that group. 
# Plots of single genes across single samples: all dots same size, dot color indicates expression
# Plots of single genes across sample clusters: dot size is proportional to fraction of samples in cluster with nonzero expression, color is average expression across all nonzero samples
make_dot_plot_impr <- function(mcounts, pcounts, outprefix, clusters=NULL, genegroups=NULL, minreads = 1, maxdims = 100, xorder = NULL, yorder = NULL, smax = NULL) {

	# stop if xorder or yorder have repeats of same gene(s) or samples
	if (!is.null(xorder)) xorder = strsplit(xorder,",")[[1]]
	if (!is.null(yorder)) yorder = rev(strsplit(yorder,",")[[1]])
	if (! is.null(xorder) && length(xorder) != length(unique(xorder))) stop("--xorder list cannot contain repeated values")
	if (! is.null(yorder) && length(yorder) != length(unique(yorder))) stop("--yorder list cannot contain repeated values")

	# stop if more than maxdims genes/gene groups or maxdims samples/sample clusters requested, plot will be unintelligible
	if (is.null(genegroups)) {
		if (nrow(mcounts) > maxdims) stop("Error: dot plot would have ",nrow(mcounts)," rows (genes), which is greater than the value of --maxdims ",maxdims,". Please either reduce the number of genes or increase --maxdims.")
		if (! is.null(yorder) && setequal(rownames(mcounts),yorder) == FALSE) stop("if using --yorder to order genes, same genes must be in --yorder and in --genes/--genefile or count matrices")
	} else {
		ggroups = unique(genegroups[,1])
		print(ggroups)
		print(yorder)
		if (length(ggroups) > maxdims) stop("Error: dot plot would have ",length(ggroups)," rows (gene groups), which is greater than the value of --maxdims ",maxdims,". Please either reduce the number of genes or increase --maxdims.")
		if (! is.null(yorder) && setequal(ggroups,yorder) == FALSE) stop("if using --yorder to order gene groups, same gene groups must be in --yorder and in --genes/--genefile or count matrices")
	}
	if (is.null(clusters)) {
		if (ncol(mcounts) > maxdims) stop("Error: dot plot would have ",ncol(mcounts)," columns (samples), which is greater than the value of --maxdims ",maxdims,". Please either reduce the number of samples or increase --maxdims.")
		if (! is.null(xorder) && setequal(colnames(mcounts),xorder) == FALSE) stop("if using --xorder to order samples, same samples must be in --sampfile or count matrices")
	} else {
		cclusters = unique(clusters[,1])
		if (length(cclusters) > maxdims) stop("Error: dot plot would have ",length(cclusters)," columns (sample clusters), which is greater than the value of --maxdims ",maxdims,". Please either reduce the number of sample clusters or increase --maxdims.")
		if (! is.null(xorder) && setequal(cclusters,xorder) == FALSE) stop("if using --xorder to order sample clusters, same sample clusters must be in --xorder and in --sampfile or count matrices")
	}
	
	# get % maternal
	pmat = mcounts / (mcounts + pcounts)
	tot = mcounts + pcounts
	
	# filter by minreads (e.g. drop all gene/cell combos with fewer than minreads total allelic reads)
	pmat = ifelse(tot >= minreads, pmat, NaN)
	tot = ifelse(tot >= minreads, tot, 0)
			
	# get average %mat and fraction nonzero over clusters/groups
	if (is.null(clusters)) {
		pmat_mean = pmat
		fracnonzero = ifelse(tot > 0,1,0)
	} else {
		ppmat = merge(clusters,t(pmat),by=0)
		pptot = merge(clusters,t(tot),by=0)
		rownames(ppmat) = ppmat$Row.names
		rownames(pptot) = pptot$Row.names
		ppmat$Row.names = NULL
		pptot$Row.names = NULL
			
		pmat_mean = aggregate(ppmat[, 2:ncol(ppmat)], list(ppmat[,1]), mean, na.rm=TRUE)
		rownames(pmat_mean) = pmat_mean$Group.1
		pmat_mean$Group.1 = NULL
				
		ptot_mask = pptot
		ptot_mask[2:ncol(pptot)] = ifelse(ptot_mask[2:ncol(ptot_mask)] > 0,1,0)
		fracnonzero = aggregate(ptot_mask[, 2:ncol(ptot_mask)], list(ptot_mask[,1]), mean)
		rownames(fracnonzero) = fracnonzero$Group.1
		fracnonzero$Group.1 = NULL
		
		pmat_mean = t(pmat_mean)
		fracnonzero = t(fracnonzero)
	}	
	
	# get expression averages and frac_nonzero across all gene groups (if provided)
	if (! is.null(genegroups)) {	
		pmat_mean = merge(genegroups,pmat_mean,by=0)
		print(head(pmat_mean))
		rownames(pmat_mean) = pmat_mean$Row.names
		pmat_mean$Row.names = NULL
		
		pmat_mean = aggregate(pmat_mean[, 2:ncol(pmat_mean)] > 0, list(pmat_mean[,1]), mean, na.rm=TRUE)
		rownames(pmat_mean) = pmat_mean$Group.1
		pmat_mean$Group.1 = NULL

		fracnonzero = merge(genegroups,fracnonzero,by=0)
		rownames(fracnonzero) = fracnonzero$Row.names
		fracnonzero$Row.names = NULL

		fracnonzero = aggregate(fracnonzero[, 2:ncol(fracnonzero)], list(fracnonzero[,1]), mean)
		rownames(fracnonzero) = fracnonzero$Group.1
		fracnonzero$Group.1 = NULL
	}
			
	# reshape matrices for plot to make final dataset
	mydata_final = data.frame(samp = c(), gene = c(), mmean = c(), frac = c())
	for (samp in colnames(pmat_mean)) {
		mydata_tmp = data.frame(samp = rep(samp,nrow(pmat_mean)), gene = rownames(pmat_mean), mmean = unname(pmat_mean[,samp]), frac = unname(fracnonzero[,samp]))
		mydata_final = rbind(mydata_final,mydata_tmp)
	}
		
	if (is.null(xorder)) xorder = unique(mydata_final$samp)
	if (is.null(yorder)) yorder = unique(mydata_final$gene)

	xorder = data.frame(samp = xorder, x = 1:length(xorder))
	yorder = data.frame(gene = yorder, y = 1:length(yorder))
		
	mydata_final = merge(mydata_final,xorder,by="samp")
	mydata_final = merge(mydata_final,yorder,by="gene")
	
	if (opt$samesize == TRUE) {
		mydata_final$frac = 1
	}

	# make plot
	if (opt$skipplot == FALSE) {
		pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 2 + (nrow(yorder) / 4))
		pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 4 + (nrow(xorder) / 2))
		cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
		pdf(paste(outprefix,'_plot.pdf',sep=''), width = pwidth, height = pheight)
		p1 = ggplot(mydata_final, aes(x=x, y=y)) + geom_point(aes(size=frac, fill=mmean),pch=21) + theme_bw()
		p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize)) + scale_fill_gradientn(colours=c("royalblue", "khaki", "firebrick"), values = c(0,0.66667,1))
		p1 = p1 + scale_y_continuous(breaks=1:nrow(yorder), labels=yorder$gene, limits = c(0.5,nrow(yorder)+0.5)) + theme(panel.grid.minor.y = element_blank())
		p1 = p1 + scale_x_continuous(breaks=1:nrow(xorder), labels=xorder$samp, limits = c(0.5,nrow(xorder)+0.5)) + theme(panel.grid.minor.x = element_blank())
		if (is.null(clusters)) {
			p1 = p1 + xlab("Samples")
		} else {
			p1 = p1 + xlab("Sample clusters")
		}	
		if (is.null(genegroups)) {
			p1 = p1 + ylab("Genes")
		} else {
			p1 = p1 + ylab("Gene groups")
		}	
		p1 = p1 + labs(size = "frac. non-zero", fill = "avg. % maternal")
	
		show(p1)
		graphics.off()
	}

	if (opt$outputdatalong == TRUE) {
		write.table(mydata_final,file=paste(outprefix,"_data.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	}
}
 
# Line plot of average paternal and maternal expression, over samples (or clusters if provided)
# Can optionally add an extra step where paternal reads are fitted to a ZINB/NB distribution, then 
# simulations equivalent to two draws from this distribution are made to match ploidy with maternal counts, and make the
# mat and pat distributions more comparable (done by cluster, must have cluster)
make_lin_plot_impr <- function(mcounts, pcounts, outprefix, clusters, coef_nu = 1, coef_mu_ZINB = 1, coef_mu_NB = 1, maxdims = 100, xorder = NULL, noprelog = FALSE, postlog = FALSE, linewidth = 2, nreps = 500, addsims = FALSE, ymax = NULL, ymin = NULL, smax = NULL) {

	# stop if xorder or yorder have repeats of same gene(s) or samples
	if (!is.null(xorder)) xorder = strsplit(xorder,",")[[1]]
	if (! is.null(xorder) && length(xorder) != length(unique(xorder))) stop("--xorder list cannot contain repeated values")

	# stop if more than maxdims genes/gene groups or maxdims samples/sample clusters requested, plot will be unintelligible
	if (is.null(clusters)) stop("Plots of type 'lin' currently only work for clusters of samples (e.g. --sampfile must have 2 columns; to cheat, make column 2 identical to column 1")
	cclusters = sort(unique(clusters[,1]))
	if (length(cclusters) > maxdims) stop("Error: dot plot would have ",length(cclusters)," columns (sample clusters), which is greater than the value of --maxdims ",maxdims,"and would probably create a very crowded plot. Please either reduce the number of sample clusters or increase --maxdims.")
	if (! is.null(xorder) && setequal(cclusters,xorder) == FALSE) stop("if using --xorder to order sample clusters, must have same values as column 2 of --sampfile")
	if (! is.null(xorder)) cclusters = xorder
	
	pcounts = pcounts[,colnames(mcounts),drop=FALSE]

	# for each cluster, fit ZINB/NB and simulate drawing 2x from same distribution to adjust pat -> mat and
	# make comparable
	if (addsims == TRUE) {
		simpmeans = c()
		simpsds = c()
		for (i in 1:length(cclusters)) {
			samplist = rownames(clusters[clusters == cclusters[i],,drop=FALSE])
			psubset = unname(pcounts[,colnames(pcounts) %in% samplist])
			
			if (max(psubset) == 0) {
				simpmeans[i] = 0
				simpsds[i] = 0
			} else {
				# fit distribution
				pfit = fitdist(psubset,"father")
				nobs = length(psubset)			
			
				if (pfit[[4]] == "ZINB") {
					# estimate maternal params for ZINB based on the coefficients
					mat_nu = 1 - ((1-pfit[[1]])/coef_nu)
					mat_mu = pfit[[2]] / coef_mu_ZINB
					mat_sigma = pfit[[3]]

					# simulate nreps draws of nobs observations from the ZINB
					totS = 0; totSS = 0
					for (n in 1:nreps) {
						rand = rZINBI(nobs, mu = mat_mu, nu = mat_nu, sigma = mat_sigma) + rZINBI(nobs, mu = mat_mu, nu = mat_nu, sigma = mat_sigma)						
						rand = sum(rand) / nobs
						if (noprelog == FALSE) rand = log2(rand + 1)
						totS = totS + rand
						totSS = totSS + (rand^2)
					}
					simpmeans[i] = totS / nreps
					simpsds[i] = sqrt((totSS / nreps) - (simpmeans[i])^2)
				} else {
					# estimate maternal params for NB based on the coefficients
					mat_mu = pfit[[2]] / coef_mu_ZINB
					mat_sigma = pfit[[3]]

					# simulate nreps draws of nobs observations from the NB
					totS = 0; totSS = 0
					for (n in 1:nreps) {
						rand = rNBI(nobs, mu = mat_mu, sigma = mat_sigma) + rNBI(nobs, mu = mat_mu, sigma = mat_sigma)
						rand = sum(rand) / nobs
						if (noprelog == FALSE) rand = log2(rand + 1)
						totS = totS + rand
						totSS = totSS + (rand^2)
					}
					simpmeans[i] = totS / nreps
					simpsds[i] = sqrt((totSS / nreps) - (simpmeans[i])^2)
				}	
			}
		}
	}
	
	# get actual count distributions	
	if (noprelog == FALSE) {
		mcounts = log2(mcounts + 1)
		pcounts = log2(pcounts + 1)
	}

	clusters = clusters[colnames(mcounts),,drop=FALSE]
	mnonzero = aggregate(t(mcounts), list(clusters[,1]), FUN= function(x) sum(x>0)); rownames(mnonzero) = mnonzero[,1]; mnonzero[,1] = NULL
	pnonzero = aggregate(t(pcounts), list(clusters[,1]), FUN= function(x) sum(x>0)); rownames(pnonzero) = pnonzero[,1]; pnonzero[,1] = NULL
	totals = aggregate(t(mcounts), list(clusters[,1]), FUN=length); rownames(totals) = totals[,1]; totals[,1] = NULL
	mnonzero = mnonzero / totals
	pnonzero = pnonzero / totals
	mcounts = aggregate(t(mcounts), list(clusters[,1]), FUN=mean); rownames(mcounts) = mcounts[,1]; mcounts[,1] = NULL
	pcounts = aggregate(t(pcounts), list(clusters[,1]), FUN=mean); rownames(pcounts) = pcounts[,1]; pcounts[,1] = NULL
		
	if (postlog == TRUE) {
		mcounts = log2(mcounts + 1)
		pcounts = log2(pcounts + 1)
	}
	
	nx = nrow(mcounts)
	if (! is.null(xorder)) {
		xnames = xorder
	} else {
		xnames = rownames(mcounts)
	}
	
	mcounts = mcounts[xnames,]
	pcounts = pcounts[xnames,]
	if (! is.null(clusters)) {
		mnonzero = mnonzero[xnames,]
		pnonzero = pnonzero[xnames,]
	}
		
	# plot mat and pat average expression
	if (addsims == TRUE) {
		mydata = data.frame(samp = rep(xnames,3), parent = c(rep("Mat",nx),rep("Pat",nx),rep("symPat",nx)), counts = c(mcounts,pcounts,simpmeans), x=c(1:nx,1:nx,1:nx), stringsAsFactors=FALSE)
		if (is.null(ymax)) {
			ymax = max(mydata$counts)
		}
		if (is.null(ymin)) {
			ymin = min(mydata$counts)
		}

		if (! is.null(clusters)) {
			mydata$nonzerofrac = c(mnonzero,pnonzero,rep(0,length(pnonzero)))
		} else {
			mydata$nonzerofrac = rep(1,nrow(mydata))
		}
		if (is.null(smax)) smax = max(mydata$nonzerofrac)
	
		if (opt$skipplot == FALSE) {
			pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 8)
			pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 18)
			cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
			pdf(paste(outprefix,'_plot.pdf',sep=''), width = pwidth, height = pheight)
			p1 = ggplot(mydata, aes(x=x, y=counts)) + geom_line(size = linewidth, aes(color=parent, linetype = parent)) + geom_point(aes(size=nonzerofrac, fill=parent, color=parent),pch=21) + theme_bw()
			p1 = p1 + scale_colour_manual(values = c('#F8766D','#00BFC4','#00BFC4'))
			p1 = p1 + scale_fill_manual(values = c('#F8766D','#00BFC4','#00BFC4'))
			p1 = p1 + scale_linetype_manual(values = c(1,1,2))
			if (noprelog == FALSE && postlog == FALSE) p1 = p1 + labs(fill = "avg. counts", size = "frac. non-zero") + ylab("avg. counts")
			if (noprelog == TRUE && postlog == TRUE) p1 = p1 + labs(fill = "log2(avg. counts)", size = "frac. non-zero") + ylab("log2(avg. counts)")
			if (noprelog == TRUE && postlog == FALSE) p1 = p1 + labs(fill = "avg. counts", size = "frac. non-zero") + ylab("avg. counts")
			p1 = p1 + ylim(ymin,ymax)
			p1 = p1 + scale_x_continuous(breaks=1:nx, labels=xnames, limits = c(0.5,nx+0.5)) + theme(panel.grid.minor.x = element_blank())
			p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax))

			show(p1)
			graphics.off()
		}
	} else {
		mydata = data.frame(samp = rep(xnames,2), parent = c(rep("Mat",nx),rep("Pat",nx)), counts = c(mcounts,pcounts), x=c(1:nx,1:nx), stringsAsFactors=FALSE)
		if (is.null(ymax)) {
			ymax = max(mydata$counts)
		}
		if (is.null(ymin)) {
			ymin = min(mydata$counts)
		}

		if (! is.null(clusters)) {
			mydata$nonzerofrac = c(mnonzero,pnonzero)
		} else {
			mydata$nonzerofrac = rep(1,nrow(mydata))
		}
		smax = max(mydata$nonzerofrac)
	
		if (opt$skipplot == FALSE) {
			pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 8)
			pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 18)
			cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
			pdf(paste(outprefix,'_plot.pdf',sep=''), width = pwidth, height = pheight)
			p1 = ggplot(mydata, aes(x=x, y=counts)) + geom_line(size = linewidth, aes(color=parent)) + geom_point(aes(size=nonzerofrac, fill=parent, color=parent),pch=21) + theme_bw()
			if (noprelog == FALSE && postlog == FALSE) p1 = p1 + labs(fill = "avg. counts", size = "frac. non-zero") + ylab("avg. counts")
			if (noprelog == TRUE && postlog == TRUE) p1 = p1 + labs(fill = "log2(avg. counts)", size = "frac. non-zero") + ylab("log2(avg. counts)")
			if (noprelog == TRUE && postlog == FALSE) p1 = p1 + labs(fill = "avg. counts", size = "frac. non-zero") + ylab("avg. counts")
			p1 = p1 + ylim(ymin,ymax)
			p1 = p1 + scale_x_continuous(breaks=1:nx, labels=xnames, limits = c(0.5,nx+0.5)) + theme(panel.grid.minor.x = element_blank())
			p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax))

			show(p1)
			graphics.off()
		}
	}

	if (opt$outputdatalong == TRUE || opt$outputdatawide == TRUE) {
		write.table(mydata,file=paste(outprefix,"_data.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	}	
}


# Line plot of average total expression, over samples (or clusters if provided)
make_lin_plot_expr <- function(expr_matrix, outprefix, clusters, maxdims = 100, xorder = NULL, noprelog = FALSE, postlog = FALSE, includezeros = FALSE, linewidth = 2, color='black', expr_type = 'counts', ymax = NULL, ymin = NULL, smax = NULL) {

	# stop if xorder has repeats of same gene(s) or samples
	if (!is.null(xorder)) xorder = strsplit(xorder,",")[[1]]
	if (! is.null(xorder) && length(xorder) != length(unique(xorder))) stop("--xorder list cannot contain repeated values")

	# stop if more than maxdims genes/gene groups or maxdims samples/sample clusters requested, plot will be unintelligible
	if (is.null(clusters)) stop("Plots of type 'lin' currently only work for clusters of samples (e.g. --sampfile must have 2 columns; to cheat, make column 2 identical to column 1")
	cclusters = sort(unique(clusters[,1]))
	if (length(cclusters) > maxdims) stop("Error: dot plot would have ",length(cclusters)," columns (sample clusters), which is greater than the value of --maxdims ",maxdims,"and would probably create a very crowded plot. Please either reduce the number of sample clusters or increase --maxdims.")
	if (! is.null(xorder) && setequal(cclusters,xorder) == FALSE) stop("if using --xorder to order sample clusters, must have same values as column 2 of --sampfile")
	if (! is.null(xorder)) cclusters = xorder
	
	# take log unless not requested
	if (noprelog == FALSE) {
		expr_matrix = log2(expr_matrix + 1)
	}
				
	# get actual count distributions	
	clusters = clusters[colnames(expr_matrix),,drop=FALSE]
	
	nonzero = aggregate(t(expr_matrix), list(clusters[,1]), FUN= function(x) sum(x>0)); rownames(nonzero) = nonzero[,1]; nonzero[,1] = NULL
	totals = aggregate(t(expr_matrix), list(clusters[,1]), FUN=length); rownames(totals) = totals[,1]; totals[,1] = NULL
	frac_nonzero = unname(nonzero / totals)
	frac_nonzero = frac_nonzero[order(match(rownames(frac_nonzero), cclusters)), , drop = FALSE]
	frac_nonzero = as.matrix(frac_nonzero)	
	if (includezeros == TRUE) {
		meancounts = unname(aggregate(t(expr_matrix), list(clusters[,1]), FUN=mean)); rownames(meancounts) = meancounts[,1]; meancounts[,1] = NULL
		meancounts = meancounts[order(match(rownames(meancounts), cclusters)), , drop = FALSE]
		meancounts = as.matrix(meancounts)	
	} else {
		a = cbind(t(expr_matrix),clusters)
		a = a[a[,1] > 0,]
		meancounts = unname(aggregate(a[,1], list(a[,2]), FUN=mean)); rownames(meancounts) = meancounts[,1]; meancounts[,1] = NULL		
		cclustersdf = data.frame(clusterid = cclusters)

		meancounts = merge(x = cclustersdf, y = meancounts, by.x = 1, by.y = 0, all.x = TRUE)
		meancounts[is.na(meancounts)] = 0	
		rownames(meancounts) = meancounts$clusterid
		meancounts$clusterid = NULL
		colnames(meancounts) = NULL
		meancounts = meancounts[order(match(rownames(meancounts), cclusters)), , drop = FALSE]
		meancounts = as.matrix(meancounts)
	}
						
	if (postlog == TRUE) {
		meancounts = log2(meancounts + 1)
	}
				
	# plot mat and pat average expression
	mydata = data.frame(samp = rownames(meancounts), counts = meancounts, nonzerofrac = frac_nonzero, x=1:nrow(meancounts), stringsAsFactors=FALSE)
	
	if (is.null(smax)) smax = max(mydata$nonzerofrac)
	
	if (is.null(ymax)) {
		ymax = max(meancounts)
	}
	if (is.null(ymin)) {
		ymin = min(meancounts)
	}

	if (opt$skipplot == FALSE) {
		pheight = ifelse(! is.null(opt$plotheight), opt$plotheight, 8)
		pwidth = ifelse(! is.null(opt$plotwidth), opt$plotwidth, 18)
		cat('plot dimensions are',pwidth,'x',pheight,'(w x h)\n')
		pdf(paste(outprefix,'_plot.pdf',sep=''), width = pwidth, height = pheight)
		p1 = ggplot(mydata, aes(x=x, y=counts)) + geom_line(size = linewidth, color = color) + geom_point(aes(size=nonzerofrac),pch=21, color = color, fill=color) + theme_bw()
		p1 = p1 + labs(size = "frac. non-zero")
		if (expr_type == "CPM") {
			if (noprelog == FALSE && postlog == FALSE) p1 = p1 + ylab("avg. log2(CPM)")
			if (noprelog == TRUE && postlog == TRUE) p1 = p1 + ylab("log2(avg. CPM)")
			if (noprelog == TRUE && postlog == FALSE) p1 = p1 + ylab("avg. CPM")
		} else if (expr_type == "counts") {
			if (noprelog == FALSE && postlog == FALSE) p1 = p1 + ylab("avg. log2(norm. read counts)")
			if (noprelog == TRUE && postlog == TRUE) p1 = p1 + ylab("log2(avg. norm. read counts)")
			if (noprelog == TRUE && postlog == FALSE) p1 = p1 + ylab("avg. norm. read counts")
		} else {
			if (noprelog == FALSE && postlog == FALSE) p1 = p1 + ylab("avg. log2(expression)")
			if (noprelog == TRUE && postlog == TRUE) p1 = p1 + ylab("log2(avg. expression)")
			if (noprelog == TRUE && postlog == FALSE) p1 = p1 + ylab("avg. expression")
		} 
		p1 = p1 + ylim(ymin,ymax)
		p1 = p1 + scale_x_continuous(breaks=1:nrow(meancounts), labels=rownames(mydata), limits = c(0.5,nrow(meancounts)+0.5)) + theme(panel.grid.minor.x = element_blank())
		p1 = p1 + scale_size_continuous(range = c(0,opt$dotsize), limits=c(0,smax))

		show(p1)
		graphics.off()
	}

	if (opt$outputdatalong == TRUE || opt$outputdatawide == TRUE) {
		write.table(mydata,file=paste(outprefix,"_data.txt",sep=''),col.names = NA,row.names = T,quote = F,sep='\t')
	}	
}



# HOUSEKEEPING
# ---------------------
# This section checks that inputs are ok, and then outputs summary of params to user
cat("Checking all input files/params...\n")

# Check all inputs ok
if (opt$noprelog == FALSE && opt$postlog == TRUE) stop("can't take log2 both before and after taking mean (either drop --postlog or add --noprelog)")

# Check --genes or --genefile OK
if (is.null(opt$genes) && is.null(opt$genefile)) {
	stop("must provide gene(s) for plot using either --genes or --genefile")
}

if (! is.null(opt$genes) && ! is.null(opt$genefile)) {
	stop("must use either --genes or --genefile, not both")
}

if (! is.null(opt$genes) && is.null(opt$genefile)) {
	genelist = strsplit(opt$genes,",")[[1]]
	if (length(genelist) != length(unique(genelist))) stop("--genes list cannot contain repeated values")
	genegroups = NULL; ngroups = 0
	ngenes = length(genelist)
	genell = NULL
} else {
	if(isFALSE(file.exists(opt$genefile))) stop("cannot open --genefile file ",opt$genefile)
	genell = read.table(opt$genefile, row.names=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
	ngenes = nrow(genell)
	genelist = rownames(genell)
	if (length(genelist) != length(unique(genelist))) stop("--genes list cannot contain repeated values")

	if (ncol(genell) == 1) {
		colnames(genell) = c('groups')
		ngroups = length(unique(genell$groups))
	} else if (ncol(genell) == 0) {
		genell = NULL
	} else {
		stop("--genelist file contains ",ncol(genell)," columns, but should only have a maximum of two (gene IDs, and optionally gene groups, which are only used if plottype == 'dot').")
	}
}

# Check --sampfile OK
if (opt$plottype %in% c('bar','nuc','dot','lin') && ! is.null(opt$sampfile)) {
	if(isFALSE(file.exists(opt$sampfile))) stop("cannot open --sampfile file ",opt$sampfile)	
	clusters = read.table(opt$sampfile, row.names=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
	if (ncol(clusters) > 1) stop("--sampfile file should have 1 or 2 columns; file",opt$sampfile,"has",ncol(clusters)+1,"columns.")
	if (length(rownames(clusters)) != length(unique(rownames(clusters)))) stop("samples list provided in --sampfile cannot contain repeated values")
} else if (! is.null(opt$sampfile)) {
	cat("Warning: --sampfile file provided, but is not used for plot type",opt$plottype,"and will be ignored.\n")
} else if (is.null(opt$sampfile) && opt$plottype == 'lin') {
	stop("--sampfile is required for plot type 'lin'")
}

# Check input count or CPM file(s) OK
if (opt$plottype %in% c('gof','cmp','nuc')) {
	if (is.null(opt$mcounts)) stop("for plot type '",opt$plottype,"', the --mcounts argument is required.")
	if (is.null(opt$pcounts)) stop("for plot type '",opt$plottype,"', the --pcounts argument is required.")
	if(isFALSE(file.exists(opt$mcounts))) stop("cannot open --mcounts file ",opt$mcounts)
	if(isFALSE(file.exists(opt$pcounts))) stop("cannot open --pcounts file ",opt$pcounts)
	
	if (! is.null(opt$acounts)) cat("Warning: for plot type '",opt$plottype,"', the --acounts argument is ignored.\n")
	if (! is.null(opt$CPM)) cat("Warning: for plot type '",opt$plottype,"', the --CPM argument is ignored.\n")
} else if (opt$plottype %in% c('dot','bar','lin')) {
	if (is.null(opt$acounts) && is.null(opt$CPM) && is.null(opt$mcounts) && is.null(opt$pcounts)) stop("for plot type '",opt$plottype,"', either --acounts OR --CPM (expression plot) or --mcounts AND --pcounts (imprinting plot) must be provided.")
	
	if (! is.null(opt$acounts) || ! is.null(opt$CPM)) {
		if (! is.null(opt$mcounts) || ! is.null(opt$pcounts)) stop("for plot type '",opt$plottype,"', must provide either --acounts/--CPM OR --mcounts + --pcounts, not both.") 
		mmode = "expr"
		if (! is.null(opt$acounts) && is.null(opt$CPM)) {
			expr_datafile = opt$acounts; expr_type = "counts"
			if(isFALSE(file.exists(expr_datafile))) stop("cannot open --acounts file ",expr_datafile)	
		} else {
			expr_datafile = opt$CPM; expr_type = "CPM"
			if(isFALSE(file.exists(expr_datafile))) stop("cannot open --CPM file ",expr_datafile)	
		}
	} else {
		if (is.null(opt$mcounts)) stop("for plot type '",opt$plottype,"' in 'imprinting mode', both the --mcounts and --pcounts arguments are required.")
		if (is.null(opt$pcounts)) stop("for plot type '",opt$plottype,"' in 'imprinting mode', both the --mcounts and --pcounts arguments are required.")
		if(isFALSE(file.exists(opt$mcounts))) stop("cannot open --mcounts file ",opt$mcounts)
		if(isFALSE(file.exists(opt$pcounts))) stop("cannot open --pcounts file ",opt$pcounts)	
		mmode = "impr"			
	}
}	

# Output summary of params to user, for each plot type
if (opt$plottype == 'gof') {	
	if (is.null(opt$ASE_params)) stop("for plot type '",opt$plottype,"', the --ASE_params argument is required.")
	if(isFALSE(file.exists(opt$ASE_params))) stop("cannot open --ASE_params file ",opt$ASE_params)
	
	if (ngenes != 1) stop("for plot type '",opt$plottype,"', only one gene can be provided via --genes or --genefile. ",ngenes," genes provided by user, exiting.")
	gene = genelist[1]
	
	cat("\n")
	cat("Making plot showing goodness-of-fit of ZINB models, for gene:",gene,"\n")
	cat("-----------------------\n")
	cat("Current directory:",getwd(),"\n")
	cat("Maternal normalized counts:",opt$mcounts,"\n")
	cat("Paternal normalized counts:",opt$pcounts,"\n")
	cat("ASE analysis params file:",opt$ASE_params,"\n")
	if (isFALSE(opt$skipplot)) cat("Saving plot to:",paste(opt$outprefix,'_plot.pdf',sep=''),"\n")
	if (isTRUE(opt$outputdatalong) || isTRUE(opt$outputdatawide)) cat("Saving processed data for plot to:",paste(opt$outprefix,'_data.txt',sep=''),"\n")
	cat("-----------------------\n")
	if (! is.null(opt$xmax)) cat("Maximum value for x-axis:",opt$xmax,"\n")
	if (! is.null(opt$xmax)) cat("-----------------------\n")
} else if (opt$plottype == 'cmp') {
	if (ngenes != 1) stop("for plot type '",opt$plottype,"', only one gene can be provided via --genes or --genefile. ",ngenes," genes provided by user, exiting.")
	gene = genelist[1]
	
	cat("\n")
	cat("Comparing ZINB and NB distribution fits to allelic count data, for gene:",opt$genes,"\n")
	cat("-----------------------\n")
	cat("Current directory:",getwd(),"\n")
	cat("Maternal normalized counts:",opt$mcounts,"\n")
	cat("Paternal normalized counts:",opt$pcounts,"\n")
	if (isFALSE(opt$skipplot)) cat("Saving plot to:",paste(opt$outprefix,'_plot.pdf',sep=''),"\n")
	if (isTRUE(opt$outputdatalong) || isTRUE(opt$outputdatawide)) cat("Saving processed data for plot to:",paste(opt$outprefix,'_data.txt',sep=''),"\n")
	cat("-----------------------\n")
	if (! is.null(opt$xmax)) cat("Maximum value for x-axis:",opt$xmax,"\n")
	if (! is.null(opt$xmax)) cat("-----------------------\n")
} else if (opt$plottype == 'bar') {
	if (ngenes != 1) stop("for plot type '",opt$plottype,"', only one gene can be provided via --genes or --genefile. ",ngenes," genes provided by user, exiting.")
	gene = genelist[1]
	
	cat("\n")
	cat("Making plot of allelic counts per nucleus (bar chart), for gene:",opt$genes,"\n")
	cat("-----------------------\n")
	cat("Current directory:",getwd(),"\n")
	if (mmode == "impr") {
		cat("Maternal normalized counts:",opt$mcounts,"\n")
		cat("Paternal normalized counts:",opt$pcounts,"\n")
	} else {
		cat("Expression data (",expr_type,"):",expr_datafile,"\n")
	}
	if (isFALSE(opt$skipplot)) cat("Saving plot to:",paste(opt$outprefix,'_plot.pdf',sep=''),"\n")
	if (isTRUE(opt$outputdatalong) || isTRUE(opt$outputdatawide)) cat("Saving processed data for plot to:",paste(opt$outprefix,'_data.txt',sep=''),"\n")
	cat("-----------------------\n")
	if (! is.null(opt$sampfile)) cat("Nuclei samples/clusters provided in file:",opt$sampfile,"\n")
	if (! is.null(opt$sampfile)) cat("-----------------------\n")
} else if (opt$plottype == 'nuc') {
	cat("\n")
	cat("Making plot of allelic counts per nucleus (nuc plot) for",ngenes,"genes\n")
	cat("-----------------------\n")
	cat("Current directory:",getwd(),"\n")
	cat("Maternal normalized counts:",opt$mcounts,"\n")
	cat("Paternal normalized counts:",opt$pcounts,"\n")
	if (isFALSE(opt$skipplot)) cat("Saving plot to:",paste(opt$outprefix,'_plot.pdf',sep=''),"\n")
	if (isTRUE(opt$outputdatalong) || isTRUE(opt$outputdatawide)) cat("Saving processed data for plot to:",paste(opt$outprefix,'_data.txt',sep=''),"\n")
	cat("-----------------------\n")
	if (! is.null(opt$sampfile)) cat("Nuclei samples/clusters provided in file:",opt$sampfile,"\n")
	if (! is.null(opt$sampfile)) cat("-----------------------\n")
} else if (opt$plottype == 'dot') {
	cat("\n")
	if (mmode == "expr") {
		if (is.null(genell)) cat("Making dot plot of expression over different nuclei or nuclei clusters for",ngenes,"genes\n")
		if (! is.null(genell)) cat("Making dot plot of expression over different nuclei or nuclei clusters for",ngroups,"gene groups\n")
		cat("-----------------------\n")
		cat("Current directory:",getwd(),"\n")
		cat("Expression data (",expr_type,"):",expr_datafile,"\n")
		if (isFALSE(opt$skipplot)) cat("Saving plot to:",paste(opt$outprefix,'_plot.pdf',sep=''),"\n")
		if (isTRUE(opt$outputdatalong) || isTRUE(opt$outputdatawide)) cat("Saving processed data for plot to:",paste(opt$outprefix,'_data.txt',sep=''),"\n")
	} else {
		if (is.null(genell)) cat("Making dot plot of allelic bias over different nuclei or nuclei clusters for",ngenes,"genes\n")
		if (! is.null(genell)) cat("Making dot plot of allelic bias over different nuclei or nuclei clusters for",ngroups,"gene groups\n")
		cat("-----------------------\n")
		cat("Current directory:",getwd(),"\n")
		cat("Maternal normalized counts:",opt$mcounts,"\n")
		cat("Paternal normalized counts:",opt$pcounts,"\n")
		if (isFALSE(opt$skipplot)) cat("Saving plot to:",paste(opt$outprefix,'_plot.pdf',sep=''),"\n")
		if (isTRUE(opt$outputdatalong) || isTRUE(opt$outputdatawide)) cat("Saving processed data for plot to:",paste(opt$outprefix,'_data.txt',sep=''),"\n")
	}
	cat("-----------------------\n")
	if (! is.null(opt$sampfile)) cat("Nuclei samples/clusters provided in file:",opt$sampfile,"\n")
	if (! is.null(opt$sampfile)) cat("-----------------------\n")
} else if (opt$plottype == 'lin') {
	if (ngenes != 1) stop("for plot type '",opt$plottype,"', only one gene can be provided via --genes or --genefile. ",ngenes," genes provided by user, exiting.")
	gene = genelist[1]

	cat("\n")
	if (mmode == "expr") {
		cat("Making connected line plot of expression over different nuclei or nuclei clusters for gene:",opt$genes,"\n")
		cat("-----------------------\n")
		cat("Current directory:",getwd(),"\n")
		cat("Expression data (",expr_type,"):",expr_datafile,"\n")
		if (isFALSE(opt$skipplot)) cat("Saving plot to:",paste(opt$outprefix,'_plot.pdf',sep=''),"\n")
		if (isTRUE(opt$outputdatalong) || isTRUE(opt$outputdatawide)) cat("Saving processed data for plot to:",paste(opt$outprefix,'_data.txt',sep=''),"\n")
	} else {
		cat("Making connected line plot of allelic bias over different nuclei or nuclei clusters for gene:",opt$genes,"\n")
		cat("-----------------------\n")
		cat("Current directory:",getwd(),"\n")
		cat("Maternal normalized counts:",opt$mcounts,"\n")
		cat("Paternal normalized counts:",opt$pcounts,"\n")
		if (isFALSE(opt$skipplot)) cat("Saving plot to:",paste(opt$outprefix,'_plot.pdf',sep=''),"\n")
		if (isTRUE(opt$outputdatalong) || isTRUE(opt$outputdatawide)) cat("Saving processed data for plot to:",paste(opt$outprefix,'_data.txt',sep=''),"\n")
	}
	cat("-----------------------\n")
	if (! is.null(opt$sampfile)) cat("Nuclei samples/clusters provided in file:",opt$sampfile,"\n")
	if (! is.null(opt$sampfile)) cat("-----------------------\n")
} else {
	stop("\nPlot type \"",opt$plottype,"\" is not currently supported. Supported plot types are 'gof', 'bar', 'nuc', and 'dot'.")
}

# load all libraries
cat("Loading all required libraries...\n")
for (lib in liblist) {
	suppressPackageStartupMessages(library(lib, character.only = TRUE))
}


# MAIN CODE
# ---------------------

# read in input files and check ok
cat("Reading in input files...\n")
if (opt$plottype %in% c('gof','nuc','cmp') || ((opt$plottype %in% c('dot','bar','lin') && mmode == 'impr'))) {
	mcounts_norm = data.matrix(read.table(opt$mcounts, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))
	pcounts_norm = data.matrix(read.table(opt$pcounts, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))

	# make sure they have same dimensions and contain the same values
	if (all(dim(mcounts_norm) == dim(pcounts_norm)) == FALSE) stop("Input maternal and paternal count matrices didn't have the same dimensions (mcounts has dimensions",dim(mcounts_norm),"while pcounts has dimensions",dim(pcounts_norm),")\n")
	if (setequal(rownames(mcounts_norm),rownames(pcounts_norm)) == FALSE) stop("Row names not equal in pcounts and mcounts")
	if (setequal(colnames(mcounts_norm),colnames(pcounts_norm)) == FALSE) stop("Column names not equal in pcounts and mcounts")

	# make sure no duplicate rows or columns
	if (length(rownames(mcounts_norm)) != length(unique(rownames(mcounts_norm)))) stop("input --mcounts file contains the same gene ID more than once in the first column, please check your inputs.")
	if (length(colnames(mcounts_norm)) != length(unique(colnames(mcounts_norm)))) stop("input --mcounts file contains the same sample ID more than once in the first row, please check your inputs.")

	# sort pcounts so that row and column order match mcounts
	pcounts_norm = pcounts_norm[match(rownames(mcounts_norm), rownames(pcounts_norm)),]
	pcounts_norm = pcounts_norm[,match(colnames(mcounts_norm), colnames(pcounts_norm))]

	cat("Input count matrices have",nrow(mcounts_norm),"rows (genes) x",ncol(mcounts_norm),"columns (samples/cells)\n")
} else {
	expr_data = data.matrix(read.table(expr_datafile, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE))	
	if (expr_type == "CPM") cat("Input CPM matrix has",nrow(expr_data),"rows (genes) x",ncol(expr_data),"columns (samples/cells)\n")
	if (expr_type == "counts") cat("Input count matrix has",nrow(expr_data),"rows (genes) x",ncol(expr_data),"columns (samples/cells)\n")
}

if (opt$plottype == 'gof' || opt$plottype == 'lin') {
	if (! is.null(opt$ASE_params)) {
		ASE_params = read.table(opt$ASE_params, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)	
		coef_nu = ASE_params[which(ASE_params[,1] == "coef_nu"),2]
		coef_mu_ZINB = ASE_params[which(ASE_params[,1] == "coef_mu_ZINB"),2]
		coef_mu_NB = ASE_params[which(ASE_params[,1] == "coef_mu_NB"),2]

		if (length(coef_nu) == 0) stop("--ASE_params file ",opt$ASE_params," did not contain required field 'coef_nu'. Make sure you're using the file output by single_cell_ASE_analysis.R.")
		if (length(coef_mu_ZINB) == 0) stop("--ASE_params file ",opt$ASE_params," did not contain required field 'coef_mu_ZINB'. Make sure you're using the file output by single_cell_ASE_analysis.R.")
		if (length(coef_mu_NB) == 0) stop("--ASE_params file ",opt$ASE_params," did not contain required field 'coef_mu_NB'. Make sure you're using the file output by single_cell_ASE_analysis.R.")	
	} else {
		coef_nu = opt$coef_nu
		coef_mu_ZINB = opt$coef_mu_ZINB
		coef_mu_NB = opt$coef_mu_NB
	}
}

if (! is.null(opt$sampfile)) {
	# ensure all nuclei in this file are present in the count matrices/expr_data, and subset accordingly
	if (opt$plottype %in% c('gof','nuc','cmp') || ((opt$plottype %in% c('dot','bar','lin') && mmode == 'impr'))) {
		mcounts_norm = mcounts_norm[,colnames(mcounts_norm) %in% rownames(clusters)]
		if (ncol(mcounts_norm) < nrow(clusters)) stop("not all samples in --sampfile were found in input count matrices, aborting (does your file have a header? if so, removing it may fix this error)")
		pcounts_norm = pcounts_norm[,colnames(pcounts_norm) %in% rownames(clusters)]
		cat("Retaining only the",ncol(mcounts_norm),"samples present in --sampfile\n")
	} else {
		expr_data = expr_data[,colnames(expr_data) %in% rownames(clusters)]
		if (ncol(expr_data) < nrow(clusters)) stop("not all samples in --sampfile were found in input matrix, aborting (does your file have a header? if so, removing it may fix this error)")
		cat("Retaining only the",ncol(expr_data),"samples present in --sampfile\n")
	}
	if (ncol(clusters) == 1) colnames(clusters) = c("cluster")
	if (ncol(clusters) == 0) clusters = NULL
} else {
	clusters = NULL
}

# get data for gene(s) to use for plot
if (opt$plottype %in% c('gof','nuc','cmp') || ((opt$plottype %in% c('dot','bar','lin') && mmode == 'impr'))) {
	i = which(rownames(mcounts_norm) %in% genelist)
} else {
	i = which(rownames(expr_data) %in% genelist)
}

if (length(i) != length(genelist)) {
	missing_genes = 0; genelist_touse = c()

	for (gene in genelist) {
		if (opt$plottype %in% c('gof','nuc','cmp') || ((opt$plottype %in% c('dot','bar','lin') && mmode == 'impr'))) {
			j = which(rownames(mcounts_norm) == gene)
		} else {
			j = which(rownames(expr_data) == gene)
		}
		if (length(j) == 0) {
			if (opt$allowmissinggenes == FALSE) {
				stop("Gene ",gene," not found in input count or CPM matrices")
			} else {
				missing_genes = missing_genes + 1
			}
		} else {
			genelist_touse = c(genelist_touse, gene)
		}
	}	
	if (opt$allowmissinggenes == TRUE) cat("Warning: a total of",missing_genes,"genes present in --genes/--genefile were not found in matrix and have been discarded.\n")
	genelist = genelist_touse
	if (! is.null(genell)) genell = subset(genell, rownames(genell) %in% genelist)
}

if (opt$plottype %in% c('gof','nuc','cmp') || ((opt$plottype %in% c('dot','bar','lin') && mmode == 'impr'))) {
	mcounts = mcounts_norm[i,, drop = FALSE]
	if (ngenes > 1) mcounts = mcounts[match(genelist, rownames(mcounts)),]
	i = which(rownames(pcounts_norm) %in% genelist)
	pcounts = pcounts_norm[i,, drop = FALSE]
	if (ngenes > 1) pcounts = pcounts[match(genelist, rownames(pcounts)),]
	
	if (all(is.na(mcounts))) {
		stop("Gene(s) provided are present in the input maternal count matrix, but have all missing values. Aborting.")
	}
	if (all(is.na(pcounts))) {
		stop("Gene(s) provided are present in the input paternal count matrix, but have all missing values. Aborting.")
	}
} else {
	exprs = expr_data[i,, drop = FALSE]
	if (ngenes > 1) exprs = exprs[match(genelist, rownames(exprs)),]
	if (all(is.na(exprs))) {
		stop("Gene(s) provided are present in the input expression matrix, but have all missing values. Aborting.")
	}
}

# make plot
cat("Making plot...")
if (opt$plottype == "gof") {
	make_gof_plot(mcounts, pcounts, opt$outprefix, coef_nu, coef_mu_ZINB, coef_mu_NB, xmax = opt$xmax)
} else if (opt$plottype == 'cmp') {
	make_cmp_plot(mcounts, pcounts, opt$outprefix, xmax = opt$xmax)
} else if (opt$plottype == "bar") {
	if (mmode == "expr") {
		make_bar_plot_expr(exprs, opt$outprefix, clusters, xorder = opt$xorder)
	} else {
		make_bar_plot_impr(mcounts, pcounts, opt$outprefix, clusters, xorder = opt$xorder)
	}				
} else if (opt$plottype == "nuc") {
	make_nuc_plot(mcounts, pcounts, opt$outprefix, clusters, opt$minreads)
} else if (opt$plottype == "dot") {
	if (mmode == "expr") {
		make_dot_plot_expr(exprs, opt$outprefix, clusters = clusters, genegroups = genell, maxdims = opt$maxdims, xorder = opt$xorder, yorder = opt$yorder, expr_type = expr_type, noprelog = opt$noprelog, includezeros = opt$includezeros, fmax = opt$fillupper, smax = opt$sizeupper)
	} else {
		make_dot_plot_impr(mcounts, pcounts, opt$outprefix, clusters = clusters, genegroups = genell, minreads = opt$minreads, maxdims = opt$maxdims, xorder = opt$xorder, yorder = opt$yorder)
	}
} else if (opt$plottype == "lin") {
	if (mmode == "expr") {
		make_lin_plot_expr(exprs, opt$outprefix, clusters = clusters, maxdims = opt$maxdims, xorder = opt$xorder, expr_type = expr_type, postlog = opt$postlog, noprelog = opt$noprelog, includezeros = opt$includezeros, linewidth = opt$linewidth, color = 'grey40', ymin = opt$ymin, ymax = opt$ymax, smax = opt$sizeupper)
	} else {
		make_lin_plot_impr(mcounts, pcounts, opt$outprefix, clusters, coef_nu = coef_nu, coef_mu_ZINB = coef_mu_ZINB, coef_mu_NB = coef_mu_NB, maxdims = opt$maxdims, xorder = opt$xorder, postlog = opt$postlog, noprelog = opt$noprelog, linewidth = opt$linewidth, nreps = opt$nreps, addsims = opt$addsims, ymin = opt$ymin, ymax = opt$ymax, smax = opt$sizeupper)
	}
}

cat("DONE\n")





