#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(gamlss))
suppressPackageStartupMessages(library(maxLik))

# version 1.0 (07/31/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 07/31/2019
# -------------------------

# Description:
# -------------------------
# This script simulates single-cell counts from endosperm, where there are two maternal and
# one paternal genomes. 

# Maternal counts are simulated as m ~ ZINB() + ZINB() (sum of two ZINB distributions)
# while paternal counts are simulated as p ~ ZINB() (single ZINB distribution).

# ZINB is parameterized by nu, mu and sigma, where nu is the probability of a 'drop-out' zero,
# mu is the average of the nonzero counts, and sigma is a variance/dispersion parameter.

# User specifies desired average total expression, sigma, log2fc(m/p), and number of observations
# per condition. Can specify multiple values for each (comma-separated); script will do all
# possible combinations.

# User must also provide either mu or nu; other param will be calculated by script to provide
# desired log2fc and total mean expression. For a given total mean_expr m and log2fc r, 
# script will calculate:

# mean_m = (2^r + 1) / (2^r * m)
# mean_p = m - mean_m

# Since:

# mean_m = mu_m*(1-nu_m)
# mean_p = mu_p*(1-nu_p)

# as long as either nu or mu are known, the other parameter value can be obtained that will
# satisfy m and r. User must provide mu_m or nu_m separately from mu_p or nu_p, though the
# values can be the same. Must provide same number of mu_m and mu_p (or nu_m and nu_p) values,
# script will iterate through both in parallel instead of doing all combinations.

# Note that if more than one nu/mu provided, these are looped through in parallel with 
# mean_expr, and so there must be as many as mean_expr

# Notes:
# -------------------------
# TODO

# Usage: single_cell_ASE_analysis.R [options] expr_matrix.txt outprefix

parser = ArgumentParser()
parser$add_argument("dataset", nargs=1, help = "dataset with true values of mean_m, mu_m, etc. (output *_fits.txt from single_cell_ASE_analysis)")
parser$add_argument("outfile", nargs=1, help = "prefix for output files")
parser$add_argument("--seed", default = 123456, help = "Set random number generator seed; runs with same seed and other inputs should always produce same results.", type="integer")
parser$add_argument("--n_reps", default = 100, help = "Number of iterations per condition (combination of params below)", type="integer")
parser$add_argument("--ncells", default = 100, help = "Number of cells/nuclei to simulate.", type="integer")
parser$add_argument("--log2fc", default = 1, help = "Target log2(m/p) values for simulations", type="character")
parser$add_argument("--mean_expr", default = 1, help = "Target average -total- expression levels (maternal + paternal)", type="character")
parser$add_argument("--nu_m", help = "Override estimate of nu_m based on dataset and use this value instead; one value per value in --mean_expr (but can provide same value to --mean_expr multiple times) (not used if --fix_mu)", type="character")
parser$add_argument("--nu_p", help = "Override estimate of nu_p based on dataset and use this value instead; one value per value in --mean_expr (but can provide same value to --mean_expr multiple times) (only used if --fix_nu)", type="character")
parser$add_argument("--mu_m", help = "Override estimate of mu_m based on dataset and use this value instead; one value per value in --mean_expr (but can provide same value to --mean_expr multiple times) (not used if --fix_nu or no fix)", type="character")
parser$add_argument("--fix_nu", default=FALSE, action="store_true", help = "Use this so that simulated maternal and paternal reads will share nu (takes coef_nu into account)")
parser$add_argument("--fix_mu", default=FALSE, action="store_true", help = "Use this so that simulated maternal and paternal reads will share mu (takes coef_nu into account)")
parser$add_argument("--fix_sigma", default=FALSE, action="store_true", help = "Use this so that simulated maternal and paternal reads will share sigma (will use maternal sigma)")
parser$add_argument("--sigma_pctls", default = '0.5', help = "Use value(s) of sigma corresponding to these pctls (0.5 = median) of values in the real dataset (based on obs with similar mean expression to desired value); default median", type="character")
parser$add_argument("--sigma_m", help = "Use these value(s) for sigma_m (must also provide sigma_p) (overrides --sigma_pctls))", type="character")
parser$add_argument("--sigma_p", help = "Use these value(s) for sigma_p (must also provide sigma_m) (overrides --sigma_pctls))", type="character")
parser$add_argument("--coef_nu", default = 1, help = "Coefficient describing relationship between mat and pat nu params", type="double")
parser$add_argument("--coef_mu_ZINB", default = 1, help = "Coefficient describing relationship between mat and pat mu params in ZINB", type="double")

opt <- parser$parse_args()

dataset = read.table(opt$dataset, header=TRUE, row.names=1, sep="\t", stringsAsFactors = FALSE)

set.seed(opt$seed)
n_reps = opt$n_reps
n_obs = opt$ncells
log2fc = as.numeric(strsplit(opt$log2fc,',')[[1]])
mean_expr = as.numeric(strsplit(opt$mean_expr,',')[[1]])
fix_nu = opt$fix_nu; fix_mu = opt$fix_mu
if (fix_nu == TRUE & fix_mu == TRUE) {
	stop("Error, you can't fix both mu and nu (nothing left to vary to hit log2fc and expression targets)")
}

coef_nu = opt$coef_nu
coef_mu_ZINB = opt$coef_mu_ZINB
outfile = opt$outfile

if (! is.null(opt$nu_m)) {
	nu_m_override = as.numeric(strsplit(opt$nu_m,',')[[1]])
	if (length(nu_m_override) != length(mean_expr)) { stop("Error: must provide as many values to --nu_m as --mean_expr") }
}

if (! is.null(opt$nu_p)) {
	nu_p_override = as.numeric(strsplit(opt$nu_p,',')[[1]])
	if (length(nu_p_override) != length(mean_expr)) { stop("Error: must provide as many values to --nu_p as --mean_expr") }
}

if (! is.null(opt$mu_m)) {
	mu_m_override = as.numeric(strsplit(opt$mu_m,',')[[1]])
	if (length(mu_m_override) != length(mean_expr)) { stop("Error: must provide as many values to --mu_m as --mean_expr") }
}

if (! is.null(opt$sigma_m) | ! is.null(opt$sigma_p)) {
	if (is.null(opt$sigma_m)) { stop("If providing --sigma_p, must also provide --sigma_m") }
	if (is.null(opt$sigma_p)) { stop("If providing --sigma_m, must also provide --sigma_p") }
	sigma_m_override = as.numeric(strsplit(opt$sigma_m,',')[[1]])
	sigma_p_override = as.numeric(strsplit(opt$sigma_p,',')[[1]])
	if (length(sigma_m_override) != length(sigma_p_override)) {
		stop("Must provide same number of values to --sigma_m and --sigma_p")
	}
	num_sigmas = length(sigma_m_override)
} else {
	sigma_pctls = as.numeric(strsplit(opt$sigma_pctls,',')[[1]])
	num_sigmas = length(sigma_pctls)
}

args <- commandArgs(trailingOnly = TRUE)	
if (length(opt) <= 1) {
	cat("Usage: single_cell_simulate_reads.R [options] outprefix\n")
	cat("----------------------\n")
	cat("outfile : name for output files\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("TODO\n")
	cat("----------------------\n")
}

cat("\nRunning single_cell_simulate_reads.R v.1.0 (08/25/2019)\n")
cat("-----------------------\n")
cat("Current directory:",getwd(),"\n")
cat("Name for output file:",outfile,"\n")
cat("-----------------------\n")
cat("Additional options:\n")
cat("Number of iterations per single set of conditions below:",n_reps,"\n")
cat("Number of cells to simulate:",opt$ncells,"\n")
cat("Target total expression level:",opt$mean_expr,"\n")
cat("Target log2(m/p) ratio:",opt$log2fc,"\n")
if (! is.null(opt$mu_m)) {
	cat("Value(s) of mu_m:",opt$mu_m,"\n")
}
if (! is.null(opt$nu_m)) {
	cat("Value(s) of nu_m:",opt$nu_m,"\n")
}
if (! is.null(opt$nu_p)) {
	cat("Value(s) of nu_p:",opt$nu_p,"\n")
}
if (! is.null(opt$sigma_m) | ! is.null(opt$sigma_p)) {
	cat("Value(s) of sigma_m:",opt$sigma_m,"\n")
	cat("Value(s) of sigma_p:",opt$sigma_p,"\n")
} else {
	cat("Using sigma values corresponding to these percentiles:",opt$sigma_pctls,"\n")
}
cat("coef_nu:",opt$coef_nu,"\n")
cat("coef_mu_ZINB:",opt$coef_mu_ZINB,"\n")
cat("-----------------------\n")


# FUNCTIONS
# ----------------
# All functions now live in single_cell_ASE_src.R
source("single_cell_ASE_src.R")


# MAIN
# ----------------------

# data frames to store results of simulations
totreps = length(log2fc) * length(mean_expr) * num_sigmas * n_reps
sim_params <- matrix(data=NA, nrow = totreps, ncol = 11)
sim_params <- as.data.frame(sim_params, stringsAsFactors = FALSE)
colnames(sim_params) = c("tot_mean","exp_mean_m","exp_mean_p","log2fc_mp","true_mu_m","true_nu_m","true_sigma_m","true_mu_p","true_nu_p","true_sigma_p","rep")

sim_results <- matrix(data=NA, nrow = totreps, ncol = 35)
sim_results <- as.data.frame(sim_results, stringsAsFactors = FALSE)	
colnames(sim_results) = c("fracnonzero_m","median_nonzero_m","mean_m","fracnonzero_p","median_nonzero_p","mean_p","log2fc_m_over_p","mu_m","nu_m","sigma_m","fit_m","mu_p","nu_p","sigma_p","fit_p","fit_H0","nu_H0","mu_H0","sigma_H0","logL_H1","logL_H0","pval","nu_H0_mod","mu_H0_mod","sigma_H0_mod","logL_H0_mod","pval_mod","logL_H0_2","pval_H0_2","logL_H0_2_mod","pval_H0_2_mod","logL_H0_3","pval_H0_3","logL_H0_3_mod","pval_H0_3_mod")

# cache and sort dataset so that we can quickly look up appropriate nu values
ss_m = dataset[dataset$fit_m == "ZINB",names(dataset) %in% c("mean_m","mu_m","nu_m","sigma_m")] 
ss_m = ss_m[order(ss_m$mean_m),]
ss_p = dataset[dataset$fit_p == "ZINB",names(dataset) %in% c("mean_p","mu_p","nu_p","sigma_p")] 
ss_p = ss_p[order(ss_p$mean_p),]

pseudocount = 1/n_obs
searchlen = floor(max(20,dim(dataset)[1] / 500))
coef_mu_NB = coef_mu_ZINB

cc = 1
#options(show.error.messages = FALSE)		# note - turn this off for debugging, it's annoying
for (r in log2fc) {
	# calculate value of d = multiplier in mu_p = d*coef_mu_ZINB*mu_m, (1-nu_p)=d*coef_nu*(1-nu_m)
	if (fix_nu == TRUE | fix_mu == TRUE) {
		d = (2^(1-r))/(coef_mu_ZINB*coef_nu)
	} else {
		d = sqrt((2^(1-r))/(coef_nu*coef_mu_ZINB))
	}

	for (i in 1:length(mean_expr)) {
	
		m = mean_expr[i]
	
		# get expected mean_m and mean_p based on these values of total expression and log2fc
		true_mean_m = ((2^r)*m)/((2^r)+1)
		true_mean_p = m - true_mean_m
		
		# to get single unique combination of params nu_m, nu_p, mu_m, mu_p, just need one
		# of the four - get mu_m as average of entries in example dataset that have mean_m
		# close to desired true_mean_m
		idx_start = which.min(abs(ss_m$mean_m - true_mean_m))		# find values in dataset with similar mean
		idx_end = tail(which(abs(ss_m$mean_m - true_mean_m) == abs(ss_m$mean_m[idx_start] - true_mean_m)),1)
		
		idx_start = idx_start - (searchlen / 2)
		idx_end = idx_end + (searchlen / 2)
		
		if (idx_end - idx_start > searchlen) {
			dif = floor((searchlen - (idx_end - idx_start))/2)
			idx_start = max(1,idx_start-dif); idx_end = min(idx_end+dif,length(ss_m$nu_m))
		}
		mu_m_estim = mean(ss_m[max(1,idx_start):min(idx_end,length(ss_m$mu_m)),]$mu_m)
		nu_m_estim = mean(ss_m[max(1,idx_start):min(idx_end,length(ss_m$nu_m)),]$nu_m)
		if (is.null(opt$sigma_m)) {
			sigma_m_touse = quantile(ss_m[max(1,idx_start):min(idx_end,length(ss_m$mu_m)),]$sigma_m,sigma_pctls)
		} else {
			sigma_m_touse = sigma_m_override
		}
				
		idx_start = which.min(abs(ss_p$mean_p - true_mean_p))		# find values in dataset with similar mean
		idx_end = tail(which(abs(ss_p$mean_p - true_mean_p) == abs(ss_p$mean_p[idx_start] - true_mean_p)),1)
		
		idx_start = idx_start - (searchlen / 2)
		idx_end = idx_end + (searchlen / 2)

		if (idx_end - idx_start > searchlen) {
			dif = floor((searchlen - (idx_end - idx_start))/2)
			idx_start = max(1,idx_start-dif); idx_end = min(idx_end+dif,length(ss_p$nu_p))
		}
		nu_p_estim = mean(ss_p[max(1,idx_start):min(idx_end,length(ss_p$nu_p)),]$nu_p)
		if (opt$fix_sigma == TRUE) {
			sigma_p_touse = sigma_m_touse
		} else {
			if (is.null(opt$sigma_p)) {
				sigma_p_touse = quantile(ss_p[max(1,idx_start):min(idx_end,length(ss_p$mu_p)),]$sigma_p,sigma_pctls)
			} else {
				sigma_p_touse = sigma_p_override
			}
		}
		
		# override with user-provided values if provided
		if (! is.null(opt$nu_m)) {
			nu_m_estim = nu_m_override[i]
		}
		if (! is.null(opt$mu_m)) {
			mu_m_estim = mu_m_override[i]
		}
		if (! is.null(opt$nu_p)) {
			nu_p_estim = nu_p_override[i]
		}
				
		# get all other params based on one fixed from above (depending on fix_nu/mu)
		if (fix_nu == TRUE) {
			nu_m = max(0.00001,min(1-0.00001,(nu_m_estim + nu_p_estim)/2))
			mu_m = max(0.00001,true_mean_m / (2*(1-nu_m)))
			mu_p = max(0.00001,d*coef_mu_ZINB*mu_m)
			nu_p = max(0.00001,min(1-0.00001,1-coef_nu*(1-nu_m)))
		} else if (fix_mu == TRUE) {
			mu_m = max(0.00001,mu_m_estim)
			nu_m = max(0.00001,min(1-0.00001,1-(true_mean_m / (2*mu_m))))
			mu_p = max(0.00001,coef_mu_ZINB*mu_m)
			nu_p = max(0.00001,min(1-0.00001,1-(d*coef_nu*(1-nu_m))))
		} else {
			nu_m = max(0.00001,min(1-0.00001,nu_m_estim))
			mu_m = max(0.00001,true_mean_m / (2*(1-nu_m)))
			mu_p = max(0.00001,d*coef_mu_ZINB*mu_m)
			nu_p = max(0.00001,min(1-0.00001,1-(d*coef_nu*(1-nu_m))))
		}
		
		# OK, now run simulations
		for (k in 1:num_sigmas) {
			sigma_m = sigma_m_touse[k]; sigma_p = sigma_p_touse[k]			
			for (n in 1:n_reps) {
				cat("\r",paste0("Running simulation ",cc," of ",totreps))

				sim_params[cc,] = c(m, true_mean_m, true_mean_p, r, mu_m, nu_m, sigma_m, mu_p, nu_p, sigma_p, n)
		
				sim_m = rZINBI(n_obs, mu_m, sigma_m, nu_m) + rZINBI(n_obs, mu_m, sigma_m, nu_m)
				sim_p = rZINBI(n_obs, mu_p, sigma_p, nu_p)
				
				# re-draw until we get something better than all-zero counts
				while (mean(c(sim_m,sim_p)) == 0) {
					sim_m = rZINBI(n_obs, mu_m, sigma_m, nu_m) + rZINBI(n_obs, mu_m, sigma_m, nu_m)
					sim_p = rZINBI(n_obs, mu_p, sigma_p, nu_p)
				}				
				
				# add pseudocount to largest value in each set of counts
				idx = which.max(sim_m)
				sim_m[idx] = sim_m[idx] + 1
				idx = which.max(sim_p)
				sim_p[idx] = sim_p[idx] + 1
			
				# first get mean, etc. of this new random draw
				sim_results$fracnonzero_m[cc] = length(sim_m[sim_m != 0]) / n_obs
				sim_results$median_nonzero_m[cc] = median(sim_m[sim_m != 0])
				sim_results$mean_m[cc] = mean(sim_m)
				sim_results$fracnonzero_p[cc] = length(sim_p[sim_p != 0]) / n_obs
				sim_results$median_nonzero_p[cc] = median(sim_p[sim_p != 0])
				sim_results$mean_p[cc] = mean(sim_p)
				sim_results$log2fc_m_over_p[cc] = log2(sim_results$mean_m[cc] / sim_results$mean_p[cc])
				
				# estimate params and best fit
				res_m = fitdist(sim_m, "mother")
				sim_results$nu_m[cc] = res_m[[1]]
				sim_results$mu_m[cc] = res_m[[2]]
				sim_results$sigma_m[cc] = res_m[[3]]
				sim_results$fit_m[cc] = res_m[[4]]
				logL_m = res_m[[5]]
		
				res_p = fitdist(sim_p, "father")
				sim_results$nu_p[cc] = res_p[[1]]
				sim_results$mu_p[cc] = res_p[[2]]
				sim_results$sigma_p[cc] = res_p[[3]]
				sim_results$fit_p[cc] = res_p[[4]]
				logL_p = res_p[[5]]
				
				# do LRTs
				res = testH0(sim_m, sim_p, sim_results[cc,], coef_nu = coef_nu, coef_mu_ZINB = coef_mu_ZINB, coef_mu_NB = coef_mu_NB)
				sim_results[cc] = res
				cc = cc + 1		
			}
		}
	}
}
options(show.error.messages = TRUE)

# append the parameters used to each df
sim_fullres = merge(sim_params, sim_results, by=0); sim_fullres = sim_fullres[order(as.numeric(sim_fullres$Row.names)),]; sim_fullres$Row.names = NULL; rownames(sim_fullres) = 1:length(rownames(sim_fullres))

# output all results
write.table(sim_fullres,file=outfile,col.names = NA,row.names = T,quote = F,sep='\t')

