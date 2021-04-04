#!/usr/bin/env Rscript

# version 1.0 (09/18/2019)
# -------------------------
# Version history:
# v.1.0: initial build - 09/18/2019
# -------------------------

# Description:
# -------------------------
# Contains functions used by single_cell_ASE_analysis.R, single_cell_simulate_reads.R
# and single_cell_ASE_bycluster.R.

# normcounts uses code from DESingle (credit to Zhun Miao et al.) to normalize all
# three count matrices (acounts = all counts, mcounts = maternal counts, pcounts = paternal counts)
# based on expression in acounts
normCounts <- function(acounts, mcounts, pcounts) {
	geneNum <- nrow(acounts)
	nObs <- ncol(acounts)

	GEOmean <- rep(NA,geneNum)
	for (i in 1:geneNum)
	{
		gene_NZ <- acounts[i,acounts[i,] > 0]
		GEOmean[i] <- exp(sum(log(gene_NZ), na.rm=TRUE) / length(gene_NZ))
	}

	S <- rep(NA, nObs)
	acounts_norm <- acounts
	mcounts_norm <- mcounts
	pcounts_norm <- pcounts
	
	for (j in 1:nObs)
	{
		sample_j <- acounts[,j]/GEOmean
		S[j] <- median(sample_j[which(sample_j != 0)])
		if (sum(acounts[,j]) == 0) {
			acounts_norm[,j] = acounts[,j]
		} else {
			acounts_norm[,j] <- acounts[,j]/S[j]
		}
		if (sum(mcounts[,j]) == 0) {
			mcounts_norm[,j] = mcounts[,j]
		} else {
			mcounts_norm[,j] <- mcounts[,j]/S[j]
		}
		if (sum(pcounts[,j]) == 0) {
			pcounts_norm[,j] = pcounts[,j]
		} else {
			pcounts_norm[,j] <- pcounts[,j]/S[j]
		}
	}

	acounts_norm <- ceiling(acounts_norm)
	mcounts_norm <- ceiling(mcounts_norm)
	pcounts_norm <- ceiling(pcounts_norm)
	
	res = list(acounts_norm, mcounts_norm, pcounts_norm)
	return(res)
}


# fits counts to a regular ZINB distribution, with starting point for ML: nu = pred_nu, mu = pred_mu, sigma = pred_sigma
fitZINB1 <- function(counts, pred_nu = max(0.0001,min(1-0.0001,length(counts[counts==0])/length(counts))), pred_mu = max(0.0001,mean(c(counts[counts!=0]))), pred_sigma = 0.5) {

	if (is.na(pred_mu)) { pred_mu = 0.0001 }
	
	pred_nu2 = max(0.0001,min(1-0.0001,length(counts[counts==0])/length(counts)))
	pred_mu2 = max(0.0001,mean(c(counts[counts!=0])))
	pred_sigma2 = 0.5
	if (is.na(pred_mu2)) { pred_mu2 = 0.0001 }
	
	# constraints (mu and sigma must both be > 0, nu must be between 0 and 1)
	A = matrix(rbind(c(1,0,0),c(-1,0,0),c(0,1,0),c(0,0,1)),4,3)
	B = c(-0.000001,1-0.00001,-0.000001,-0.000001)

	# find estimates for all three params nu, mu and sigma
	opt_me <- function(param){
		nu <- param[1]
		mu <- param[2]
		sigma <- param[3]
		logLH <- sum(dZINBI(counts, mu = mu, sigma = sigma, nu = nu, log = TRUE))
		logLH
	}
	
	res1 <- try(maxLik(logLik = opt_me, start = c(nu = pred_nu, mu = pred_mu, sigma = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res2 <- try(maxLik(logLik = opt_me, start = c(nu = 0.5, mu = pred_mu, sigma = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res3 <- try(maxLik(logLik = opt_me, start = c(nu = 0.9999, mu = pred_mu, sigma = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)

	if ('try-error' %in% class(res1)) { logLH1=NA } else { logLH1 = res1$maximum }
	if ('try-error' %in% class(res2)) { logLH2=NA } else { logLH2 = res2$maximum }
	if ('try-error' %in% class(res3)) { logLH3=NA } else { logLH3 = res3$maximum }
	ll = c(logLH1,logLH2,logLH3); rr = list(res1, res2, res3)

	if (pred_nu2 != pred_nu | pred_mu2 != pred_mu | pred_sigma2 != pred_sigma) {
		res4 <- try(maxLik(logLik = opt_me, start = c(nu = pred_nu2, mu = pred_mu2, sigma = pred_sigma2), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
		if ('try-error' %in% class(res4)) { logLH4=NA } else { logLH4 = res4$maximum }
		ll = c(ll,logLH4); rr[[4]] = res4
	}
		
	if (length(which.max(ll)) == 0) {
		cat("Optimizing log-likelihood in fitZINB1() failed for the following counts:\n")
		print(unname(counts))
		logLH = NA; sigma_pred = NA; mu_pred = NA; nu_pred = NA
	} else {
		res = rr[[which.max(ll)]]
		logLH = res$maximum; nu_pred = unname(res$estimate[1]); mu_pred = unname(res$estimate[2]); sigma_pred = unname(res$estimate[3])
	}
	
	return(c(logLH,nu_pred,mu_pred,sigma_pred))
}

# returns probability of each value in counts under the convolution of two ZINB distributions
# with params nu, mu and sigma
getProbZinbConv <- function(counts, nu, mu, sigma)	{
	res = ifelse(counts==0,(nu^2) + 2*nu*(1-nu)*dNBI(0, mu=mu, sigma=sigma) + ((1-nu)^2)*dNBI(0, mu = 2*mu, sigma = 0.5*sigma),2*nu*(1-nu)*dNBI(counts,mu=mu, sigma=sigma) + ((1-nu)^2)*dNBI(counts, mu = 2*mu, sigma = 0.5*sigma))
	return(res)
}

# fits a convolution of two identically distributed ZINB() distributions
# with params mu, nu and sigma
fitZINB2 <- function(counts, pred_nu = max(0.0001,min(1-0.0001,sqrt(length(counts[counts==0])/length(counts)))), pred_mu = max(0.0001,mean(c(counts[counts!=0]))), pred_sigma = 0.5) {

	if (is.na(pred_mu)) { pred_mu = 0.0001 }

	pred_nu2 = max(0.0001,min(1-0.0001,sqrt(length(counts[counts==0])/length(counts))))
	pred_mu2 = max(0.0001,mean(c(counts[counts!=0])))
	pred_sigma2 = 0.5
	if (is.na(pred_mu2)) { pred_mu2 = 0.0001 }

	# constraints (mu and sigma must both be > 0, nu must be between 0 and 1)
	A = matrix(rbind(c(1,0,0),c(-1,0,0),c(0,1,0),c(0,0,1)),4,3)
	B = c(-0.000001,1-0.00001,-0.000001,-0.000001)
	
	# optimize maternal distribution
	opt_me <- function(param){
		nu <- param[1]
		mu <- param[2]
		sigma <- param[3]
		
		logLH = sum(log(getProbZinbConv(counts,nu,mu,sigma)))
		if (logLH == -Inf) {
			logLH = -1e+10
		}
		logLH
	}
	
	res1 <- try(maxLik(logLik = opt_me, start = c(nu = pred_nu, mu = pred_mu, sigma = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res2 <- try(maxLik(logLik = opt_me, start = c(nu = 0.5, mu = pred_mu, sigma = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res3 <- try(maxLik(logLik = opt_me, start = c(nu = 0.9999, mu = pred_mu, sigma = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	
	if ('try-error' %in% class(res1)) { logLH1=NA } else { logLH1 = res1$maximum }
	if ('try-error' %in% class(res2)) { logLH2=NA } else { logLH2 = res2$maximum }
	if ('try-error' %in% class(res3)) { logLH3=NA } else { logLH3 = res3$maximum }
	ll = c(logLH1,logLH2,logLH3); rr = list(res1, res2, res3)
	
	if (pred_nu2 != pred_nu | pred_mu2 != pred_mu | pred_sigma2 != pred_sigma) {
		res4 <- try(maxLik(logLik = opt_me, start = c(nu = pred_nu2, mu = pred_mu2, sigma = pred_sigma2), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
		if ('try-error' %in% class(res4)) { logLH4=NA } else { logLH4 = res4$maximum }
		ll = c(ll,logLH4); rr[[4]] = res4
	}

	if (length(which.max(ll)) == 0) {
		cat("Optimizing log-likelihood in fitZINB2() failed for the following counts:\n")
		print(unname(counts))
		logLH = NA; sigma_pred = NA; mu_pred = NA; nu_pred = NA
	} else {
		res = rr[[which.max(ll)]]
		logLH = res$maximum; nu_pred = unname(res$estimate[1]); mu_pred = unname(res$estimate[2]); sigma_pred = unname(res$estimate[3])
	}
		
	return(c(logLH,nu_pred,mu_pred,sigma_pred))
}

# fits maternal counts to m ~ ZINB(nu,mu,sigma) + ZINB(nu,mu,sigma) and paternal counts
# to p ~ ZINB(nu,mu,sigma) - in other words, forces them to all have the same params
fitZINBboth <- function(mcounts, pcounts, coef_nu = 1, coef_mu_ZINB = 1, pred_nu = max(0.0001,min(1-0.0001,(sqrt(length(mcounts[mcounts==0])/length(mcounts))+length(pcounts[pcounts==0])/length(pcounts))/2)), pred_mu = max(0.0001,mean(c(mcounts[mcounts!=0],pcounts[pcounts!=0]))), pred_sigma = 0.5) {

	if (is.na(pred_mu)) { pred_mu = 0.0001 }

	# constraints (mu and sigma must both be > 0, nu must be between 0 and 1)
	A = matrix(rbind(c(1,0,0),c(-1,0,0),c(0,1,0),c(0,0,1)),4,3)
	B = c(-0.000001,1-0.00001,-0.000001,-0.000001)
	
	opt_me <- function(param){
		nu <- param[1]
		mu <- param[2]
		sigma <- param[3]
		
		logLH_m = sum(log(getProbZinbConv(mcounts,nu,mu,sigma)))
		logLH_p = sum(dZINBI(pcounts, mu = max(0.000001,mu*coef_mu_ZINB), sigma = sigma, nu = max(0.000001,min(1-0.00001,1-((1-nu)*coef_nu))), log = TRUE))
		
		logLH = logLH_m + logLH_p

		if (logLH == -Inf) {
			logLH = -1e+10
		}	
		logLH	
	}
	
	res1 <- try(maxLik(logLik = opt_me, start = c(nu = pred_nu, mu = pred_mu, sigma = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res2 <- try(maxLik(logLik = opt_me, start = c(nu = 0.5, mu = pred_mu, sigma = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res3 <- try(maxLik(logLik = opt_me, start = c(nu = 0.9999, mu = pred_mu, sigma = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	
	if ('try-error' %in% class(res1)) { logLH1=NA } else { logLH1 = res1$maximum }
	if ('try-error' %in% class(res2)) { logLH2=NA } else { logLH2 = res2$maximum }
	if ('try-error' %in% class(res3)) { logLH3=NA } else { logLH3 = res3$maximum }
	ll = c(logLH1,logLH2,logLH3); rr = list(res1, res2, res3)
		
	if (length(which.max(ll)) == 0) {
		cat("Optimizing log-likelihood in fitZINBboth() failed for the following counts:\n")
		print(unname(mcounts))
		print(unname(pcounts))
		logLH = NA; sigma_pred = NA; mu_pred = NA; nu_pred = NA
	} else {
		res = rr[[which.max(ll)]]
		logLH = res$maximum; nu_pred = unname(res$estimate[1]); mu_pred = unname(res$estimate[2]); sigma_pred = unname(res$estimate[3])
	}
	
	return(c(logLH,nu_pred,mu_pred,sigma_pred))
}

# same as fitZINBboth(), but only constrains mu, not nu or sigma (so nu_m and nu_p can be different)
fitZINBboth2 <- function(mcounts, pcounts, coef_mu_ZINB = 1, pred_nu_m = max(0.0001,min(1-0.0001,sqrt(length(mcounts[mcounts==0])/length(mcounts)))), pred_nu_p = max(0.0001,min(1-0.0001,length(pcounts[pcounts==0])/length(pcounts))), pred_mu = max(0.0001,mean(c(mcounts[mcounts!=0],pcounts[pcounts!=0]))), pred_sigma = 0.5) {

	if (is.na(pred_mu)) { pred_mu = 0.0001 }

	# constraints (mu and sigma must both be > 0, nu must be between 0 and 1)
	A = matrix(rbind(c(1,0,0,0,0),c(-1,0,0,0,0),c(0,1,0,0,0),c(0,-1,0,0,0),c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1)),7,5)
	B = c(-0.000001,1-0.00001,-0.000001,1-0.00001,-0.000001,-0.000001,-0.000001)
	
	opt_me <- function(param){
		nu_m <- param[1]
		nu_p <- param[2]
		mu <- param[3]
		sigma_m <- param[4]
		sigma_p <- param[5]
		
		logLH_m = sum(log(getProbZinbConv(mcounts,nu_m,mu,sigma_m)))
		logLH_p = sum(dZINBI(pcounts, mu = max(0.000001,mu*coef_mu_ZINB), sigma = sigma_p, nu = nu_p, log = TRUE))
		
		logLH = logLH_m + logLH_p

		if (logLH == -Inf) {
			logLH = -1e+10
		}	
		logLH	
	}

	res1 <- try(maxLik(logLik = opt_me, start = c(nu_m = pred_nu_m, nu_p = pred_nu_p, mu = pred_mu, sigma_m = pred_sigma, sigma_p = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res2 <- try(maxLik(logLik = opt_me, start = c(nu_m = 0.5, nu_p = 0.5, mu = pred_mu, sigma_m = pred_sigma, sigma_p = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res3 <- try(maxLik(logLik = opt_me, start = c(nu_m = 0.9999, nu_p = 0.9999, mu = pred_mu, sigma_m = pred_sigma, sigma_p = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	
	if ('try-error' %in% class(res1)) { logLH1=NA } else { logLH1 = res1$maximum }
	if ('try-error' %in% class(res2)) { logLH2=NA } else { logLH2 = res2$maximum }
	if ('try-error' %in% class(res3)) { logLH3=NA } else { logLH3 = res3$maximum }
	ll = c(logLH1,logLH2,logLH3); rr = list(res1, res2, res3)
	
	if (length(which.max(ll)) == 0) {
		cat("Optimizing log-likelihood in fitZINBboth2() failed for following counts:\n")
		print(unname(mcounts))
		print(unname(pcounts))
		logLH = NA; nu_pred_m = NA; nu_pred_p = NA; mu_pred = NA; sigma_pred_m = NA; sigma_pred_p = NA
	} else {
		res = rr[[which.max(ll)]]
		logLH = res$maximum; nu_pred_m = unname(res$estimate[1]); nu_pred_p = unname(res$estimate[2]); mu_pred = unname(res$estimate[3]); sigma_pred_m = unname(res$estimate[4]); sigma_pred_p = unname(res$estimate[5])
	}
	
	return(c(logLH,nu_pred_m,nu_pred_p,mu_pred,sigma_pred_m,sigma_pred_p))
}

# same as fitZINBboth(), but only constrains nu, not mu or sigma (so mu_m and mu_p can be different)
fitZINBboth3 <- function(mcounts, pcounts, coef_nu = 1, pred_nu = max(0.0001,min(1-0.0001,(sqrt(length(mcounts[mcounts==0])/length(mcounts))+length(pcounts[pcounts==0])/length(pcounts))/2)), pred_mu_m = max(0.0001,mean(c(mcounts[mcounts!=0]))), pred_mu_p = max(0.0001,mean(c(pcounts[pcounts!=0]))), pred_sigma = 0.5) {

	if (is.na(pred_mu_m)) { pred_mu_m = 0.0001 }
	if (is.na(pred_mu_p)) { pred_mu_p = 0.0001 }

	# constraints (mu and sigma must both be > 0, nu must be between 0 and 1)
	A = matrix(rbind(c(1,0,0,0,0),c(-1,0,0,0,0),c(0,1,0,0,0),c(0,0,1,0,0),c(0,0,0,1,0),c(0,0,0,0,1)),6,5)
	B = c(-0.000001,1-0.00001,-0.000001,-0.000001,-0.000001,-0.000001)
	
	opt_me <- function(param){
		nu <- param[1]
		mu_m <- param[2]
		mu_p <- param[3]
		sigma_m <- param[4]
		sigma_p <- param[5]
		
		logLH_m = sum(log(getProbZinbConv(mcounts,nu,mu_m,sigma_m)))
		logLH_p = sum(dZINBI(pcounts, mu = mu_p, sigma = sigma_p, nu = max(0.000001,min(1-0.00001,1-((1-nu)*coef_nu))), log = TRUE))
		
		logLH = logLH_m + logLH_p

		if (logLH == -Inf) {
			logLH = -1e+10
		}	
		logLH	
	}

	res1 <- try(maxLik(logLik = opt_me, start = c(nu = pred_nu, mu_m = pred_mu_m, mu_p = max(0.0001,pred_mu_p/(1-pred_nu)), sigma_m = pred_sigma, sigma_p = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res2 <- try(maxLik(logLik = opt_me, start = c(nu = 0.5, mu_m = pred_mu_m, mu_p = pred_mu_p, sigma_m = pred_sigma, sigma_p = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res3 <- try(maxLik(logLik = opt_me, start = c(nu = 0.9999, mu_m = pred_mu_m, mu_p = pred_mu_p, sigma_m = pred_sigma, sigma_p = pred_sigma), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
		
	if ('try-error' %in% class(res1)) { logLH1=NA } else { logLH1 = res1$maximum }
	if ('try-error' %in% class(res2)) { logLH2=NA } else { logLH2 = res2$maximum }
	if ('try-error' %in% class(res3)) { logLH3=NA } else { logLH3 = res3$maximum }
	ll = c(logLH1,logLH2,logLH3); rr = list(res1, res2, res3)
	
	if (length(which.max(ll)) == 0) {
		cat("Optimizing log-likelihood in fitZINBboth3() failed for following counts:\n")
		print(unname(mcounts))
		print(unname(pcounts))
		logLH = NA; nu_pred = NA; mu_pred_m = NA; mu_pred_p = NA; sigma_pred_m = NA; sigma_pred_p = NA
	} else {
		res = rr[[which.max(ll)]]
		logLH = res$maximum; nu_pred = unname(res$estimate[1]); mu_pred_m = unname(res$estimate[2]); mu_pred_p = unname(res$estimate[3]); sigma_pred_m = unname(res$estimate[4]); sigma_pred_p = unname(res$estimate[5])
	}
	
	return(c(logLH,nu_pred,mu_pred_m,mu_pred_p,sigma_pred_m,sigma_pred_p))
}

# fits a single NB() distribution with params mu and sigma (note that since the convolution
# of two NB(mu,sigma) distributions is NB(2*mu, 0.5*sigma), fitting the convolution doesn't 
# require a separate function, just adjust the params accordingly after)
fitNB <- function(counts, pred_mu = max(0.0001,mean(counts)), pred_sigma = max(0.0001,(var(counts) - pred_mu)/(pred_mu^2))) {

	pred_mu2 = max(0.0001,mean(counts))
	pred_sigma2 = max(0.0001,(var(counts) - pred_mu2)/(pred_mu2^2))
	if (is.na(pred_mu)) { pred_mu = 0.0001 }
	if (is.na(pred_mu2)) { pred_mu2 = 0.0001 }

	# constraints (mu and sigma must both be > 0)
	A = matrix(rbind(c(1,0),c(0,1)),2,2)
	B = c(-0.000001,-0.000001)

	# find estimates for all three params nu, mu and sigma
	opt_me <- function(param){
		mu <- param[1]
		sigma <- param[2]
		logLH <- sum(dNBI(counts, mu = mu, sigma = sigma, log = TRUE))
		logLH
	}

	res1 <- try(maxLik(logLik = opt_me, start = c(sigma = pred_sigma, mu = pred_mu), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	res2 <- try(maxLik(logLik = opt_me, start = c(sigma = 100, mu = pred_mu), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
	
	if ('try-error' %in% class(res1)) { logLH1=NA } else { logLH1 = res1$maximum }
	if ('try-error' %in% class(res2)) { logLH2=NA } else { logLH2 = res2$maximum }
	ll = c(logLH1,logLH2); rr = list(res1, res2)
	
	if (pred_mu2 != pred_mu | pred_sigma2 != pred_sigma) {
		res3 <- try(maxLik(logLik = opt_me, start = c(mu = pred_mu2, sigma = pred_sigma2), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
		if ('try-error' %in% class(res3)) { logLH3=NA } else { logLH3 = res3$maximum }
		ll = c(ll,logLH3); rr[[3]] = res3
	}

	if (length(which.max(ll)) == 0) {
		cat("Optimizing log-likelihood in fitNB() failed for following counts:\n")
		print(unname(counts))
		logLH = NA; nu_pred = NA; mu_pred_m = NA; mu_pred_p = NA; sigma_pred_m = NA; sigma_pred_p = NA
	} else {
		res = rr[[which.max(ll)]]
		logLH = res$maximum; mu_pred = unname(res$estimate[1]); sigma_pred = unname(res$estimate[2])
	}
	
	return(c(logLH,mu_pred,sigma_pred))
}

fitNBboth <- function(mcounts, pcounts, coef_mu_NB = 1) {

	pred_mu = max(0.0001,(0.5*mean(mcounts)+mean(pcounts))/1.5)
	if (is.na(pred_mu)) { pred_mu = 0.0001 }
	pred_sigma = max(0.0001,(0.5*((var(mcounts) - pred_mu)/(pred_mu^2))+((var(pcounts) - pred_mu)/(pred_mu^2)))/1.5)

	# constraints (mu and sigma must both be > 0, nu must be between 0 and 1)
	A = matrix(rbind(c(1,0),c(0,1)),2,2)
	B = c(-0.000001,-0.000001)
	
	opt_me <- function(param){
		mu <- param[1]
		sigma <- param[2]
		
		# convolution of two NBs is just NB(2mu,0.5sigma)
		logLH_m = sum(dNBI(mcounts, mu = 2*mu, sigma = 0.5*sigma, log = TRUE))		
		logLH_p = sum(dNBI(pcounts, mu = mu*coef_mu_NB, sigma = sigma, log = TRUE))
		
		logLH = logLH_m + logLH_p

		if (logLH == -Inf) {
			logLH = -1e+10
		}
		
		logLH
	}

	res <- try(maxLik(logLik = opt_me, start = c(sigma = pred_sigma, mu = pred_mu), constraints=list(ineqA=A, ineqB=B)), silent=TRUE)
		
	if ('try-error' %in% class(res)) {	
		cat("Optimizing log-likelihood in fitNBboth() failed for following counts:")
		print(unname(mcounts))
		print(unname(pcounts))
		logLH = NA; sigma_pred = NA; mu_pred = NA
	} else {
		logLH = res$maximum
		mu_pred = res$estimate[1]; names(mu_pred) = NULL
		sigma_pred = res$estimate[2]; names(sigma_pred) = NULL
	}
	
	return(c(logLH,mu_pred,sigma_pred))
}

# uses the functions above to fit counts to appropriate distribution and return logLH and param estimates
fitdist <- function(counts, parent, pred_nu = max(0.0001,min(1-0.0001,length(counts[counts==0])/length(counts))), pred_mu = max(0.0001,mean(c(counts[counts!=0]))), pred_sigma = 0.5, noaic= FALSE, disttouse="both") {

	if (! disttouse %in% c('NB','ZINB','both')) stop("internal error - in calls to fitdist(), disttouse must be either 'NB', 'ZINB', or 'both'")

	if (parent != "mother" & parent != "father") {
		stop("In function fitdist(), second argument must be either 'mother' or 'father'")
	}
		
	# fit both the ZINB and NB distributions
	if (parent == "mother") {
		pred_nu = sqrt(pred_nu)
		zinb_fit = fitZINB2(counts, pred_nu = pred_nu, pred_mu = pred_mu, pred_sigma = pred_sigma)
	} else {
		zinb_fit = fitZINB1(counts, pred_nu = pred_nu, pred_mu = pred_mu, pred_sigma = pred_sigma)	
	}
	nb_fit = fitNB(counts)
	
	# decide which is best? AIC = 2*(# of model params) - 2*(max log likelihood)
	if (noaic == FALSE) {
		aic_zinb = 6 - 2*zinb_fit[1]
		aic_nb = 4 - 2*nb_fit[1]
	} else {
		aic_zinb = -1*zinb_fit[1]
		aic_nb = -1*nb_fit[1]
	}	
	
	if ((aic_zinb < aic_nb || disttouse == "ZINB") && disttouse != "NB") {
		fit = "ZINB"
		nu_estim = zinb_fit[2]
		mu_estim = zinb_fit[3]
		sigma_estim = zinb_fit[4]
		logLH = zinb_fit[1]			
	} else {
		fit = "NB"
		nu_estim = NA
		mu_estim = nb_fit[2]
		sigma_estim = nb_fit[3]	
		logLH = nb_fit[1]			
		
		# adjust estimates if parent is mother (convolution of two NB(mu,sigma) is NB(2mu,0.5sigma)
		if (parent == "mother") {
			mu_estim = 0.5*mu_estim
			sigma_estim = 2*sigma_estim
		}
	}
	res = list(nu_estim, mu_estim, sigma_estim, fit, logLH)
	return(res)
}

testH0 <- function(mcounts, pcounts, results, coef_nu = 1, coef_mu_ZINB = 1, coef_mu_NB = 1) {

	# test inputs: results should be a list of length 37
	if (dim(results)[1] != 1 || dim(results)[2] != 37) {
		stop("Incorrect input to function testH0(), please contact CLP\n")
	}

	# get log-likelihood under H1 (full model)
	logL_p = results$logL_p
	logL_m = results$logL_m
	logL_H1 = logL_m + logL_p
		
	# get log-likelihood under H0, three possibilities:
	# (1) both mcounts and pcounts fit to ZINB	
	if (results$fit_m == "ZINB" & results$fit_p == "ZINB") {

		# get logL under H0 (both nu and mu fixed)
		res1 = fitZINBboth(mcounts,pcounts)
		res2 = fitZINBboth(mcounts,pcounts,coef_nu = coef_nu, coef_mu_ZINB = coef_mu_ZINB)

		fit_H0 = "ZINB"
		logL_H0 = res1[1]; nu_H0 = res1[2]; mu_H0 = res1[3]; sigma_H0 = res1[4]
		logL_H0_mod = res2[1]; nu_H0_mod = res2[2]; mu_H0_mod = res2[3]; sigma_H0_mod = res2[4]

		if (logL_H0 > logL_H1 | logL_H0_mod > logL_H1) {		# often this means the initial fits were off; re-do them using the estimates from the joint
			cat("\nWarning: bad fit in primary tests ZINB/ZINB (logL_H0 =",logL_H0," or logL_H0(mod) =",logL_H0_mod," > logL_H1",logL_H1,"), trying again:\n")
			print(unname(mcounts))
			print(unname(pcounts))
		
			if (res1[1] > res2[1]) {
				nu_H0_m = res1[2]; nu_H0_p = res1[2]; mu_H0_m = res1[3]; mu_H0_p = res1[3]; sigma_H0_m = res1[4]; sigma_H0_p = res1[4];
			} else {
				nu_H0_m = res2[2]; nu_H0_p = max(0.0001,min(1-0.0001,1-((1-res2[2])*coef_nu))); mu_H0_m = res2[3]; mu_H0_p = res2[3]*coef_mu_ZINB; sigma_H0_m = res2[4]; sigma_H0_p = res2[4];
			}

			if (results$mean_m > 0) {
				res_m = fitZINB2(mcounts,pred_nu = nu_H0_m, pred_mu = mu_H0_m, pred_sigma = sigma_H0_m)
				if (res_m[1] > logL_m) {
					results$nu_m = res_m[2]; results$mu_m = res_m[2]; results$sigma_m = res_m[3]
					logL_m = res_m[1]
				}
			} else {
				results$nu_m = 1; results$mu_m = 0; results$sigma_m = 1; results$fit_m = "ZINB"; logL_m = 0
			}

			if (results$mean_p > 0) {
				res_p = fitZINB1(pcounts,pred_nu = nu_H0_p, pred_mu = mu_H0_p, pred_sigma = sigma_H0_p)
				if (res_p[1] > logL_p) {
					results$nu_p = res_p[2]; results$mu_p = res_p[3]; results$sigma_p = res_p[4]
					logL_p = res_p[1]
				}
			} else {
				results$nu_p = 1; results$mu_p = 0; results$sigma_p = 1; results$fit_p = "ZINB"; logL_p = 0
			}
		
			logL_H1 = logL_m + logL_p
		}
		# if still not fixed, die
		if ((logL_H0 - logL_H1) > 0.0001 | (logL_H0_mod - logL_H1) > 0.0001 ) {
			cat("Warning: could not fix bad fit in NB/NB (logL_H0 =",logL_H0," or logL_H0(mod) =",logL_H0_mod," > logL_H1",logL_H1,").\n")
		}

		pval = 1 - pchisq(2 *(logL_H1 - logL_H0), df = 3)
		pval_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_mod), df = 3)

		# if pval or pval_mod are small, also test H0_2: only nu fixed
		if (pval < 0.05 | pval_mod < 0.05) {
			res_H0_2 = fitZINBboth2(mcounts,pcounts); logL_H0_2 = res_H0_2[1]
			res_H0_3 = fitZINBboth3(mcounts,pcounts); logL_H0_3 = res_H0_3[1]			
			res_H0_2_mod = fitZINBboth2(mcounts,pcounts,coef_mu_ZINB = coef_mu_ZINB); logL_H0_2_mod = res_H0_2_mod[1]		
			res_H0_3_mod = fitZINBboth3(mcounts,pcounts,coef_nu = coef_nu); logL_H0_3_mod = res_H0_3_mod[1]
			ll = c(res_H0_2[1],res_H0_3[1],res_H0_2_mod[1],res_H0_3_mod[1])

			if (max(ll) > logL_H1) {
				# try re-doing fit without aic check (sometimes leads to very close calls)
				cat("\nWarning: bad fit in secondary tests ZINB/ZINB (logL_H0_2 =",res_H0_2[1],",logL_H0_3 =",res_H0_3[1],",logL_H0_2_mod =",res_H0_2_mod[1],", or logL_H0_3_mod =",res_H0_3_mod[1]," > logL_H1",logL_H1,"), trying again:\n")
				print(unname(mcounts))
				print(unname(pcounts))
				
				if (max(ll) == res_H0_2[1]) {
					nu_H0_m = res_H0_2[2]; nu_H0_p = res_H0_2[3]; mu_H0_m = res_H0_2[4]; mu_H0_p = res_H0_2[4]; sigma_H0_m = res_H0_2[5]; sigma_H0_p = res_H0_2[6];
				} else if (max(ll) == res_H0_3[1]) {
					nu_H0_m = res_H0_3[2]; nu_H0_p = res_H0_3[2]; mu_H0_m = res_H0_3[3]; mu_H0_p = res_H0_3[4]; sigma_H0_m = res_H0_3[5]; sigma_H0_p = res_H0_3[6];
				} else if (max(ll) == res_H0_2_mod[1]) {
					nu_H0_m = res_H0_2_mod[2]; nu_H0_p = res_H0_2_mod[3]; mu_H0_m = res_H0_2_mod[4]; mu_H0_p = res_H0_2_mod[4]*coef_mu_ZINB; sigma_H0_m = res_H0_2_mod[5]; sigma_H0_p = res_H0_2_mod[6];
				} else {
					nu_H0_m = res_H0_3_mod[2]; nu_H0_p = max(0.0001,min(1-0.0001,1-((1-res_H0_3_mod[2])*coef_nu))); mu_H0_m = res_H0_3_mod[3]; mu_H0_p = res_H0_3_mod[4]; sigma_H0_m = res_H0_3_mod[5]; sigma_H0_p = res_H0_3_mod[6];
				}
														
				if (results$mean_m > 0) {
					res_m = fitZINB2(mcounts,pred_nu = nu_H0_m, pred_mu = mu_H0_m, pred_sigma = sigma_H0_m)
					if (res_m[1] > logL_m) {
						results$nu_m = res_m[2]; results$mu_m = res_m[2]; results$sigma_m = res_m[3]
						logL_m = res_m[1]
					}
				} else {
					results$nu_m = 1; results$mu_m = 0; results$sigma_m = 1; results$fit_m = "ZINB"; logL_m = 0
				}

				if (results$mean_p > 0) {
					res_p = fitZINB1(pcounts,pred_nu = nu_H0_p, pred_mu = mu_H0_p, pred_sigma = sigma_H0_p)
					if (res_p[1] > logL_p) {
						results$nu_p = res_p[2]; results$mu_p = res_p[3]; results$sigma_p = res_p[4]
						logL_p = res_p[1]
					}
				} else {
					results$nu_p = 1; results$mu_p = 0; results$sigma_p = 1; results$fit_p = "ZINB"; logL_p = 0
				}
		
				logL_H1 = logL_m + logL_p
				
				if ((max(ll) - logL_H1) > 0.0001) {
					# still not better
					cat("Warning: could not fix bad fit in secondary tests ZINB/ZINB (logL_H0_2 =",logL_H0_2,",logL_H0_3 =",logL_H0_3,",logL_H0_2_mod =",logL_H0_2_mod,", or logL_H0_3_mod =",logL_H0_3_mod," > logL_H1",logL_H1,") for the following counts:\n")
				}							
			}
			# obtain pvalues for all tests
			pval = 1 - pchisq(2 *(logL_H1 - logL_H0), df = 3)
			pval_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_mod), df = 3)
			pval_H0_2 = 1 - pchisq(2 *(logL_H1 - logL_H0_2), df = 1)
			pval_H0_3 = 1 - pchisq(2 *(logL_H1 - logL_H0_3), df = 1)
			pval_H0_2_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_2_mod), df = 1)
			pval_H0_3_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_3_mod), df = 1)								
			
		} else {
			logL_H0_2 = NA; pval_H0_2 = NA; logL_H0_3 = NA; pval_H0_3 = NA;
			logL_H0_2_mod = NA; pval_H0_2_mod = NA; logL_H0_3_mod = NA; pval_H0_3_mod = NA;	
		}	
	} else if (results$fit_m == "NB" & results$fit_p == "NB") {
		# (2) both distributions fit to NB

		# both maternal and paternal fit to NB, again straightforward LRT
		res1 = fitNBboth(mcounts,pcounts)
		res2 = fitNBboth(mcounts,pcounts,coef_mu_NB = coef_mu_NB)
		fit_H0 = "NB"; logL_H0 = res1[1]; nu_H0 = NA; mu_H0 = res1[2]; sigma_H0 = res1[3]
		logL_H0_mod = res2[1]; nu_H0_mod = NA; mu_H0_mod = res2[2]; sigma_H0_mod = res2[3]

		if (logL_H0 > logL_H1 | logL_H0_mod > logL_H1) {		# often this means the initial fits were off; re-do them using the estimates from the joint
			cat("\nWarning: bad fit in primary tests NB/NB (logL_H0 =",logL_H0," or logL_H0(mod) =",logL_H0_mod," > logL_H1",logL_H1,"), trying again:\n")
			print(unname(mcounts))
			print(unname(pcounts))

			if (res1[1] > res2[1]) {
				mu_H0_m = 2*res1[2]; mu_H0_p = res1[2]; sigma_H0_m = res1[3]/2; sigma_H0_p = res1[3];
			} else {
				mu_H0_m = 2*res2[2]; mu_H0_p = res2[2]*coef_mu_NB; sigma_H0_m = res2[3]/2; sigma_H0_p = res2[3];
			}

			if (results$mean_m > 0) {
				res_m = fitNB(mcounts, pred_mu = mu_H0_m, pred_sigma = sigma_H0_m)
				if (res_m[1] > logL_m) {
					results$mu_m = res_m[2]; results$sigma_m = res_m[3]
					logL_m = res_m[1]
				}
			} else {
				results$nu_m = 1; results$mu_m = 0; results$sigma_m = 1; results$fit_m = "ZINB"; logL_m = 0
			}

			if (results$mean_p > 0) {
				res_p = fitNB(pcounts, pred_mu = mu_H0_p, pred_sigma = sigma_H0_p)
				if (res_p[1] > logL_p) {
					results$mu_p = res_p[2]; results$sigma_p = res_p[3]
					logL_p = res_p[1]
				}
			} else {
				results$nu_p = 1; results$mu_p = 0; results$sigma_p = 1; results$fit_p = "ZINB"; logL_p = 0
			}
		
			logL_H1 = logL_m + logL_p
		}
		# if still not fixed, die
		if ((logL_H0 - logL_H1) > 0.0001 | (logL_H0_mod - logL_H1) > 0.0001 ) {
			cat("Warning: could not fix bad fit (logL_H0 =",logL_H0," or logL_H0(mod) =",logL_H0_mod," > logL_H1",logL_H1,").\n")
		}

		# get pvalues for all tests
		pval = 1 - pchisq(2 *(logL_H1 - logL_H0), df = 2)
		pval_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_mod), df = 2)

		logL_H0_2 = NA; pval_H0_2 = NA; logL_H0_3 = NA; pval_H0_3 = NA;
		logL_H0_2_mod = NA; pval_H0_2_mod = NA; logL_H0_3_mod = NA; pval_H0_3_mod = NA;	
	} else {
		# (3) maternal counts fit ZINB, paternal fit NB or vice-versa; for H0 try both ZINB and NB and pick best fit in H0, then use matched H1 logL
		res1a = fitZINBboth(mcounts,pcounts); res1b = fitZINBboth(mcounts,pcounts,coef_nu = coef_nu, coef_mu_ZINB = coef_mu_ZINB)
		res2a = fitNBboth(mcounts,pcounts); res2b = fitNBboth(mcounts,pcounts,coef_mu_NB = coef_mu_NB)
		if (max(res1a[1],res1b[1]) > max(res2a[1],res2b[1])) {
			fit_H0 = "ZINB"
			logL_H0 = res1a[1]; nu_H0 = res1a[2]; mu_H0 = res1a[3]; sigma_H0 = res1a[4]
			logL_H0_mod = res1b[1]; nu_H0_mod = res1b[2]; mu_H0_mod = res1b[3]; sigma_H0_mod = res1b[4]

			# re-fit NB component with noaic=TRUE (this fixes some corner cases where logL_H0 is very close to logL_H1,
			# but logL_H1 is underestimated b/c NB fit, while better by AIC criteria, did not have literally higher logL than ZINB)
			if (results$fit_m == "NB") {
				res_m = fitdist(mcounts,"mother", noaic=TRUE)
				if (res_m[[5]] > logL_m) {
					results$nu_m = res_m[[1]]; results$mu_m = res_m[[2]]; results$sigma_m = res_m[[3]]; results$fit_m = res_m[[4]]
					logL_m = res_m[[5]]
				}			
			} else {
				res_p = fitdist(pcounts,"father", noaic=TRUE)
				if (res_p[[5]] > logL_p) {
					results$nu_p = res_p[[1]]; results$mu_p = res_p[[2]]; results$sigma_p = res_p[[3]]; results$fit_p = res_p[[4]]
					logL_p = res_p[[5]]
				}
			}

			logL_H1 = logL_m + logL_p

			if (logL_H0 > logL_H1 | logL_H0_mod > logL_H1) {		# often this means the initial fits were off; re-do them using the estimates from the joint
				cat("\nWarning: bad fit in primary tests ZINB/NB (logL_H0 =",logL_H0," or logL_H0(mod) =",logL_H0_mod," > logL_H1",logL_H1,"), trying again:\n")
				print(unname(mcounts))
				print(unname(pcounts))
			
				if (res1a[1] > res1b[1]) {
					nu_H0_m = res1a[2]; nu_H0_p = res1a[2]; mu_H0_m = res1a[3]; mu_H0_p = res1a[3]; sigma_H0_m = res1a[4]; sigma_H0_p = res1a[4];
				} else {
					nu_H0_m = res1b[2]; nu_H0_p = max(0.0001,min(1-0.0001,1-((1-res1b[2])*coef_nu))); mu_H0_m = res1b[3]; mu_H0_p = res1b[3]*coef_mu_ZINB; sigma_H0_m = res1b[4]; sigma_H0_p = res1b[4];
				}
				
				if (results$mean_m > 0) {
					res_m = fitdist(mcounts,"mother", pred_nu = nu_H0_m, pred_mu = mu_H0_m, pred_sigma = sigma_H0_m, noaic=TRUE)
					if (res_m[[5]] > logL_m) {
						results$nu_m = res_m[[1]]; results$mu_m = res_m[[2]]; results$sigma_m = res_m[[3]]; results$fit_m = res_m[[4]]
						logL_m = res_m[[5]]
					}
				} else {
					results$nu_m = 1; results$mu_m = 0; results$sigma_m = 1; results$fit_m = "ZINB"; logL_m = 0
				}

				if (results$mean_p > 0) {
					res_p = fitdist(pcounts,"father", pred_nu = nu_H0_p, pred_mu = mu_H0_p, pred_sigma = sigma_H0_p, noaic=TRUE)
					if (res_p[[5]] > logL_p) {
						results$nu_p = res_p[[1]]; results$mu_p = res_p[[2]]; results$sigma_p = res_p[[3]]; results$fit_p = res_p[[4]]
						logL_p = res_p[[5]]
					}
				} else {
					results$nu_p = 1; results$mu_p = 0; results$sigma_p = 1; results$fit_p = "ZINB"; logL_p = 0
				}
			
				logL_H1 = logL_m + logL_p
			}
			# if still not fixed, warn user
			if ((logL_H0 - logL_H1) > 0.0001 | (logL_H0_mod - logL_H1) > 0.0001 ) {
				cat("Warning: could not fix bad fit (logL_H0 =",logL_H0," or logL_H0(mod) =",logL_H0_mod," > logL_H1",logL_H1,").\n")
			}

			# since we allow re-fit above, in rare cases model that mcounts and pcounts were fit to change, so recalc dfs of full model
			if (results$fit_m == "ZINB" & results$fit_p == "ZINB") {
				df_full = 6
			} else if ((results$fit_m == "ZINB" & results$fit_p == "NB") | (results$fit_m == "NB" & results$fit_p == "ZINB") ) {
				df_full = 5
			} else {
				df_full = 4		# note this really shouldn't happen
			}

			pval = 1 - pchisq(2 *(logL_H1 - logL_H0), df = df_full-3)	
			pval_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_mod), df = df_full-3)

			# if pval or pval_mod are small, also test H0_2: only nu fixed
			if (pval < 0.05 | pval_mod < 0.05) {
				res_H0_2 = fitZINBboth2(mcounts,pcounts); logL_H0_2	= res_H0_2[1]
				res_H0_3 = fitZINBboth3(mcounts,pcounts); logL_H0_3	= res_H0_3[1]			
				res_H0_2_mod = fitZINBboth2(mcounts,pcounts,coef_mu_ZINB = coef_mu_ZINB); logL_H0_2_mod	= res_H0_2_mod[1]		
				res_H0_3_mod = fitZINBboth3(mcounts,pcounts,coef_nu = coef_nu); logL_H0_3_mod = res_H0_3_mod[1]
				ll = c(res_H0_2[1],res_H0_3[1],res_H0_2_mod[1],res_H0_3_mod[1])

				if (max(ll) > logL_H1) {
					# try re-doing fit without aic check and using parameter values discovered above as starting position for ML algorithm
					cat("\nWarning: bad fit in secondary tests ZINB/NB (logL_H0_2 =",res_H0_2[1],",logL_H0_3 =",res_H0_3[1],",logL_H0_2_mod =",res_H0_2_mod[1],", or logL_H0_3_mod =",res_H0_3_mod[1]," > logL_H1",logL_H1,"), trying again:\n")
					print(unname(mcounts))
					print(unname(pcounts))

					if (max(ll) == res_H0_2[1]) {
						nu_H0_m = res_H0_2[2]; nu_H0_p = res_H0_2[3]; mu_H0_m = res_H0_2[4]; mu_H0_p = res_H0_2[4]; sigma_H0_m = res_H0_2[5]; sigma_H0_p = res_H0_2[6];
					} else if (max(ll) == res_H0_3[1]) {
						nu_H0_m = res_H0_3[2]; nu_H0_p = res_H0_3[2]; mu_H0_m = res_H0_3[3]; mu_H0_p = res_H0_3[4]; sigma_H0_m = res_H0_3[5]; sigma_H0_p = res_H0_3[6];
					} else if (max(ll) == res_H0_2_mod[1]) {
						nu_H0_m = res_H0_2_mod[2]; nu_H0_p = res_H0_2_mod[3]; mu_H0_m = res_H0_2_mod[4]; mu_H0_p = res_H0_2_mod[4]*coef_mu_ZINB; sigma_H0_m = res_H0_2_mod[5]; sigma_H0_p = res_H0_2_mod[6];
					} else {
						nu_H0_m = res_H0_3_mod[2]; nu_H0_p = max(0.0001,min(1-0.0001,1-((1-res_H0_3_mod[2])*coef_nu))); mu_H0_m = res_H0_3_mod[3]; mu_H0_p = res_H0_3_mod[4]; sigma_H0_m = res_H0_3_mod[5]; sigma_H0_p = res_H0_3_mod[6];
					}
															
					if (results$mean_m > 0) {
						res_m = fitdist(mcounts,"mother", pred_nu = nu_H0_m, pred_mu = mu_H0_m, pred_sigma = sigma_H0_m, noaic=TRUE)
						if (res_m[[5]] > logL_m) {
							results$nu_m = res_m[[1]]; results$mu_m = res_m[[2]]; results$sigma_m = res_m[[3]]; results$fit_m = res_m[[4]]
							logL_m = res_m[[5]]
						}
					} else {
						results$nu_m = 1; results$mu_m = 0; results$sigma_m = 1; results$fit_m = "ZINB"; logL_m = 0
					}

					if (results$mean_p > 0) {
						res_p = fitdist(pcounts,"father", pred_nu = nu_H0_p, pred_mu = mu_H0_p, pred_sigma = sigma_H0_p, noaic=TRUE)
						if (res_p[[5]] > logL_p) {
							results$nu_p = res_p[[1]]; results$mu_p = res_p[[2]]; results$sigma_p = res_p[[3]]; results$fit_p = res_p[[4]]
							logL_p = res_p[[5]]
						}
					} else {
						results$nu_p = 1; results$mu_p = 0; results$sigma_p = 1; results$fit_p = "ZINB"; logL_p = 0
					}
			
					logL_H1 = logL_m + logL_p
					
					if ((max(ll) - logL_H1) > 0.0001) {
						# still not better
						cat("Warning: could not fix bad fit in secondary tests ZINB/NB (logL_H0_2 =",logL_H0_2,",logL_H0_3 =",logL_H0_3,",logL_H0_2_mod =",logL_H0_2_mod,", or logL_H0_3_mod =",logL_H0_3_mod," > logL_H1",logL_H1,") for the following counts:\n")
					}
					
					# since we allow re-fit above, in rare cases model that mcounts and pcounts fit to change, so recalc dfs of full model
					if (results$fit_m == "ZINB" & results$fit_p == "ZINB") {
						df_full = 6
					} else if ((results$fit_m == "ZINB" & results$fit_p == "NB") | (results$fit_m == "NB" & results$fit_p == "ZINB") ) {
						df_full = 5
					} else {
						df_full = 4		# note this really shouldn't happen
					}						
				}
				# obtain pvalue for all tests
				pval = 1 - pchisq(2 *(logL_H1 - logL_H0), df = df_full-3)
				pval_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_mod), df = df_full-3)
				pval_H0_2 = 1 - pchisq(2 *(logL_H1 - logL_H0_2), df = 1)
				pval_H0_3 = 1 - pchisq(2 *(logL_H1 - logL_H0_3), df = 1)
				pval_H0_2_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_2_mod), df = 1)
				pval_H0_3_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_3_mod), df = 1)								

			} else {
				logL_H0_2 = NA; pval_H0_2 = NA; logL_H0_3 = NA; pval_H0_3 = NA;
				logL_H0_2_mod = NA; pval_H0_2_mod = NA; logL_H0_3_mod = NA; pval_H0_3_mod = NA;	
			}
		} else {
			fit_H0 = "NB"; logL_H0 = res2a[1]; nu_H0 = NA; mu_H0 = res2a[2]; sigma_H0 = res2a[3]; 
			pval = 1 - pchisq(2 *(logL_H1 - logL_H0), df = 3)	# note here 3 dfs freedom b/c full model has 5 params (2 from NB, equiv to nu==0, 3 from ZINB) and restricted model has 2

			# test instead the 'modified' fit
			logL_H0_mod = res2b[1]; nu_H0_mod = NA; mu_H0_mod = res2b[2]; sigma_H0_mod = res2b[3]
			pval_mod = 1 - pchisq(2 *(logL_H1 - logL_H0_mod), df = 3)

			logL_H0_2 = NA; pval_H0_2 = NA; logL_H0_3 = NA; pval_H0_3 = NA;
			logL_H0_2_mod = NA; pval_H0_2_mod = NA; logL_H0_3_mod = NA; pval_H0_3_mod = NA;	
		}				
	}								
	
	# save results
	results$fit_H0 = fit_H0	
	results$nu_H0 = nu_H0
	results$mu_H0 = mu_H0
	results$sigma_H0 = sigma_H0
	results$logL_H1 = logL_H1
	results$logL_H0 = logL_H0
	results$pval = pval			
	results$nu_H0_mod = nu_H0_mod
	results$mu_H0_mod = mu_H0_mod
	results$sigma_H0_mod = sigma_H0_mod
	results$logL_H0_mod = logL_H0_mod
	results$pval_mod = pval_mod
	results$logL_H0_2 = logL_H0_2
	results$pval_H0_2 = pval_H0_2			
	results$logL_H0_2_mod = logL_H0_2_mod
	results$pval_H0_2_mod = pval_H0_2_mod			
	results$logL_H0_3 = logL_H0_3
	results$pval_H0_3 = pval_H0_3			
	results$logL_H0_3_mod = logL_H0_3_mod
	results$pval_H0_3_mod = pval_H0_3_mod	

	return(results)
}























