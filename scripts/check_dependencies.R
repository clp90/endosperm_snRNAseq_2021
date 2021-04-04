#!/usr/bin/env Rscript

# version 1.0 (04/01/2021)
# -------------------------
# Version history:
# v.1.0: initial build - 04/01/2021		CLP
# -------------------------

# Description:
# Quickly checks for all required R library dependencies for the single cell endosperm
# project. Reports any missing libraries to user.

# List of all libraries:
# 	  Package		Version		Citation
#	- argparse		2.0.1		Trevor L Davis (2019). argparse: Command Line Optional and Positional Argument Parser. R package version 2.0.1.
#	- ggplot2		3.3.2		H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#	- RColorBrewer	1.1.2		Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer
#	- viridis		0.5.1		Simon Garnier (2018). viridis: Default Color Maps from 'matplotlib'. R package version 0.5.1. https://CRAN.R-project.org/package=viridis
#	- pheatmap		1.0.12		Raivo Kolde (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12. https://CRAN.R-project.org/package=pheatmap
#	- dplyr			1.0.2		Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2020). dplyr: A Grammar of Data Manipulation. R package version 1.0.2.
#	- optparse		1.6.6		Trevor L Davis (2020). optparse: Command Line Option Parser. R package version 1.6.6. https://CRAN.R-project.org/package=optparse
#	- DEsingle		1.6.0		Zhun Miao, Ke Deng, Xiaowo Wang, Xuegong Zhang. DEsingle for detecting three types of differential expression in single-cell RNA-seq data. Bioinformatics (2018): bty332.
#	- gplots		3.1.0		Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber, Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2020). gplots: Various R Programming Tools for Plotting Data. R package version 3.1.0. https://CRAN.R-project.org/package=gplots
#	- topGO			2.38.1		Adrian Alexa and Jorg Rahnenfuhrer (2019). topGO: Enrichment Analysis for Gene Ontology. R package version 2.38.1.
#	- biomaRt		2.42.1		Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).
#	- Rgraphviz		2.30.0		Kasper Daniel Hansen, Jeff Gentry, Li Long, Robert Gentleman, Seth  Falcon, Florian Hahne and Deepayan Sarkar (2019). Rgraphviz: Provides  plotting capabilities for R graph objects. R package version 2.30.0.
#	- grid			3.6.3		R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
#	- gmodels		2.18.1		Gregory R. Warnes, Ben Bolker, Thomas Lumley, Randall C Johnson. Contributions from Randall C. Johnson are Copyright SAIC-Frederick, Inc. Funded by the Intramural Research Program, of the NIH, National Cancer Institute and Center for Cancer Research under NCI Contract NO1-CO-12400. (2018). gmodels: Various R Programming Tools for Model Fitting. R package version 2.18.1.
#	- gamlss		5.2.0		Rigby R.A. and Stasinopoulos D.M. (2005). Generalized additive models for location, scale and shape,(with discussion), Appl. Statist., 54, part 3, pp 507-554.
#	- VGAM			1.1.3		Thomas W. Yee (2010). The VGAM Package for Categorical Data Analysis. Journal of Statistical Software, 32(10), 1-34. URL http://www.jstatsoft.org/v32/i10/.
#	- maxLik		1.4.4		Arne Henningsen and Ott Toomet (2011). maxLik: A package for maximum likelihood estimation in R. Computational Statistics 26(3), 443-458. DOI 10.1007/s00180-010-0217-1.
#	- metap			1.4			Michael Dewey (2020). metap: meta-analysis of significance values. R package version 1.4.
#	- gridExtra		2.3			Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid" Graphics. R package version 2.3.
#	- vcd			1.4.8		David Meyer, Achim Zeileis, and Kurt Hornik (2020). vcd: Visualizing Categorical Data. R package version 1.4-8.
#	- scater		1.14.6		McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R.” _Bioinformatics_, *33*, 1179-1186.
#	- fpc			2.2.8		Christian Hennig (2020). fpc: Flexible Procedures for Clustering. R package version 2.2-8. https://CRAN.R-project.org/package=fpc
#	- princurve		2.1.5		Cannoodt R., princurve 2.0: Fit a Principal Curve in Arbitrary Dimension (Jun., 2018).
#	- SC3			1.14.0		Vladimir Yu. Kiselev, et al. (2017): SC3 - consensus clustering of single-cell RNA-Seq data. Nature Methods
#	- reshape2		1.4.4		Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20.
#	- Rtsne			0.15		Jesse H. Krijthe (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation
#	- data.table	1.13.0		Matt Dowle and Arun Srinivasan (2020). data.table: Extension of `data.frame`. R package version 1.13.0.
#	- tidyr			1.1.2		Hadley Wickham (2020). tidyr: Tidy Messy Data. R package version 1.1.2. https://CRAN.R-project.org/package=tidyr
#	- igraph		1.2.6		Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006.
#	- maptools		1.0.2		Roger Bivand and Nicholas Lewin-Koh (2020). maptools: Tools for Handling Spatial Objects. R package version 1.0-2. https://CRAN.R-project.org/package=maptools
#	- spatstat		1.64.1		Adrian Baddeley, Ege Rubak, Rolf Turner (2015). Spatial Point Patterns: Methodology and Applications with R. London: Chapman and Hall/CRC Press, 2015
#	- RANN			2.6.1		Sunil Arya, David Mount, Samuel E. Kemp and Gregory Jefferis (2019). RANN: Fast Nearest Neighbour Search (Wraps ANN Library) Using L2 Metric. R package version 2.6.1.
#	- pscl			1.5.5		Simon Jackman (2020). pscl: Classes and Methods for R Developed in the Political Science Computational Laboratory. United States Studies Centre, University of Sydney. Sydney, New South Wales, Australia. R package version 1.5.5.
#	- MASS			7.3.53		Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0
#	- boot			1.3.25		Angelo Canty and Brian Ripley (2020). boot: Bootstrap R (S-Plus) Functions. R package version 1.3-25.
#	- stats			3.6.3		R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
#	- edgeR			3.28.1		Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#	- plyr			1.8.6		Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.

liblist = c("argparse","ggplot2","RColorBrewer","viridis","pheatmap","dplyr","optparse","DEsingle","gplots","topGO","biomaRt","Rgraphviz","grid","gmodels","gamlss","VGAM","maxLik","metap","gridExtra","vcd","scater","fpc","princurve","SC3","reshape2","Rtsne","data.table","tidyr","princurve","igraph","maptools","spatstat","RANN","pscl","MASS","boot","stats","edgeR","plyr")

if(length(setdiff(liblist, rownames(installed.packages()))) > 0) {
	missinglibs = setdiff(liblist, rownames(installed.packages()))
	missinglibs = paste(as.character(missinglibs), sep="' '", collapse="\n - ")
	stop("the following required R libraries are not installed:\n - ",missinglibs,"\nPlease install them, then re-run this script.")
}

cat("All required R libraries were detected.\n")

