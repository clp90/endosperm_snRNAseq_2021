#!/usr/bin/env python

# Note: much of this was inspired by example code provided here:
# https://python-graph-gallery.com/405-dendrogram-with-heatmap-and-coloured-leaves/

''' 
-------------------------
Usage: plot_SC3_heatmap.py SC3_matrix.txt SC3_clusters.txt outprefix [options]

This script takes the similarity matrix output by SC3 (see single_cell_cluster_SC3.R)
and a set of metadata (cell type, stage, etc.), as well as a list of clusters called by
SC3 and (possibly) some information about imprinting, and plots a nice heatmap. It's
nicer than the default SC3 heatmaps. Allows more color customization. Less dumb default
color choices. A few other quality-of-life improvements (see options).

Sample inputs:
- SC3_matrix.txt
nuc_ID	P3_2C	P3_2D	P8_2C	P3_2H
P3_2C	1		0.417	0.367	0.45
P3_2D	0.417	1		0.483	0.65
P8_2C	0.367	0.483	1		0.6
P3_2H	0.45	0.65	0.6		1

- SC3_clusters.txt
nuc_ID	cluster
P3_1D	5
P3_2A	1
P3_2C	1
P3_2D	1

- metadata.txt (optional)
nuc_ID	peak
P3_4C	3N
P3_4D	3N
P3_4E	3N
P3_5B	3N

- imprinting_info.txt (optional)
nuc_ID	AT5G11560	AT5G26210
P3_4C		
P3_4D	0	
P3_4E		
P3_5B		

* Note that the first column in each file is assumed to correspond to cell/nuclei library IDs,
and that the name of the variable (here "nuc_ID") doesn't matter.

All input files must have the same number of rows, and SC3_matrix.txt must have the same number
of rows as columns.

v.1.0	09/05/2018
by Colette Picard

Version history:
v.1.0 - 09/05/2018

-------------------------
'''
 
import sys, os, argparse, re, math, numpy as np
np.set_printoptions(threshold='nan')

if len(sys.argv) == 1:
	print "-------------------------"
	print "plot_SC3_heatmap v1.0		by Colette L. Picard, 09/06/2018"
	print "-------------------------"
	print """This script takes the similarity matrix output by SC3 (see single_cell_cluster_SC3.R)
and a set of metadata (cell type, stage, etc.), as well as a list of clusters called by
SC3 and (possibly) some information about imprinting, and plots a nice heatmap. It's
nicer than the default SC3 heatmaps. Allows more color customization. Less dumb default
color choices. A few other quality-of-life improvements (see options).
"""
	print "\nUsage: plot_SC3_heatmap.py SC3_matrix.txt SC3_clusters.txt outprefix [options]"
	print "Options:"
	print "none at the moment"
	print "-------------------------"
	sys.exit(1)

# read in arguments
parser = argparse.ArgumentParser()
parser.add_argument('similarity', help = 'SC3 "similarity" values for each library to all other libraries (NxN matrix)')
parser.add_argument('clusters', help = 'Cluster ID assigned to each library by SC3')
parser.add_argument('outprefix', help = 'Prefix for output files (e.g. embryo1)')
parser.add_argument('--metadata', default="", help = 'Include metadata from each column in this file in plot')
parser.add_argument('--imprinting', default="", help = 'Including imprinting (% maternal) for these genes in plot')
parser.add_argument('--nocluscolors', default = False, action="store_true", help = 'Do not include a row color-coding clusters')
parser.add_argument('--nospaces', default = False, action="store_true", help = 'Dont try to add spaces between clusters')
parser.add_argument('--log2', default = False, action="store_true", help = 'Take log2 of values before plotting (adds pseudocount of 1 before log2)')
parser.add_argument('--log10', default = False, action="store_true", help = 'Take log10 of values before plotting (adds pseudocount of 1 before log2)')
parser.add_argument('--maincolorimpr', default = False, action="store_true", help = 'Use primary color scheme for imprinting (useful if displaying % maternal instead of similarity matrix values)')
parser.add_argument('--maincolorexpr', default = False, action="store_true", help = 'Use primary color scheme for log2(fpkm)')
args = parser.parse_args()

similarity = args.similarity
clusters = args.clusters
outprefix = args.outprefix

print ""
print "Running plot_SC3_heatmap v1.0		by Colette L. Picard, 09/06/2018"
print "-------------------------"
print "SC3 similarity matrix:",similarity
print "SC3 cluster assignments:",clusters
print "Prefix for output files:",outprefix
print "-------------------------"
print "Additional parameters:"
if args.metadata == "":
	print "No metadata included"
else:
	print "Metadata to include in plot:",args.metadata
if args.imprinting == "":
	print "No imprinting data included"
else:
	print "Imprinting data to include in plot:",args.imprinting
print "-------------------------"
print ""

#-------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------

colorpalettes = ["Set2", "Paired", "hls"]		# color palettes will be used in this order, if
											# other color schemes are not provided

def addspaces(mat_orig_order, clusters, numins = 0):
	# adds spaces between rows and columns in different clusters, for neatness
	# returns a list of cell/nuclei IDs from the similarity matrix (in order of
	# similarity matrix) with numins "gaps" or "blanks" inserted between clusters,
	# as well as before the entire matrix
	
	# Example:
	# mat_orig_order:
	#	[ "P3_4F" "P3_4E" "P3_4C" "P3_4D" ]
	
	# clusters:
	#	geneID  clus        
	#	P3_4C	2
	#	P3_4D	2
	#	P3_4E	1
	#	P3_4F	1
	
	# Returns (with numins == 2):
	# [ "blank_1" "blank_2" "P3_4F" "P3_4E" "blank_3" "blank_4" "P3_4C" "P3_4D" ]

	# create dict mapping cell/nuclei IDs to cluster #
	cluskeys = clus.index.tolist()
	clusvalues = clus.cluster.tolist()
	clusdict = dict(zip(cluskeys,clusvalues))
	
	if len(mat_orig_order) != len(clusdict.keys()):
		print len(mat_orig_order)
		print len(clusdict.keys())
		print "Error: number of cells/nuclei in similarity matrix not equal to number of cells/nuclei in cluster file"
		sys.exit(1)
	
	# insert numins rows/columns between each cluster (numins is calculated as
	# max(1, ceil(# of cells/nuclei in matrix divided by 50)) unless user-specified
	if numins == 0:
		numins = int(max(1,math.ceil(float(len(mat_orig_order))/50)))
	mat_new_order = []
	nn = 0
	for j in range(0,numins):
		mat_new_order.append("blank_"+str(nn))
		nn+=1
		
	for i in range(0,len(mat_orig_order)):
		if i < len(mat_orig_order)-1:
			if clusdict[mat_orig_order[i]] == clusdict[mat_orig_order[i+1]]:
				mat_new_order.append(mat_orig_order[i])
			else:
				mat_new_order.append(mat_orig_order[i])
				for j in range(0,numins):
					mat_new_order.append("blank_"+str(nn))
					nn+=1
		else:
			mat_new_order.append(curlist[i])

	return mat_new_order
		

#-------------------------------------------------------------
# MAIN
#-------------------------------------------------------------

numcolorsused = 0		# number of times a color palette has been used
import matplotlib.colors
#heatmapdefaultcolorscheme = matplotlib.colors.LinearSegmentedColormap.from_list('test',[(0,"#2E5BB1"),(0.5,"#FBFBCE"),(1,"#C81C0D")])
heatmapdefaultcolorscheme = "coolwarm"
imprcolormap = matplotlib.colors.LinearSegmentedColormap.from_list('test',[(0,"#3258B2"),(0.67,"#F6E993"),(1,"#B11711")])
maincolorexpr = matplotlib.colors.LinearSegmentedColormap.from_list('test',[(0,"#FFFFFF"),(0.2,"#CBE4E9"),(0.4,"#AAD2DA"),(0.6,"#8BBFCA"),(0.8,"#6CADBB"),(1,"#529BAC")])

if args.maincolorimpr == True:
	heatmapdefaultcolorscheme = imprcolormap
	
if args.maincolorexpr == True:
	heatmapdefaultcolorscheme = maincolorexpr

# Open all input files, making sure they can be opened
try:
	SC3_mat = open(similarity, 'r') 
except IOError:
	print 'Could not open SC3 similarity matrix file',similarity,'...'
	sys.exit(2)

try:
	SC3_clus = open(clusters, 'r') 
except IOError:
	print 'Could not open SC3 clusters file',clusters,'...'
	sys.exit(2)

if args.metadata != "":
	try:
		metadata = open(args.metadata, 'r') 
	except IOError:
		print 'Could not open metadata file',args.metadata,'...'
		sys.exit(2)

if args.imprinting != "":
	try:
		imprinting = open(args.imprinting, 'r') 
	except IOError:
		print 'Could not open imprinting file',args.imprinting,'...'
		sys.exit(2)

# Load all required modules
try:
	import warnings
except ImportError, e:
	print 'Could not import module "warnings"; is it installed?'
	sys.exit(3)	

warnings.filterwarnings('ignore')

# note: loading seaborn currently gives the following warning (ignored here because of prev. line):
# ShimWarning: The `IPython.html` package has been deprecated since IPython 4.0. You should import from `notebook` instead. ...
# This is not an issue.
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.pyplot as plt
except ImportError, e:
	print 'Could not import module "matplotlib"; is it installed?'
	sys.exit(3)	
	
try:
	import seaborn as sns
except ImportError, e:
	print 'Could not import module "seaborn"; is it installed?'
	sys.exit(3)	
	
try:
	import pandas as pd
except ImportError, e:
	print 'Could not import module "pandas"; is it installed?'
	sys.exit(3)	
	
# Read inputs

# Similarity matrix
df = pd.read_table(similarity)
col1 = list(df.columns.values)[0]
df = df.set_index(col1)

# Take log2 if requested
if args.log2 == True:
	df = np.log2(df+1)
if args.log10 == True:
	df = np.log10(df+1)

# Clusters
clus = pd.read_table(clusters)
col1 = list(clus.columns.values)[0]
clus = clus.set_index(col1)

# Set up basic plot: matrix + color bar indicating clusters

# assign a color to each cluster using a preset color palette
lut = dict(zip(set(clus.cluster), sns.color_palette(colorpalettes[numcolorsused],len(set(clus.cluster)))))
numcolorsused+=1

# add spaces between each cluster to emphasize
#curlist = df.index.tolist()
curlist = df.columns.values

if args.nospaces == False:
	newlist = addspaces(curlist, clus)
	
	df = df.reindex(newlist)	# add in new rows
	df = df.T					# transpose
	df = df.reindex(newlist)	# add in new rows (matrix is symmetrical)

	clus_spaced = clus.reindex(newlist)
	clus = clus_spaced
else:
	newlist = curlist
	
# get list of colors corresponding to each cell/nucleus in matrix, with spacing
collist = [lut[clus.cluster[i]] if not np.isnan(clus.cluster[i]) else (1,1,1) for i in range(0, len(clus.index))]
liblist = clus.index.tolist()

# make new dataframe with each cell/nuc ID and the resulting color
d = {'idx': liblist, 'cluster': collist}
cluscolors = pd.DataFrame(data=d)
cluscolors = cluscolors.set_index('idx')

if args.imprinting == "" and args.metadata == "":
	# if no other arguments, output plot of just matrix + clusters	
	# make plot, adding legend for clusters
	if args.nocluscolors == True:
		with sns.axes_style("white"):
			g=sns.clustermap(df, col_cluster=False, row_cluster = False, cmap=heatmapdefaultcolorscheme, xticklabels = False, yticklabels = False)
	else:	
		with sns.axes_style("white"):
			g=sns.clustermap(df, col_cluster=False, row_cluster = False, cmap=heatmapdefaultcolorscheme, xticklabels = False, yticklabels = False, col_colors=cluscolors)

		for label in lut.keys():
			g.ax_col_dendrogram.bar(0, 0, color=lut[label],label=label, linewidth=0)

		g.ax_col_dendrogram.legend(loc="center", ncol=4, title="cluster")

	plt.show()

if args.metadata != "":
	# metadata provided; add to plot
	md = pd.read_table(metadata)
	col1 = list(md.columns.values)[0]
	md = md.set_index(col1)
	
	# get a list of all factors
	factorlist = list(md.columns.values)
#	print "Adding",len(factorlist),"factors to plot:"
#	for factor in factorlist:
#		print factor

	# map each factor to a color scheme
	factorcolors = pd.DataFrame()
	colorschemes = {}
	for factor in factorlist:
		curfactor = getattr(md, factor)
	
		numcolortouse = numcolorsused % 3
		lut = dict(zip(set(curfactor), sns.color_palette(colorpalettes[numcolortouse],len(set(curfactor)))))
		numcolorsused+=1
		colorschemes[factor] = lut
		curfactor_spaced = curfactor.reindex(newlist)
		
		collist = [lut[curfactor_spaced[i]] if not pd.isnull(curfactor_spaced[i]) else (1,1,1) for i in range(0, len(curfactor_spaced))]
		
		d = {'idx': newlist, factor: collist}
		colors = pd.DataFrame(data=d)
		colors = colors.set_index('idx')
		
		if factorcolors.empty:
			factorcolors = colors
		else:
			factorcolors = factorcolors.merge(colors, left_index=True, right_index=True, how='left')
				
	if args.imprinting == "":	
		# make plot
		if args.nocluscolors == False:
			factorcolors = factorcolors.merge(cluscolors, left_index=True, right_index=True, how='left')
		
		with sns.axes_style("white"):
			g=sns.clustermap(df, col_cluster=False, row_cluster = False, cmap=heatmapdefaultcolorscheme, xticklabels = False, yticklabels = False, col_colors=factorcolors, figsize=(20,20))

		for ff in colorschemes:
#			print "Adding legend for",ff
			for label in colorschemes[ff].keys():
				g.ax_col_dendrogram.bar(0, 0, color=colorschemes[ff][label],label=label, linewidth=0)
			g.ax_col_dendrogram.legend(loc="center", ncol=4, title=ff)

		plt.show()

if args.imprinting != "":
	# imprinting data provided; add to plot
	impr = pd.read_table(imprinting)
	col1 = list(impr.columns.values)[0]
	impr = impr.set_index(col1)

	# map values to a custom color map	
	imprcolormap = matplotlib.colors.LinearSegmentedColormap.from_list('test',[(0,"#3258B2"),(0.67,"#F6E993"),(1,"#B11711")])
	imprcolors = pd.DataFrame()
	
	for gene in list(impr.columns.values):
		curgene = getattr(impr, gene)
		curgene_spaced = curgene.reindex(newlist)
				
		collist = []
		for i in range(0,len(curgene_spaced)):
			if pd.isnull(curgene_spaced[i]):
				collist.append((1,1,1))
			else:
				collist.append(imprcolormap(curgene_spaced[i])[0:3])
		
		d = {'idx': newlist, gene: collist}
		colors = pd.DataFrame(data=d)
		colors = colors.set_index('idx')
		
		if imprcolors.empty:
			imprcolors = colors
		else:
			imprcolors = imprcolors.merge(colors, left_index=True, right_index=True, how='left')
	
		
	# make plot
	if args.metadata == "":	
		if args.nocluscolors == False:
		
			with sns.axes_style("white"):
				if args.nocluscolors == False:
					g=sns.clustermap(df, col_cluster=False, row_cluster = False, cmap=heatmapdefaultcolorscheme, xticklabels = False, yticklabels = False, row_colors = imprcolors, col_colors=cluscolors, figsize=(20,20))
					for label in lut.keys():
						g.ax_col_dendrogram.bar(0, 0, color=lut[label],label=label, linewidth=0)

					g.ax_col_dendrogram.legend(loc="center", ncol=4, title="cluster")
				else:
					g=sns.clustermap(df, col_cluster=False, row_cluster = False, cmap=heatmapdefaultcolorscheme, xticklabels = False, yticklabels = False, row_colors = imprcolors, figsize=(20,20))				

	else:
		if args.nocluscolors == False:
			factorcolors = factorcolors.merge(cluscolors, left_index=True, right_index=True, how='left')
		
		with sns.axes_style("white"):
			g=sns.clustermap(df, col_cluster=False, row_cluster = False, cmap=heatmapdefaultcolorscheme, xticklabels = False, yticklabels = False, row_colors = imprcolors, col_colors=factorcolors, figsize=(20,20))

		for ff in colorschemes:
#			print "Adding legend for",ff
			for label in colorschemes[ff].keys():
				g.ax_col_dendrogram.bar(0, 0, color=colorschemes[ff][label],label=label, linewidth=0)
			g.ax_col_dendrogram.legend(loc="center", ncol=4, title=ff)

	plt.show()	

plt.savefig(outprefix+"_plot.png")







