#!/usr/bin/env python

''' 
-------------------------
Usage: make_metagenome.py [options] SNPfile.bed genome.fa outprefix
See below for more detail.

Version history:
v.1.0	03/12/2016		- initial build
v.1.1	11/13/2018		- added option to provide a GTF file, which is then modified to
						work with the metagenome. See description below for more detail.
						- also added check that all SNP ref and alt alleles are in [atgc]
v.1.2	12/07/2018		- added option to provide SNPs for both genomes vs. reference (e.g
						if both strain A and B are not the reference genome)
						

-------------------------
'''
 
import sys, os, argparse, re, textwrap
from Bio import SeqIO
import itertools

if len(sys.argv) == 1:
	print "-------------------------"
	print "make_metagenome v1.2			by Colette L. Picard, 12/07/2018"
	print "-------------------------"
	print """Outputs a "metagenome" of the provided genome using the provided SNP file.
Currently only supports SNPs, but future versions may add support support e.g. for indels.
Will also output a "metachrom" file mapping each renamed chromosome to the original. Example:
SNPs:
Chr1	4	5	T>G
Chr2	2	3	T>A

genome:
>Chr1
ATTCTGAGGTTA
>Chr2
GATTACA

Output metagenome:
>Chr1A
ATTCTGAGGTTA
>Chr2A
GATTACA
>Chr1B
ATTCGGAGGTTA
>Chr2B
GAATACA

Output metachrom:
chr	strain1	strain2
Chr1	Chr1A	Chr1B
Chr2	Chr2A	Chr2B

Will append "A" to names of chromosomes in genome and "B" to SNP-substituted second genome by default, change
this with the --refchar and --altchar options. --refname is name of reference strain, --altname is
name of alt strain (used in metachrom).

Added in v.1.2: if both A and B are not the reference strain, you can provide SNPs for both A -> ref
and B -> ref by providing the A -> ref SNPs using --ASNPs. By default, if --ASNPs is not provided,
A is assumed to be the reference strain, and the SNP file provided in the main args is A (ref) -> B.

Added in v.1.1: can now also provide a GTF file, which will be modified so that it can be used
with the newly-created metagenome. GTFs store gene annotations (exons, introns, etc.) and since 
this script only uses SNP substitutions, no annotations are changed so all the script needs to do 
is append two copies of the GTF together; one with chromosomes from --refname and one with chromosomes
from --altname. For the example above, the following input GTF file:

Input GTF:
Chr1	SRC	gene	2	8	.	+	.	gene_id "gene_1"; transcript_id "gene_1"
Chr1	SRC	exon	3	6	.	+	.	gene_id "gene_1"; transcript_id "gene_1"
Chr2	SRC	gene	1	4	.	-	.	gene_id "gene_2"; transcript_id "gene_2"

is converted to:
Chr1A	SRC	gene	2	8	.	+	.	gene_id "gene_1"; transcript_id "gene_1"
Chr1A	SRC	exon	3	6	.	+	.	gene_id "gene_1"; transcript_id "gene_1"
Chr2A	SRC	gene	1	4	.	-	.	gene_id "gene_2"; transcript_id "gene_2"
Chr1B	SRC	gene	2	8	.	+	.	gene_id "gene_1"; transcript_id "gene_1"
Chr1B	SRC	exon	3	6	.	+	.	gene_id "gene_1"; transcript_id "gene_1"
Chr2B	SRC	gene	1	4	.	-	.	gene_id "gene_2"; transcript_id "gene_2"

Full list of inputs:
snpfile = list of SNPs in BED format (chr start end SNP, where SNP = e.g. C>T if ref == C and alt == T)
genome = full genome in fasta format
outprefix = prefix for output files (one containing the .fa sequence of the metagenome and one containing metachrom info)
--refchar = single character to append to chromosome name for chromosomes from ref file (e.g. Chr1 -> Chr1C) [default A]
--altchar = single character to append to chromosome name for chromosomes from alt file (e.g. Chr1 -> Chr1V) [default B]
--refname = name of reference strain (e.g. Col, default strain1)
--altname = name of alternate strain (e.g. Cvi, default strain2)
--GTF = name of a GTF annotation file to modify in order to make it compatible with the generated metagenome
--ASNPs = if strain A is not the reference, provide ref -> A SNPs here (main snpfile will be ref -> B)
"""
	print "\nUsage: make_metagenome.py [options] SNPfile.bed genome.fa outprefix"
	print "-------------------------"
	sys.exit(1)

# Read in arguments from user
#-------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('snpfile', help = 'list of SNPs in BED format (chr start end SNP, where SNP is e.g. C>T for ref == C, alt == T)')
parser.add_argument('genome', help = 'genome sequence in .fa (FASTA) format')
parser.add_argument('outprefix', help = 'output file prefix (output genome will be in .fa format, metachrom in .txt format)')
parser.add_argument('--refchar', default="A", help = 'a single character that will be appended to chromosome names for ref strain')
parser.add_argument('--altchar', default="B", help = 'a single character that will be appended to chromosome names for alt strain')
parser.add_argument('--refname', default="strain1", help = 'name of reference strain (e.g. Col)')
parser.add_argument('--altname', default="strain2", help = 'name of alternate strain (e.g. Cvi)')
parser.add_argument('--GTF', default="", help = 'name of a GTF annotation file to modify in order to make it compatible with the generated metagenome')
parser.add_argument('--ASNPs', default="", help = 'if strain A is not the reference, name of file containing ref -> A SNPs (main snpfile will be ref -> B)')

args = parser.parse_args()
snpfile = args.snpfile
genome = args.genome
outprefix = args.outprefix
refchar = args.refchar; altchar = args.altchar
refname = args.refname; altname = args.altname

print ""
print "Running make_metagenome v1.0		by Colette L. Picard, 03/12/2016"
print "-------------------------"
print "Creating metagenome from genome file:",genome
print "Using SNPs between",refname,"and",altname,"in:",snpfile
print "Outputting resulting metagenome to:",outprefix+".fa"
print "Outputting resulting metachrom file to:",outprefix+"_metachrom.txt"
if args.GTF != "":
	print "Also making modified GTF file compatible with metagenome, starting from:",args.GTF
print "-------------------------"

allowed_alleles = ["A", "T", "G", "C", "a", "t", "g", "c"]

if args.ASNPs != "":
	ASNPs = {}
	totsnps = 0
	
	print "\nReading in ref -> A SNP file",args.ASNPs
	try:
		f = open(args.ASNPs, 'r') 
	except IOError, e:
		print e
		print 'Could not open ref -> A SNP file',args.ASNPs,'...'
		sys.exit(2)

	line = f.readline()
	while line:
		ll = line.strip().split('\t')	# split line by tabs
		assert len(ll) == 4				# should have 4 fields (see desc. above of inputs)
		# ll[0] = chr, ll[1] = start, ll[2] = end, ll[3] = ref>alt, ll[3][0] = ref, ll[3][2] = alt
		# should be no duplicate positions, but check just in case and exit if not true
		if not (ll[0].lower() in ASNPs):
			ASNPs[ll[0].lower()] = {}
		# ll[1] = 0-based start position, ignore ll[2]
		if not (int(ll[1]) in ASNPs[ll[0].lower()]):
			if ll[3][0] in allowed_alleles and ll[3][2] in allowed_alleles:
				ASNPs[ll[0].lower()][int(ll[1])] = [ll[3][0],ll[3][2]]
				totsnps+=1
			else:
				print "Error: line of SNP file has allele not corresponding to A, T, G, or C:"
				print line
				sys.exit(1)
		else:
			print "Error: there is more than one SNP at position",ll[1],"on",ll[0],"please check your SNP file for errors."
			sys.exit(1)
		line = f.readline()
			
	f.close()

	print "Total number of ref -> A SNPs:",totsnps

	
# Read in SNP file
#-------------------------------------------------------------
# Input the SNP file as a dict, so SNPs = {chr: {pos : [ref,alt]}}.
# SNPs[chr1][12345] will return a 2-element list with list[0] = ref. allele, list[1] = alt.
SNPs = {}
totsnps = 0
	
if args.ASNPs != "":
	print "\nReading in ref -> B SNPs file",snpfile
else:
	print "\nReading in A (ref) -> B SNPs file",snpfile

try:
	f = open(snpfile, 'r') 
except IOError, e:
	print e
	if args.ASNPs != "":
		print "Could not open ref -> B SNPs file",snpfile
	else:
		print "Could not open A (ref) -> B SNPs file",snpfile
	sys.exit(2)

line = f.readline()
while line:
	ll = line.strip().split('\t')	# split line by tabs
	assert len(ll) == 4				# should have 4 fields (see desc. above of inputs)
	# ll[0] = chr, ll[1] = start, ll[2] = end, ll[3] = ref>alt, ll[3][0] = ref, ll[3][2] = alt
	# should be no duplicate positions, but check just in case and exit if not true
	if not (ll[0].lower() in SNPs):
		SNPs[ll[0].lower()] = {}
	# ll[1] = 0-based start position, ignore ll[2]
	if not (int(ll[1]) in SNPs[ll[0].lower()]):
		if ll[3][0] in allowed_alleles and ll[3][2] in allowed_alleles:
			SNPs[ll[0].lower()][int(ll[1])] = [ll[3][0],ll[3][2]]
			totsnps+=1
		else:
			print "Error: line of SNP file has allele not corresponding to A, T, G, or C:"
			print line
			sys.exit(1)
	else:
		print "Error: there is more than one SNP at position",ll[1],"on",ll[0],"please check your SNP file for errors."
		sys.exit(1)
	line = f.readline()
			
f.close()

if args.ASNPs != "":
	print "Total number of ref -> B SNPs:",totsnps
else:
	print "Total number of A (ref) -> B SNPs:",totsnps

# Iterate through chromosomes in provided file
#-------------------------------------------------------------
# open input file
print "-------------------------"
print "Loading genome",genome
try:
	g = open(genome, 'r')
	og = open(outprefix+".fa", 'w')
	om = open(outprefix+"_metachrom.txt", 'w')
except IOError, e:
	print e
	print 'Could not open one of the input or output files'
	sys.exit(2)

om.write("chr	"+refname+"	"+altname+"\n")

# read in genome sequence as a generator to save space
#print "Loading chrs in generator..."
seqs = SeqIO.parse(g,'fasta')
chrlist = {}

num_snp_sub = 0
for chr in seqs:
	# change chromosome names and add to metachrom file
	chrname = chr.id.lower(); chrseq = str(chr.seq); chrnametowrite=chr.id
	print "Processing", chrnametowrite+"..."
	
	refchrname = chrnametowrite + refchar
	altchrname = chrnametowrite + altchar
	om.write(chrnametowrite+"	"+refchrname+"	"+altchrname+"\n")
	
	# list of chromosomes in this genome
	if chrname not in chrlist:
		chrlist[chrname] = [refchrname, altchrname]

#	print "Writing ref chromosome to metagenome..."
	if args.ASNPs == "":
		# write reference chromosome sequence to file
		og.write(">"+refchrname+"\n")
		og.write(chrseq+"\n")
	else:
		# write strain A chromosome sequence to file
		chrseq2 = list(chrseq)
		if chrname in ASNPs:
			for pos in ASNPs[chrname]:
				if pos > len(chrseq2):
					print "Error: provided SNP on",chrname,"position",pos,"is beyond the end of the chromosome (chromosome is",len(chrseq2),"bp long)"
					sys.exit(1)
				if ASNPs[chrname][pos][0] != chrseq2[pos]:
					print "Error: provided SNP on",chrname,"position",pos,"has different ref value ("+ASNPs[chrname][pos][0]+") than genome sequence ("+chrseq2[pos]+")"
					sys.exit(1)
				chrseq2[pos] = ASNPs[chrname][pos][1]
				num_snp_sub+=1
	
		chrseq2 = ''.join(chrseq2)
		og.write(">"+refchrname+"\n")
		og.write(chrseq2+"\n")

	# write strain B chromosome sequence to file
#	print "Substituting in SNPs..."
	chrseq = list(chrseq)
	if chrname in SNPs:
		for pos in SNPs[chrname]:
			if pos > len(chrseq):
				print "Error: provided SNP on",chrname,"position",pos,"is beyond the end of the chromosome (chromosome is",len(chrseq),"bp long)"
				sys.exit(1)
			if SNPs[chrname][pos][0] != chrseq[pos]:
				print "Error: provided SNP on",chrname,"position",pos,"has different ref value ("+SNPs[chrname][pos][0]+") than genome sequence ("+chrseq[pos]+")"
				sys.exit(1)
			chrseq[pos] = SNPs[chrname][pos][1]
			num_snp_sub+=1
	
#	print "Writing alt chromosome to metagenome..."
	chrseq = ''.join(chrseq)
	og.write(">"+altchrname+"\n")
	og.write(chrseq+"\n")

print "-------------------------"
print "DONE.",num_snp_sub,"SNP substitutions made to alt genome"
g.close()


# If GTF file provided, output initially as two temp files w/ chromosomes replaced,
# then append strain2 GTF to strain1 GTF
#-------------------------------------------------------------
# open input file
if args.GTF != "":
	print "-------------------------"
	print "Creating \"meta-GTF\" file from",args.GTF
	try:
		g = open(args.GTF, 'r')
		o1 = open(outprefix+"_tmp1.gtf", 'w')
		o2 = open(outprefix+"_tmp2.gtf", 'w')
	except IOError, e:
		print e
		print 'Could not open one of the input or output files'
		sys.exit(2)
	
	line = g.readline()
	while line:
		ll = line.strip().split('\t')	# split line by tabs
		if ll[0].lower() not in chrlist:
			print "Error: chromosome",ll[0].lower(),"found in GTF file but not in genome"
			sys.exit(1)
		
		refchrname = chrlist[ll[0].lower()][0]
		altchrname = chrlist[ll[0].lower()][1]
	
		ll[0] = refchrname	
		o1.write('\t'.join(ll)+'\n')	
		ll[0] = altchrname
		o2.write('\t'.join(ll)+'\n')
	
		line = g.readline()
	
	g.close(); o1.close(); o2.close()

	# append o2 to o1
	with open(outprefix+"_metagtf.gtf", 'w') as outfile:
		for fname in [outprefix+"_tmp1.gtf",outprefix+"_tmp2.gtf"]:
			with open(fname) as infile:
				for line in infile:
					outfile.write(line)

	# delete tmp files
	os.remove(outprefix+"_tmp1.gtf")
	os.remove(outprefix+"_tmp2.gtf")



