#!/usr/bin/env python2

''' 
-------------------------
Usage: assign_to_allele.py SNPs.bed reads.sam prefix [options]
e.g. assign_to_allele.py Col_Cvi_SNPs.bed reads1.sam sample1

This script goes through each read pair in the provided SAM file and assigns
it to either Ref, Alt, none or both (conflicting) based on whether or not
the paired read overlaps with known SNPs, provided in the file SNPs.txt.
The SNP.bed file is expected to have the following 4 fields (tab-del):
chromosome
start (0-based)
end (start+1, 1-based position)
ref>alt (e.g. T>A, where T = ref base and A = alt base, no spaces)

Output files are:
prefix_`refname'.sam	(reference reads)
prefix_`altname'.sam	(alternate reads)
prefix_none.sam	(neither, no SNPs overlapping or multiple conflicting SNPs).

v.1.4	08/08/2017
by Colette Picard

Version history:
v.1.0 - initial build

v.1.1 - added XR: and XA: tags to all reads, where XR = number of overlapping 
SNPs with the ref allele, and XA = number of overlapping SNPs with the alt allele

v.1.2 - (06/25/2015 by CLP) added *_snp_report output file that tracks # of reads overlapping a snp
that were assigned to ref or alt. also added XL flag to reads that flags all
the SNPs that they overlapped ("none" if no overlap).

v.1.3 - (03/06/2016 by CLP) added option to analyze bisulfite-treated reads, where
C->T transitions (first read of pair on forward strand) and G->A transitions (first
read of pair is on reverse strand) - note that since some SNPs imitate methylation 
(e.g. C->T looks like an unmethylated cytosine in alt but is actually a SNP), the 
bismark methylation string must also be updated according to results of allele 
assignment (see updateMeString()):
	-> C->T SNPs on forward strand and G->A SNPs on the reverse will cause
	the bismark_methylation_extractor to believe that that position is fully unmethylated
	if the read is assigned to alt
	-> T->C SNPs on the forward strand and A->G SNPs on the reverse will have bismark
	reporting no cytosine there, when alt can have cytosine
	-> "none" reads overlapping a CT or GA SNP should have that position ignored, since
	we don't know which strain it came from
- Also fixed issue where unpaired reads may not have been recognized properly
- Also changed previous XR flag, which stored SNPs with ref allele, to XD to avoid
	repeated flags (bismark outputs an XR flag)

v.1.4 - (08/08/2017) fixed issue where read pairs with names containing /1 and /2 for
first and second read of pair, respectively, were never recognized as pairs due to not
sharing the same read name. Script now drops any part of name after a forward slash ('/').

-------------------------
'''
 
import sys, os, argparse, re

if len(sys.argv) == 1:
	print "-------------------------"
	print "assign_to_allele v1.4		by Colette L. Picard, 08/08/2017"
	print "-------------------------"
	print """Given a set of reads (paired or single), assigns reads to either ref or alt
based on the set of SNPs provided. Reads must be sorted by position, and the
SNP file must be in .bed format with the fourth field = ref>alt (e.g. A>T, where A 
is the reference base and T is the alt base). The first field in the SAM file (NAME)
must be the same for pairs.

Uses the CIGAR string to align SNPs in the mismatch to the read sequence. The CIGAR
string characters M,I,D,N,= and X are recognized. All others are unimplemented and
will cause the script to error.

Output files: 
- prefix_`refname'.sam : all reads that were assigned to "ref"
- prefix_`altname'.sam : all reads that were assigned to "alt"
- prefix_none.sam : all read pairs that could not be assigned, or were conflicted
(different SNPs within read supported different parents of origin)
- prefix_snp_report.bed : outputs the list of SNPs that were provided by user,
with two extra fields counting the # of reads assigned to "ref" that overlapped
that SNP, and the # of reads assigned to "alt" that overlapped that SNP.

Output SAM files have two additional SAM tags: XR = number of overlapping SNPs with
the ref allele, XA = # with the alt allele. For reads assigned to a parent of origin,
there is also the XL flag = list of all SNPs that the read overlapped.

If paired reads are included in the alignment file, and one read in the pair is assigned 
to ref or alt while its mate overlaps no snps, the second read will be assigned to the 
same parent of origin as its mate. This behavior is by default restricted to only proper pairs, 
but can be relaxed to any pair which maps to the same chromosome by using the -l option. 
NOTE: this script does account for multiple alignments, provided the RNEXT field is unique. 
Unmapped reads are skipped.

v.1.3 added support for bisulfite-converted reads. Alignments mapped to the forward
strand must ignore C->T SNPs while alignments on the reverse must ignore G->A. Since 
some SNPs imitate methylation (e.g. C->T looks like an unmethylated cytosine in alt but
is actually a SNP), the bismark methylation string must also be updated according to
results of allele assignment (see updateMeString()). 
"""
	print "\nUsage: assign_to_allele.py SNPs.bed reads.sam prefix [options]"
	print "Options:"
	print "--refname ref = name of the reference strain (to use for file naming) - e.g. Col"
	print "--altname alt = name of the alternate strain (to use for file naming) - e.g. Cvi"
	print "--relaxed = flag that enables all pairs to be treated as paired - default restricted to proper pairs"
	print "--chgflag = flag that enables pairs treated as singletons to have SAM flag changed to singleton"
	print "--bisulfite = flag that enables 'bisulfite mode'; C-T and G-A SNPs are ignored where appropriate"
	print "--allow_sub = flag that enables T and A to sub for C and G respectively in SNPs that are not C>T or G>A"
	print "--updatemestr = flag that enables updating of methylation string to fix SNPs that look like methylation info"
	print "-------------------------"
	sys.exit(1)

# read in arguments
parser = argparse.ArgumentParser()
parser.add_argument('snpfile', help = 'SNP file in .bed format (e.g. chr1	12345	12346	T>C)')
parser.add_argument('infile', help = 'SAM file of alignments, sorted by position (e.g. by samtools sort)')
parser.add_argument('prefix', help = 'Prefix for output files (e.g. embryo1)')
parser.add_argument('--refname', default = "ref", type = str,
	help = 'Name of the reference strain for file naming, default "ref"')
parser.add_argument('--altname', default = "alt", type = str,
	help = 'Name of the alternate strain for file naming, default "alt"')
parser.add_argument('--relaxed', default = False, action="store_true",
	help = 'Allow all pairs mapped to the same chromosome to inherit their mate\'s allele assignment (default is proper pairs only)')
parser.add_argument('--chgflag', default = False, action="store_true",
	help = 'For reads that were treated as singletons for the purposes of this script, change SAM flag to reflect this.')
parser.add_argument('--bisulfite', default = False, action="store_true",
	help = 'Enables "bisulfite mode"; C-T and G-A SNPs are ignored where appropriate')
parser.add_argument('--allow_sub', default = False, action="store_true",
	help = 'In bisulfite mode, enables T and A to sub for C and G respectively in SNPs that are not CT or GA')
parser.add_argument('--updatemestr', default = False, action="store_true",
	help = 'In bisulfite mode, enables updating of methylation string to allow for SNPs that look like methylation info')

args = parser.parse_args()

snp_file = args.snpfile
reads = args.infile
prefix = args.prefix
refname = args.refname
altname = args.altname
relaxed = args.relaxed
chgflag = args.chgflag
bisulfite = args.bisulfite
allow_sub = args.allow_sub
updatemestr = args.updatemestr

print "Running assign_to_allele v1.4		by Colette L. Picard, 08/08/2017"
print "-------------------------"
print "Assigning reads from:",reads,"to",refname,"or",altname,"according to SNP file",snp_file
print "Saving to output files:"
print "\t"+prefix+"_"+refname+".sam"
print "\t"+prefix+"_"+altname+".sam"
print "\t"+prefix+"_none.sam (includes conflicted reads, which are also saved to following file)"
print "\t"+prefix+"_confl.sam"
if relaxed:
	print "Proper pair not required for read to inherit mate's assignment, only map to same chr"
else:
	print "Proper pair required for read to inherit mate's assignment"
if chgflag:
	print "SAM flag of paired reads treated as singletons by this script will be set to unpaired"
else:
	print "SAM flags left unchanged"
if bisulfite:
	print "Using bisulfite mode: C->T SNPs on forward strand and G->A SNPs on reverse will be ignored"
if bisulfite and updatemestr:
	print "In bisulfite mode, updating methylation string to mask tricky SNPs"
print "-------------------------"

#-------------------------------------------------------------

# given a cigar string and its starting position (SAM file POS field), returns
# a list of indices starting with POS,... and equal to the length of the read,
# which correspond to the indices in the reference that the read was aligned to
# Example below from http://genome.sph.umich.edu/wiki/SAM:
# ex1.RefPos:     1  2  3  4  5  6  7     8  9 10 11 12 13 14 15 16 17 18 19
#	  Reference:  C  C  A  T  A  C  T     G  A  A  C  T  G  A  C  T  A  A  C
#     Read:                   A  C  T  A  G  A  A     T  G  G  C  T
# This read has POS: 5 and CIGAR: 3M1I3M1D5M
# Want to expand this to: [5,6,7,-,8,9,10,12,13,14,15,16]	* '-' means skip
# ex2.RefPos:     1  2  3  4  5  6  7     8  9 10 11 12 13 14 15 16 17 18 19
#	  Reference:  C  C  A  T  A  C  T     G  A  A  C  T  G  A  C  T  A  A  C
#     Read:                   A  C  T  A  G  .  .  .  .  .  .  C  T	 A  A  C
# This read has POS: 5 and CIGAR: 3M1I1M6N5M (e.g. spans two exons, where intron is AACTGA)
# Want to expand this to: [5,6,7,-,8,15,16,17,18,19]
def expandCIGAR(cigar,pos):
#COM	print "finding indices for cigar",cigar,"starting at",pos
	indices = []
	for match in re.finditer(r'(\d+[MIDNSHPX=])', cigar):
#COM		print "current indices:",indices
#COM		print "current match:",match.group(1)
		type = match.group(1)[-1]
		len = int(match.group(1)[:-1])
#COM		print "current type:",type,"and pos",pos,"and len",len
		if type == 'M' or type == 'X' or type == '=':				# run of aligned bases of length len
			indices = indices + range(pos,pos+len,1); pos = pos + len
		elif type == 'I':											# read has insertion of length len, ignore this base in seq
			indices = indices + ['-']*len
		elif type == 'D' or type == 'N':							# read has deletion of length len, skip ahead in ref
			pos = pos + len
		else:
			print "Error: only CIGAR strings containing M,I,D,X,= or N recognized"
			print "Cigar string is",cigar
			sys.exit(1)
	return indices
		

def assignToParent(read, bisulfite = False):
	# given a read (and the global list of SNPs, which is accessible without providing as param),
	# assigns the read to ref, alt, none or conflicted (both ref and alt SNPs) based on SNPs in SNP file.
	# read consists of: 0=name, 1=SAM flag, 2=chr, 3=startpos, 5=cigar, 7=start_mate, 8=insert_size, 9=sequence
	# accounts for following in CIGAR string: M (aligned - match or mismatch) N (skipped) D (deletion) I (insertion)
	# Additionally, if bisulfite is set to true, will look for XG field in SAM record, which is where bismark
	# reports whether GA or CT conversions are expected (due to strand), and will ignore SNPs masked by those conversions
#	print "read:",read.strip()

	r = read.strip().split('\t')	# parse read
	chr=r[2].lower()
	if not chr in SNPs:				# no SNPs on this chromosome/scaffold, skip
		tagread = read.strip() + '\tXD:i:0\tXA:i:0\n'
		return ["none",tagread]
		
	# if in bisulfite mode, look for XG field which will indicate whether CT or GA conversions occurred on that strand
	if bisulfite:
		if len(read.strip().split('XG:Z:')) == 1:
			print "Error: when using --bisulfite option, SAM files must be generated by bismark, which adds required SAM field XG. Could not find field XG."
			sys.exit(1)
		if len(read.strip().split('XG:Z:')[1]) < 2:
			print "Error: when using --bisulfite option, bismark XG field must be followed by 'CT' or 'GA', but is here followed by:",read.strip().split('XG:Z:')[1]
		toIgnore = read.strip().split('XG:Z:')[1][0:2]
		if toIgnore != "GA" and toIgnore != "CT":
			print "Error: when using --bisulfite option, bismark XG field must be followed by 'CT' or 'GA', but is here followed by:",read.strip().split('XG:Z:')[1]
	else:
		toIgnore = "none"
	
	# otherwise, count # of overlapping SNPs with ref or alt allele
	ref_snps = 0
	alt_snps = 0
	
	# get CIGAR string and translate using expandCIGAR to list of indices to iterate over
	seq = r[9]; seq_start = int(r[3]); cigar = r[5]
	cigar_iter = expandCIGAR(cigar,seq_start)			# list of positions in reference corresponding to positions in read

	try:
		assert len(cigar_iter) == len(seq)				# list should always be same length as read
	except AssertionError:
		print "Internal Error: expandCIGAR returned index list of size",len(cigar_iter),"when length of read is",len(read)
		print "Error for read:",read
		sys.exit(1)
	
	# find overlapping snps, keep track of them
	snp_list = []
	for i in range(0,len(seq)):				# loop through sequence
		refpos = cigar_iter[i]				# position in the reference according to CIGAR string
		if refpos in SNPs[chr]:
			snp_list.append(chr+"."+str(refpos))
			this_snp = SNPs[chr][refpos][0]+SNPs[chr][refpos][1]
#			print "position:",chr+":"+str(refpos)
#			print toIgnore
#			print "Matching SNP is",this_snp
			if toIgnore == "none" or (toIgnore == "CT" and this_snp != "CT" and this_snp != "TC") or (toIgnore == "GA" and this_snp != "GA" and this_snp != "AG"):			
				if seq[i] == SNPs[chr][refpos][0]:
					ref_snps+=1
				elif seq[i] == SNPs[chr][refpos][1]:
					alt_snps+=1
				else:			# bases do not match, allow conversions to be possible? e.g. maybe an observed T is originally a C?
					if bisulfite and allow_sub and toIgnore == "CT":	# allowing substitutions, on forward strand
						if SNPs[chr][refpos][0] == "C" and seq[i] == "T":
#							print "Allowing observed T to substitute for allele C -> assigned to ref allele"
							ref_snps+=1
						elif SNPs[chr][refpos][1] == "C" and seq[i] == "T":
#							print "Allowing observed T to substitute for allele C -> assigned to alt allele"
							alt_snps+=1
					elif bisulfite and allow_sub and toIgnore == "GA":	# allowing substitutions, on reverse strand
						if SNPs[chr][refpos][0] == "G" and seq[i] == "A":
#							print "Allowing observed G to substitute for allele A -> assigned to ref allele"
							ref_snps+=1
						elif SNPs[chr][refpos][1] == "G" and seq[i] == "A":
#							print "Allowing observed G to substitute for allele A -> assigned to alt allele"
							alt_snps+=1
				
	# done looking for SNPs in read, tag read with resulting # of ref snps and alt snps
	# and return call + tagged read	
	if ref_snps > 0 and alt_snps == 0:
		# add 1 to count of ref reads that were assigned based on these snps
		tagread = read.strip() + '\tXD:i:' + str(ref_snps) + '\tXA:i:'  + str(alt_snps) + '\tXL:Z:'
		if len(snp_list) == 0:
			tagread=tagread+"none"
		else:
			for i in range(0,len(snp_list)):
				if i == 0:
					tagread=tagread+snp_list[i]
				else:
					tagread=tagread+","+snp_list[i]			
		return ["ref",tagread + '\n', snp_list]
	elif ref_snps == 0 and alt_snps > 0:
		# add 1 to count of alt reads that were assigned based on these snps
		tagread = read.strip() + '\tXD:i:' + str(ref_snps) + '\tXA:i:'  + str(alt_snps) + '\tXL:Z:'
		if len(snp_list) == 0:
			tagread=tagread+"none"
		else:
			for i in range(0,len(snp_list)):
				if i == 0:
					tagread=tagread+snp_list[i]
				else:
					tagread=tagread+","+snp_list[i]			
		return ["alt",tagread + '\n', snp_list]
	elif ref_snps > 0 and alt_snps > 0:
		tagread = read.strip() + '\tXD:i:' + str(ref_snps) + '\tXA:i:'  + str(alt_snps) + '\n'
		return ["confl",tagread,[]]
	else:
		tagread = read.strip() + '\tXD:i:' + str(ref_snps) + '\tXA:i:'  + str(alt_snps) + '\n'
		return ["none",tagread,[]]

def editSAMflag(read):
	# given a read from a pair, change SAM flag to make it a singleton instead
	r = read.strip().split('\t')	# parse read
	flag = bin(int(r[1]))[2:]; adjflag = '0'*(max(9-len(flag),0))+flag
	if adjflag[-3] == '0' and adjflag[-5] == '0':
		newflag = '0'
	elif adjflag[-3] == '0' and adjflag[-5] == '1':
		newflag = '16'
	elif adjflag[-3] == '1' and adjflag[-5] == '0':
		newflag = '4'

	newread = r[0]+'\t'+newflag
	for i in range(2,len(r)):
		newread = newread+'\t'+r[i]
		
	return newread+'\n'
	
def updateMeString(call, read):
	# since some SNPs can imitate methylation (e.g. a C->T SNP, while not used to assign
	# parent of origin, will be recognized as an unmethylated cytosine by bismark when it 
	# shouldn't, if the read was assigned to alt, but if it was assigned to ref it should count as unmethylated)
	# This script walks through the sequence of the read and methylation string and updates
	# the methylation string using SNPs and the strain it was assigned to.
	# Example:
	#refseq: AACTGAAATAATGAGTCTCTGCAGCGCGGATT
	#	seq: AACTAAAATAATATATCTCTACAACACAAATT
	# 	mes: ....x.......h.h.....x..x.z.zx...
	#	ref:              A           G
	#	alt:              T           A
	# 	des: ....x.......h.h.....x..x...zx... <- what the methylation string should be (if assigned to alt)

	# Example2:
	#refseq: CACCATAAACAAAAAGAATATGAACCAACAAA
	#	seq: TATTATATATAAAAAGAATTTGAATTAATAAA
	# 	mes: h.hh.....h..............hh..h...
	#	ref:        A           A     
	#	alt:        T           C    
	# 	des: h.hh.....h.........x....hh..h... <- what the methylation string should be (if assigned to alt)
	
	# Read above is assigned to alt b/c of A->T SNP. However, at the AC SNP, a T was introduced
	# in alt instead of a C because it's an unmethylated C. This should be included in the mestr.
	
	# Note that for reads overlapping one of the GA/AG or CT/TC SNPs but not assigned to a parent,
	# we can't know whether it came from ref or alt, so to be safe, those positions are censored
	# for those reads.
	
#	print "Updating methylation string for:"
#	print read.strip()
#	print "Read was assigned to:",call
	r = read.strip().split('\t')	# parse read
	chr=r[2].lower()
	if not chr in SNPs or call == "ref":				# no SNPs on this chromosome/scaffold, or was assigned to ref so no changed to make
		return read
	status = read.strip().split('XG:Z:')[1][0:2]
#	print "Status is",status

	# get CIGAR string and translate using expandCIGAR to list of indices to iterate over
	seq = r[9]; seq_start = int(r[3]); cigar = r[5]
	mestr = list(read.strip().split('XM:Z:')[1].split('\t')[0])
	oldmestr = ''.join(mestr)
#	print "Sequence is:",seq
#	print "Mestr is:",''.join(mestr)
	if len(seq) != len(mestr):
		print "Error: methylation string and sequence not same length in following read:"; print read.strip(); print "Sequence:",seq; print "Methylation str:",mestr
		sys.exit(1)
		
	cigar_iter = expandCIGAR(cigar,seq_start)			# list of positions in reference corresponding to positions in read
	for i in range(0,len(seq)):							# loop through sequence		
		refpos = cigar_iter[i]							# position in the reference according to CIGAR string
		if refpos in SNPs[chr]:
			this_snp = SNPs[chr][refpos][0]+SNPs[chr][refpos][1]
#			print "Found SNP:",this_snp
			if status == "CT" and this_snp == "CT" and call == "alt":
				mestr[i] = "."
			elif status == "GA" and this_snp == "GA" and call == "alt":
				mestr[i] = "."
			elif status == "CT" and this_snp[1] == "C" and call == "alt":
				# this site has a cytosine only in the alt strain, get its context:
				try:
					if seq[i+1] == "G":
						context = "z"
					elif seq[i+2] == "G":
						context = "x"
					else:
						context = "h"
					if seq[i] == "C":
						context = context.upper()		# cytosine at this position is methylated
					mestr[i] = context
				except IndexError:
					# cannot determine context b/c sequence is out of range, proceed
					pass
			elif status == "GA" and this_snp[1] == "G" and call == "alt":
				try:
					# this site has a cytosine only in the alt strain, get its context:
					if seq[i-1] == "C":
						context = "z"
					elif seq[i-2] == "C":
						context = "x"
					else:
						context = "h"
					if seq[i] == "G":
						context = context.upper()		# cytosine at this position is methylated
					mestr[i] = context
				except IndexError:
					# cannot determine context b/c sequence is out of range, proceed
					pass
			elif call == "none" and status == "CT" and this_snp == "CT":
				# can't know whether or not this came from ref or alt, to be safe, drop methylation status at this position
				# (note that you don't need to test for this_snp == "TC" because bismark already outputs that as a position with no cytosine
				mestr[i] = "."
			elif call == "none" and status == "GA" and this_snp == "GA":
				mestr[i] = "."
								
	mestr = ''.join(mestr)
	if oldmestr != mestr:
		oldread = read.split('XM:Z:')
		rep = oldread[1].split('\t')
		rep[0] = mestr
		rep = '\t'.join(rep)
		oldread[1] = rep
		return 'XM:Z:'.join(oldread)
	else:
		return read


#-------------------------------------------------------------

# Input the SNP file as a dict, so SNPs = {chr: {pos : [ref,alt]}}.
# SNPs[chr1][12345] will return a 4-element list with list[0] = ref. allele, list[1] = alt.
# list[3] = # of reads assigned to ref that overlapped this SNP, list[4] = # of reads assigned to alt
global SNPs
SNPs = {}
	
print "Reading in SNP file",snp_file
try:
	f = open(snp_file, 'r') 
except IOError, e:
	print e
	print 'Could not open SNP file',snp_file,'...'
	sys.exit(2)

line = f.readline()
while line:
	ll = line.strip().split('\t')	# split line by tabs
	assert len(ll) == 4				# should have 4 fields (see desc. above of inputs)
	# ll[0] = chr, ll[1] = start, ll[2] = end, ll[3] = ref>alt, ll[3][0] = ref, ll[3][2] = alt
	# should be no duplicate positions, but check just in case and exit if not true
	if not (ll[0].lower() in SNPs):
		SNPs[ll[0].lower()] = {}
	# ll[2] = end is correct (1-based) position, ignore ll[1]
	if not (int(ll[2]) in SNPs[ll[0].lower()]):
		SNPs[ll[0].lower()][int(ll[2])] = [ll[3][0],ll[3][2],0,0]
	else:
		print "Error: there is more than one SNP at position",ll[1],"on",ll[0],"please check your SNP file for errors."
		sys.exit(1)
	line = f.readline()
			
f.close()

print "# of SNPs on each chromosome or scaffold:"
for chr in SNPs:
	print chr+"\t"+str(len(SNPs[chr]))

#-------------------------------------------------------------
# Open output files and initialize counters
try:
	out_ref = open(prefix+"_"+refname+".sam", 'w')
	out_alt = open(prefix+"_"+altname+".sam", 'w')
	out_none = open(prefix+"_none.sam", 'w')
	out_confl = open(prefix+"_confl.sam", 'w')
except IOError, e:
	print e
	print 'Could not create output files.'
	sys.exit(2)
	
try:
	f = open(reads, 'r')
except IOError, e:
	print e
	print 'Could not open input SAM file',reads,'...'
	sys.exit(2)

pairs_ref = pairs_alt = pairs_none = pairs_confl = pairs = 0
single_ref = single_alt = single_none = single_confl = single = 0

#-------------------------------------------------------------
# Go through all reads one by one, check status and assign
print "-------------------------"
line = f.readline()
# skip all lines starting with "@" (standard SAM header), add header to all output files
while line.startswith("@"):
	out_ref.write(line)
	out_alt.write(line)
	out_none.write(line)
	out_confl.write(line)
	line = f.readline()
print "Done writing headers to all output files."
print "Now going through all reads from",reads
	
cache = {}			# store information about already assigned reads here so they can find their mates
					# dict = {name+pos: ["ref" or "alt" or "none",line]}
					# clear cache after every chromosome/scaffold
maxCacheSize = 0	# store max size of dict
				
# get chr, start position of first read
r = line.strip().split('\t')
oldChr = r[2].lower(); oldPos = int(r[3])
	
while line:

	r = line.strip().split('\t')
	chr = r[2].lower(); curpos = int(r[3])
	if chr != oldChr:		
		# new chromosome, no mates remain to be found, output all remaining in cache as singletons
		# and clear the cache for the next chromosome
		print "Done reading reads on",oldChr
		unmatched=0
		for key in cache:
			unmatched+=1
			single+=1
			call = cache[key][0]
			if chgflag:
				record = editSAMflag(cache[key][1])
			else:
				record = cache[key][1]	
				
			if bisulfite and updatemestr:
				record = updateMeString(call, record)
				
			if call == "ref":
				for s in cache[key][2]:			# add 1 to ref count for all SNPs overlapping this read
					ss = s.split(".")
					SNPs[ss[0]][int(ss[1])][2]+=1
				single_ref+=1; out_ref.write(record)
			elif call == "alt":
				for s in cache[key][2]:			# add 1 to alt count for all SNPs overlapping this read
					ss = s.split(".")
					SNPs[ss[0]][int(ss[1])][3]+=1
				single_alt+=1; out_alt.write(record)
			elif call == "none":
				single_none+=1; out_none.write(record)
			elif call == "confl":
				single_confl+=1; out_none.write(record)
				out_confl.write(record)
			else:
				print "Internal Error 1: assignToParent returned illegal assignment",call,"of this read:"
				print line
				sys.exit(1)
		if unmatched > 0:
			print str(unmatched),"reads on",oldChr,"remained unmatched in the cache."
		cache.clear()
	else:
		# still on same chromosome, confirm that positions are increasing (i.e. file is sorted)
		try:
			assert curpos >= oldPos
		except AssertionError:
			print "Error: input .sam file is not sorted. Exiting."
			sys.exit(1)
			
	# now, use SAM flag to determine if the read is mapped
	flag = bin(int(r[1]))[2:]; adjflag = '0'*(max(9-len(flag),0))+flag
	if adjflag[-3] == '0':							# read is mapped
	
		# assign read to parent of origin
		call = assignToParent(line, bisulfite)			# returns tuple, first element = call (ref/alt/none), second = read with tags added
		
		# for singletons, save to proper output file, increment proper counter and done
		if adjflag[-1] == '0' or (relaxed and adjflag[-1] == '1') or (not relaxed and adjflag[-1] == '1' and adjflag[-2] == '0'):
		# only treat reads as paired if they are both mapped (if relaxed) or
		# if they are proper pairs (if not relaxed)
#			print "treating as singleton"
			single+=1
			if chgflag:
				record = editSAMflag(call[1])
			else:
				record = call[1]
				
			if bisulfite and updatemestr:
				record = updateMeString(call[0], record)
				
			# print record to correct output file
			if call[0] == "ref":
				for s in call[2]:			# add 1 to ref count for all SNPs overlapping this read
					ss = s.split(".")
					SNPs[ss[0]][int(ss[1])][2]+=1
				single_ref+=1; out_ref.write(record)
			elif call[0] == "alt":
				for s in call[2]:			# add 1 to alt count for all SNPs overlapping this read
					ss = s.split(".")
					SNPs[ss[0]][int(ss[1])][3]+=1
				single_alt+=1; out_alt.write(record)
			elif call[0] == "none":
				single_none+=1; out_none.write(record)
			elif call[0] == "confl":
				single_confl+=1; out_none.write(record)
				out_confl.write(record)
			else:
				print "Internal Error 2: assignToParent returned illegal assignment",call,"of this read:"
				print call
				print line
				sys.exit(1)

		# if the read has a mate that is mapped, try to pair with mate, key is name+pos (name+matepos for mate)
		else:
#			print "looking for pair of"
#			print line.strip()
			name = r[0].split('/')[0]; pos = r[3]; matepos = r[7]		# mate should have name+matepos+pos in cache
			# look for potential mate
			key = name+';'+matepos+';'+pos
			if key in cache:
#				print "Mate is:",cache[key][1].strip()
				# mate was in the cache, determine where it should be assigned then remove from cache
				pairs+=1
				matecall = cache[key][0]		# assignment of the mate
				
				read1 = call[1]
				read2 = cache[key][1]
					
				if call[0] == "none" and matecall == "none":
					if bisulfite and updatemestr:
						read1 = updateMeString("none", read1)
						read2 = updateMeString("none", read2)
					pairs_none+=1; out_none.write(read1+read2)
				elif (call[0] == "ref" and matecall == "none") or (call[0] == "none" and matecall == "ref") or (call[0] == "ref" and matecall == "ref"):
					all_snps = set(call[2] + cache[key][2])		# get all snps overlapped by either mate in the pair without double counting
					for s in all_snps:			# add 1 to ref count for all SNPs overlapping this read pair
						ss = s.split(".")
						SNPs[ss[0]][int(ss[1])][2]+=1					
					if bisulfite and updatemestr:
						read1 = updateMeString("ref", read1)
						read2 = updateMeString("ref", read2)
					pairs_ref+=1; out_ref.write(read1+read2)
				elif (call[0] == "alt" and matecall == "none") or (call[0] == "none" and matecall == "alt") or (call[0] == "alt" and matecall == "alt"):
					all_snps = set(call[2] + cache[key][2])		# get all snps overlapped by either mate in the pair without double counting
					for s in all_snps:			# add 1 to alt count for all SNPs overlapping this read pair
						ss = s.split(".")
						SNPs[ss[0]][int(ss[1])][3]+=1
					if bisulfite and updatemestr:
						read1 = updateMeString("alt", read1)
						read2 = updateMeString("alt", read2)
					pairs_alt+=1; out_alt.write(read1+read2)
				else:
					if bisulfite and updatemestr:
						read1 = updateMeString("none", read1)
						read2 = updateMeString("none", read2)
					pairs_confl+=1; out_none.write(read1+read2)
					out_confl.write(read2+read1)
				
				# delete the mate from the cache to save space
				del cache[key]
				
			else:
				# mate was not in cache, add this read to cache instead
				cache[name+';'+pos+';'+matepos] = call
				if len(cache) > maxCacheSize:
					maxCacheSize = len(cache)
		
	line = f.readline()
	oldPos = curpos; oldChr = chr
	
# done looping through reads, clear what remains in the cache and finish
print "Done reading reads on",oldChr
unmatched = 0
for key in cache:
	single+=1
	unmatched+=1
	call = cache[key][0]
	if chgflag:
		record = editSAMflag(cache[key][1])
	else:
		record = cache[key][1]		
	#COMprint "outputting read",key,"to",call
	if bisulfite and updatemestr:
		record = updateMeString(call, record)
	if call == "ref":
		for s in cache[key][2]:			# add 1 to ref count for all SNPs overlapping this read
			ss = s.split(".")
			SNPs[ss[0]][int(ss[1])][2]+=1
		single_ref+=1; out_ref.write(record)
	elif call == "alt":
		for s in cache[key][2]:			# add 1 to alt count for all SNPs overlapping this read
			ss = s.split(".")
			SNPs[ss[0]][int(ss[1])][3]+=1
		single_alt+=1; out_alt.write(record)
	elif call == "none":
		single_none+=1; out_none.write(record)
	elif call == "confl":
		single_confl+=1; out_none.write(record)
		out_confl.write(record)
	else:
		print "Internal Error 3: assignToParent returned illegal assignment",call,"of this read:"
		print line
		sys.exit(1)
if unmatched > 0:
	print str(unmatched),"reads on",oldChr,"remained unmatched in the cache."
cache.clear()

# also output list of SNPs with added # ref and # alt info
try:
	out_snps = open(prefix+"_snp_report.bed", 'w')
except IOError, e:
	print e
	print 'Could not create output file',prefix+"_snp_report.bed"
	sys.exit(2)

for chr in SNPs:
	for pos in SNPs[chr]:
		out_snps.write(chr+"\t"+str(pos-1)+"\t"+str(pos)+"\t"+SNPs[chr][pos][0]+">"+SNPs[chr][pos][1]+"\t"+str(SNPs[chr][pos][2])+"\t"+str(SNPs[chr][pos][3])+"\n")
	
out_snps.close()

if pairs > 0:
	ppairs_ref = float(pairs_ref) / pairs * 100
	ppairs_alt = float(pairs_alt) / pairs * 100
	ppairs_confl = float(pairs_confl) / pairs * 100
	ppairs_none = float(pairs_none) / pairs * 100
if single > 0:
	psingle_ref = float(single_ref) / single * 100
	psingle_alt = float(single_alt) / single * 100
	psingle_confl = float(single_confl) / single * 100
	psingle_none = float(single_none) / single * 100
	
print "-------------------------"
print "Done assigning reads. Final tally:"
if relaxed:
	print "Number of read pairs treated as pairs (both mates mapped, on the same chromosome):",pairs
else:
	print "Number of read pairs treated as pairs (proper pairs):", pairs
if pairs > 0:
	print "# pairs assigned to %s: %.0f (%0.2f%%)" % ( refname, pairs_ref, ppairs_ref )
	print "# pairs assigned to %s: %.0f (%0.2f%%)" % ( altname, pairs_alt, ppairs_alt )
	print "# conflicted pairs (saved to %s_none.sam): %.0f (%0.2f%%)" % ( prefix, pairs_confl, ppairs_confl )
	print "# pairs that could not be assigned: %.0f (%0.2f%%)" % ( pairs_none, ppairs_none )
else:
	print "# pairs assigned to %s: 0" % ( refname )
	print "# pairs assigned to %s: 0" % ( altname )
	print "# conflicted pairs (saved to %s_none.sam): 0" % ( prefix )
	print "# pairs that could not be assigned: 0"

print "\nNumber of singleton reads (including pairs treated as singles):", single
if single > 0:
	print "# singletons assigned to %s: %.0f (%0.2f%%)" % ( refname, single_ref, psingle_ref )
	print "# singletons assigned to %s: %.0f (%0.2f%%)" % ( altname, single_alt, psingle_alt )
	print "# conflicted singletons (saved to %s_none.sam): %.0f (%0.2f%%)" % ( prefix, single_confl, psingle_confl )
	print "# singletons that could not be assigned: %.0f (%0.2f%%)" % ( single_none, psingle_none )
else:
	print "# singletons assigned to %s: 0" % ( refname )
	print "# singletons assigned to %s: 0" % ( altname )
	print "# conflicted singletons (saved to %s_none.sam): 0" % ( prefix )
	print "# singletons that could not be assigned: 0"

print ""
print "Max size of cache:",maxCacheSize

out_ref.close()
out_alt.close()
out_none.close()

if single_none == 0 and pairs_none == 0 and pairs_confl == 0 and single_confl == 0:
	os.remove(prefix+"_none.sam")
	print prefix+"_none.sam contained no records and was deleted"
if pairs_alt == 0 and single_alt == 0:
	os.remove(prefix+"_"+altname+".sam")
	print prefix+"_"+altname+".sam contained no records and was deleted"
if pairs_ref == 0 and single_ref == 0:
	os.remove(prefix+"_"+refname+".sam")
	print prefix+"_"+refname+".sam contained no records and was deleted"






