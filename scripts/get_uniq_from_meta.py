#!/usr/bin/env python

''' 
-------------------------
Usage: get_uniq_from_meta.py reads.sam metachrom.txt outfile.sam

v.1.0	05/17/2018
by Colette Picard

Version history:
v.1.0 - initial build 05/17/2018
v.1.1 - 11/11/2018
	- fixed a very rare corner case - if two reads mapped to the exact same
	location but on two different chromosomes AND those chromosomes were
	not in the metachrom file, it caused a keyerror. Now checks that both
	chromosomes involved are in the metachrom before proceeding (censored if not)
v.1.2 - 12/17/2018
	- script now also adjusts the NH field (number of alignments) for all reads
	back to 1, so htseq-count or other programs won't see these reads as non-uniquely
	mapping.

-------------------------
'''
 
import sys, os, argparse

if len(sys.argv) == 1:
	print "-------------------------"
	print "get_uniq_from_meta v1.2		by Colette L. Picard, 12/17/2018"
	print "-------------------------"
	print """Simple script to extract unique alignments to the regular genome from a list
of alignments from mapping to a metagenome with tophat. When mapping to the metagenome,
uniquely mapped reads (mapQ == 255) are easily extracted and represent reads that can
be specifically aligned to one or the other of the two genomes in the metagenome. However,
reads that align exactly twice, to the same position in both genomes, represent reads
that would map uniquely to either genome alone, and these are harder to extract.

This script reports a single alignment (or pair of alignments for paired-end reads)
for each read that aligns uniquely to both genomes in the metagenome. In the process,
the SAM flag (field 2) is updated to reflect that the alignment is unique, the chromosome
is changed to the "normal" version instead of the strain-specific metagenome version,
and the mapq value is updated to 255 (mapQ == 3 indicates exactly 2 best alignments) for
compatibility with STAR (which outputs unique alignments with mapQ = 255).

NOTE I: input file must have reads sorted by NAME (samtools sort -n)
NOTE II: metachrom.txt is expected to have a header; it is skipped
"""
	print "\nUsage:  get_uniq_from_meta.py reads.sam metachrom.txt outfile.sam"
	print "-------------------------"
	sys.exit(1)

# read in arguments
parser = argparse.ArgumentParser()
parser.add_argument('reads', help = 'Alignments in SAM format')
parser.add_argument('metachrom', help = 'file listing chromosome names both for regular genome and metagenome')
parser.add_argument('outfile', help = 'name of output file')
args = parser.parse_args()

reads = args.reads
fmetachrom = args.metachrom
outfile = args.outfile

print "Running get_uniq_from_meta v1.2		by Colette L. Picard, 12/17/2018"
print "-------------------------"

#-------------------------------------------------------------

def updateSAMflag(oldflag):
	# updates SAM flag to reflect that read is unique best hit, not multimap
	flag = bin(int(oldflag))[2:]; adjflag = '0'*(max(9-len(flag),0))+flag
	newflag=int('0b0' + adjflag[1:], 2)
		
	return str(newflag)
	
#-------------------------------------------------------------

# Input the metachrom file as a dict, so metachrom = {S1chr: normchr}.
metachrom = {}
otherchrom = {}
fullchrlist = []
	
try:
	f = open(fmetachrom, 'r') 
except IOError, e:
	print e
	print 'Could not open file',fmetachrom,'...'
	sys.exit(2)

line = f.readline()		# skip header
line = f.readline()
while line:
	ll = line.strip().split('\t')	# split line by tabs
	assert len(ll) == 3				# should have 3 fields
	# ll[0] = normchr, ll[1] = S1chr, ll[2] = S2chr
	if not (ll[1] in metachrom):
		metachrom[ll[1]] = ll[0]
	else:
		print "Error: chromosome",ll[0],"appears twice"
		sys.exit(1)
	if not (ll[1] in otherchrom):
		otherchrom[ll[1]] = ll[2]
	else:
		print "Error: chromosome",ll[0],"appears twice"
		sys.exit(1)
	if not (ll[1] in fullchrlist):
		fullchrlist.append(ll[1])
	if not (ll[2] in fullchrlist):
		fullchrlist.append(ll[2])
	line = f.readline()
				
f.close()

#-------------------------------------------------------------
# Open output files and initialize counters
try:
	out = open(outfile, 'w')
except IOError, e:
	print e
	print 'Could not create output file.'
	sys.exit(2)
	
try:
	f = open(reads, 'r')
except IOError, e:
	print e
	print 'Could not open input SAM file',reads,'...'
	sys.exit(2)

total_alignments = 0
total_pairs = 0
total_pairs_censored = 0
total_single = 0
total_single_censored = 0
unique_alignments_single = 0
unique_alignments_paired = 0

#-------------------------------------------------------------
# Go through all reads one by one, check status and assign
line = f.readline()
# skip all lines starting with "@" (standard SAM header), add header to all output files
while line.startswith("@"):
	if line.startswith("@SQ"):
		ll = line.strip().split('\t')
		if ll[1][3:] in metachrom:
			ll[1] = 'SN:' + metachrom[ll[1][3:]]
			out.write('\t'.join(ll)+'\n')
	else:
		out.write(line)		
	line = f.readline()
	
print "Now going through all reads from",reads
	
# get name of first read
r = line.strip().split('\t')
oldName = r[0]
curreads = []
	
while line:
	r = line.strip().split('\t')
	name = r[0]
	if name != oldName:		
		# new read found; go through curreads and output if appropriate
#		print "HAVE SET"
#		print len(curreads)
		
		total_alignments+=1
		if len(curreads) == 1:
			print "Error: a single alignment was found for read",curreads[0].split('\t')[0],", are you sure your reads are sorted by name?"
			sys.exit(1)
			
		# not a paired alignment
		if len(curreads) == 2:
#			print "Set is singletons"
			r1 = curreads[0].strip().split('\t')
			r2 = curreads[1].strip().split('\t')
			
			if r1[2] in fullchrlist and r2[2] in fullchrlist:
				total_single+=1
					
				if r1[2] != r2[2]:
					if r1[2] in otherchrom:
						read1 = r1; read2 = r2
					else:
						read1 = r2; read2 = r1
							
					if read1[3] == read2[3]:
						if read2[2] == otherchrom[read1[2]]:
	#						print "OUTPUT as singleton"
							read1[2] = metachrom[read1[2]]
							read1[4] = '255'
							read1[1] = updateSAMflag(read1[1])
							read1.append("XP:Z:none")
							
							# get index of NH:i: field so we can change it to 1 (must be past 11th field of SAM file)
							index = [idx for idx, s in enumerate(read1[10:]) if 'NH:i:' in s][0] + 10
							read1[index] = "NH:i:1"		
							
							# similarly, fix HI field (hit index; numbers multimappers from 1 to NH)
							index = [idx for idx, s in enumerate(read1[10:]) if 'HI:i:' in s][0] + 10
							read1[index] = "HI:i:1"							

							# write completed alignment to file
							out.write('\t'.join(read1)+'\n')
							unique_alignments_single+=1
	#				else:
	#					print "Not included"
	#					print curreads
			else:
				total_single_censored+=1
				total_alignments-=1

		if len(curreads) == 3:
			print "Error: three alignments were found for read",curreads[0].split('\t')[0],", are you sure your reads are sorted by name?"
			sys.exit(1)

		if len(curreads) == 4:
#			print "Set is pairs"
			r1 = curreads[0].strip().split('\t')
			r2 = curreads[1].strip().split('\t')
			r3 = curreads[2].strip().split('\t')
			r4 = curreads[3].strip().split('\t')
		
			if r1[2] in fullchrlist and r2[2] in fullchrlist and r3[2] in fullchrlist and r4[2] in fullchrlist:
				total_pairs+=1

				chrlist = [x[2] for x in [r1,r2,r3,r4]]
				poslist = [x[3] for x in [r1,r2,r3,r4]]
			
				actual = zip(chrlist,poslist)
			
				chrset = set(chrlist)
				posset = set(poslist)
			
				allcombos = [(x,y) for x in chrset for y in posset]
			
				if len(chrset) == 2 and len(posset) == 2 and set(actual) == set(allcombos):
					# read mapped to same location in both genomes; output one
	#				print "OUTPUT as pair"
					unique_alignments_paired+=1
					if r1[2] in metachrom:
						r1[2] = metachrom[r1[2]]
						r1[4] = '255'
						r1[1] = updateSAMflag(r1[1])
						index = [idx for idx, s in enumerate(r1[10:]) if 'NH:i:' in s][0] + 10
						r1[index] = "NH:i:1"							
						index = [idx for idx, s in enumerate(r1[10:]) if 'HI:i:' in s][0] + 10
						r1[index] = "HI:i:1"							
						r1.append("XP:Z:none")
						out.write('\t'.join(r1)+'\n')
					if r2[2] in metachrom:
						r2[2] = metachrom[r2[2]]
						r2[4] = '255'
						r2[1] = updateSAMflag(r2[1])
						index = [idx for idx, s in enumerate(r2[10:]) if 'NH:i:' in s][0] + 10
						r2[index] = "NH:i:1"							
						index = [idx for idx, s in enumerate(r2[10:]) if 'HI:i:' in s][0] + 10
						r2[index] = "HI:i:1"							
						r2.append("XP:Z:none")
						out.write('\t'.join(r2)+'\n')
					if r3[2] in metachrom:
						r3[2] = metachrom[r3[2]]
						r3[4] = '255'
						r3[1] = updateSAMflag(r3[1])
						index = [idx for idx, s in enumerate(r3[10:]) if 'NH:i:' in s][0] + 10
						r3[index] = "NH:i:1"							
						index = [idx for idx, s in enumerate(r3[10:]) if 'HI:i:' in s][0] + 10
						r3[index] = "HI:i:1"							
						r3.append("XP:Z:none")
						out.write('\t'.join(r3)+'\n')					
					if r4[2] in metachrom:
						r4[2] = metachrom[r4[2]]
						r4[4] = '255'
						r4[1] = updateSAMflag(r4[1])
						index = [idx for idx, s in enumerate(r4[10:]) if 'NH:i:' in s][0] + 10
						r4[index] = "NH:i:1"							
						index = [idx for idx, s in enumerate(r4[10:]) if 'HI:i:' in s][0] + 10
						r4[index] = "HI:i:1"							
						r4.append("XP:Z:none")
						out.write('\t'.join(r4)+'\n')						
	#			else:
	#				print "Not included"
	#				print curreads
			else:
				total_pairs_censored+=1
				total_alignments-=1
					
		oldName = name
		curreads = [line]
	else:
		curreads.append(line)
	
	line = f.readline()



print "-------------------------"
print "Done. Summary of results:"
print total_alignments,"total reads analyzed,",total_pairs,"pairs and",total_single,"singletons"
print "Of these,",unique_alignments_paired,"pairs mapped uniquely to the (non-meta) genome"
print "and",unique_alignments_single,"singletons mapped uniquely to the (non-meta) genome"
print "for a total of",unique_alignments_single+unique_alignments_paired,"unique alignments"
print "\nAdditionally",total_pairs_censored,"pairs and",total_single_censored,"singletons were censored (at least one alignment was on a chromosome not in the provided metachrom)"

out.close()





