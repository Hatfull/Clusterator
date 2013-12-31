#! /usr/bin/env python

import re
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

Entrez.email = "gfh@pitt.edu"

Usage = """
Clusterator.py - version 1.2
Created by Graham Hatfull and Fred Hatfull
Helps to determination Cluster designation of your newly sequenced mycobacteriophage.
Enter the name of a fasta file: i.e. 
Usage: Clusterator.py yourphagename.fasta
Results will be returned as a file named yourphagename.xml
Now try again!
"""

def getblastrecords(InputFile): 
	"""This function takes an input phagename.fasta file, searches against a database
	(mycobacteriophages471 in this version), and returns the blast records
	"""
	OutputFile = InputFile[:-6] + ".xml"
	blastn_cline = NcbiblastnCommandline(query=InputFile, db="mycobacteriophages471", evalue=10, outfmt=5, out=OutputFile)
	stdout, stderr = blastn_cline()

	result_handle = open(OutputFile)
	blast_records = NCBIXML.parse(result_handle)
	return blast_records

def sumhsps(alignment_data):
	"""This function takes an alignment record, identifies the lengths of all 
	hits (hsps), and sums them.
	"""
	alignmentList = []
	for hsp in alignment_data.hsps:
		align_length = hsp.align_length
		alignmentList.append(align_length)
	return sum(alignmentList)

if len(sys.argv)<2:
	print Usage
	sys.exit(1)

InputFile = sys.argv[1]
blastdata = getblastrecords(InputFile)
blast_record = next(blastdata)

if len(blast_record.alignments) == 0:
	print "Congratulations! Your phage is a new singleton Mycobacteriophage"
else:
	top_alignment_data = blast_record.alignments[0]
	alignmenttotal = sumhsps(top_alignment_data)
	querylength = blast_record.query_length
	tophitname = blast_record.alignments[0].hit_def
	tophitlength = top_alignment_data.length
	percent_hit_to_query = ((float(alignmenttotal) / int(querylength)) *100)
	
	print "Top hit is to Mycobacteriophage %s" % tophitname[21:]
	print "Total length of alignments to top hit: %s bp" % alignmenttotal
	print "Percentage match to query: %.2f %%" % percent_hit_to_query
	print "Percentage match to top hit: %.2f %%" % ((float(alignmenttotal) / int(tophitlength)) *100)

	InFile = open('clustertable.csv', 'r')
	lines = InFile.readlines()
	for line in lines:
		line = line.strip('\n').upper()
		SearchStr="%s,(\w+),(\w+)" % tophitname.upper()
		Result = re.search(SearchStr, line)
		if Result is not None:
			cluster = Result.group(1)
			subcluster = Result.group(2)
			InFile.close()
		else:
			cluster = None
			subcluster = None

	if int(top_alignment_data.hsps[0].positives) == int(tophitlength):
		print "Your genome is 100% identical to top hit and is probably the same genome"
		if  len(blast_record.alignments) >1:
			second_hit = blast_record.alignments[1]
			alignmenttotal = sumhsps(second_hit)
			secondquerylength = blast_record.query_length
			secondhitname = blast_record.alignments[1].hit_def
			secondhitlength = second_hit.length
			second_percent_hit_to_query = ((float(alignmenttotal) / int(secondquerylength)) *100)
			print "The second top hit is %s: " % secondhitname[21:]
			print "Percent span match to second hit is %.2f %%:" % second_percent_hit_to_query
	else:
		if percent_hit_to_query > 50:
			print "Your phage belongs to Cluster %s" % cluster
	
		if (subcluster != "NONE") and (percent_hit_to_query > 70):
			print "Your phage matches the top hit more than 70%% of the query genome length and \
likely belongs to Subcluster %s" % subcluster
		else:
			print "Subcluster designation is uncertain, and your phage may represent a new subcluster"
	
	
