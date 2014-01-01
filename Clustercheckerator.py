#! /usr/bin/env python

import re
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

Usage = """
Clustercheckerator.py - version 1.2
Created by Graham Hatfull and Fred Hatfull
Helps to determination Cluster designation of your newly sequenced mycobacteriophage.
Enter the name of a fasta file: i.e. 
Usage: Clustercheckerator.py Yourfileofphagefastafiles.fasta
Results will be return list of cluster designations
Now try again!
"""

def getblastrecords(InputFile): 
	"""This function takes an input phagename.fasta file, searches against a database
	(mycobacteriophages471 in this version), and returns the blast records
	"""
	OutputFile = InputFile[:-6] + ".xml"
	blastn_cline = NcbiblastnCommandline(query=InputFile, db="mycobacteriophages471", evalue=10, outfmt=5, out=OutputFile)#, gapopen=5, gapextend=2, reward=2, penalty=-3)
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

def clusterlookup(blasthit):
	InFile = open('clustertable.csv', 'r')
	lines = InFile.readlines()
	for line in lines:
		line = line.strip('\n').upper()
		SearchStr="%s,(\w+),(\w+)" % blast_hit_name.upper()
		Result = re.search(SearchStr, line)
		if Result is not None:
			cluster = Result.group(1)
			subcluster = Result.group(2)
			InFile.close()
		else:
			cluster = None
			subcluster = None
	clusterassignment = cluster,subcluster
	return clusterassignment

if len(sys.argv)<2:
	print Usage
	sys.exit(1)

InputFile = sys.argv[1]
blast_records = getblastrecords(InputFile)

for blast_record in blast_records:
	queryname = blast_record.query
	if len(blast_record.alignments) == 0:
		# This was included because the blast output has funky alignments with multiple
		#iterations of the same sequence, but with no blast hits recorded.
		continue 		
	elif len(blast_record.alignments) == 1:
		blast_hit_name = blast_record.alignments[0].hit_def
		single_hit_cluster = clusterlookup(blast_hit_name)
		single_hit_name = blast_record.alignments[0].hit_def
		print "The query phage (%s) is the same as %s which is in Cluster %s, Subcluster %s" \
			% (queryname[20:], single_hit_name[21:], single_hit_cluster[0], single_hit_cluster[1])	
	else:
		second_alignment_data = blast_record.alignments[1]
		alignmenttotal = sumhsps(second_alignment_data)
		querylength = blast_record.query_length
		secondhitname = blast_record.alignments[1].hit_def
		secondhitlength = second_alignment_data.length
		percent_hit_to_query = ((float(alignmenttotal) / int(querylength)) *100)
		blast_hit_name = secondhitname[21:]
		secondhitcluster = clusterlookup(secondhitname)
		tophitname = blast_record.alignments[0].hit_def
		Query_cluster = clusterlookup(tophitname)
		print "%s second top hit: %s; spanmatch: %.2f (Cluster: %s, Subcluster: %s)" \
			% (queryname[20:], secondhitname[21:], percent_hit_to_query, secondhitcluster[0], secondhitcluster[1])
		if percent_hit_to_query < 70:
			print "There is no closely related phage and cluster assignment is as \
designated for the query: i.e. Cluster %s, Subcluster %s" % (Query_cluster[0], Query_cluster[1])
		else:
			print "Subcluster assignment: %s" % blast_hit_name















