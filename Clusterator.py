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
	blastn_cline = NcbiblastnCommandline(query=InputFile, db="mycobacteriophages471", evalue=10, outfmt=5, out=OutputFile)#, gapopen=4, gapextend=1, reward=1, penalty=-1)
	stdout, stderr = blastn_cline()

	result_handle = open(OutputFile)
	blast_records = NCBIXML.parse(result_handle)
	return blast_records

def sumhsps(alignment_data):
	"""This function takes an alignment record, identifies the lengths of all 
	hits (hsps), and sums them.
	"""
	alignment_list = []
	alignment_data.hsps.sort(key = lambda hsp: hsp.query_start)
	initial_hsp_start = alignment_data.hsps[0].query_start
	initial_hsp_end = alignment_data.hsps[0].query_end	
	for hsp in alignment_data.hsps:
		hsp_start = hsp.query_start
		hsp_end = hsp.query_end		
		if hsp_end <= initial_hsp_end:
			continue
		if int(hsp_start) <= int(initial_hsp_end):
			initial_hsp_end = hsp_end
		else:
			align_length = int(initial_hsp_end) - int(initial_hsp_start)
			alignment_list.append(align_length)
			initial_hsp_start = hsp_start
			initial_hsp_end = hsp_end
	align_length = int(initial_hsp_end) - int(initial_hsp_start) # this is OK
	alignment_list.append(align_length)
	return sum(alignment_list)
		
def clusterlookup(blasthit):
	InFile = open('clustertable.csv', 'r')
	lines = InFile.readlines()
	for line in lines:
		line = line.strip('\n').upper()
		SearchStr="^%s,(\w+),(\w+)" % blasthit.upper()
		Result = re.match(SearchStr, line)
		if Result is not None:
			cluster = Result.group(1)
			subcluster = Result.group(2)
			InFile.close()
			break
		else:
			cluster = None
			subcluster = None
	clusterassignment = cluster,subcluster
	return clusterassignment

def get_alignment_dict(blast_record):
	"""This function looks up blast records and places alignment span lengths into a dictionary"""
	alignment_length_dictionary = {}
	for each_alignment in blast_record.alignments:
		alignment_total = sumhsps(each_alignment)		
		alignment_name = each_alignment.hit_def[21:]
		querylength = blast_record.query_length
		percent_hit_to_query = ((float(alignment_total) / int(querylength)) *100)		
		alignment_length_dictionary[alignment_name] = percent_hit_to_query
	sorted_alignments = sorted(alignment_length_dictionary.items(), reverse=True, key=lambda x:x[1])
	return sorted_alignments

if __name__ == '__main__':

	if len(sys.argv)<2:
		print Usage
		sys.exit(1)

	InputFile = sys.argv[1]
	blastdata = getblastrecords(InputFile)
	blast_record = next(blastdata)

	if len(blast_record.alignments) == 0:
		print "Congratulations! Your phage is a new singleton Mycobacteriophage"
	else:
		alignments = get_alignment_dict(blast_record)	
		alignment_name, percent_hit_to_query = alignments[0]		
		top_hit_name = alignment_name
		top_hit_percent = "%.2f" % percent_hit_to_query
		top_cluster_assignment = clusterlookup(alignment_name)
		
		top_alignment_data = blast_record.alignments[0]
		top_hit_length = top_alignment_data.length
		
		if int(top_alignment_data.hsps[0].positives) == int(top_hit_length):
			print "Your genome is 100%% identical to %s and is probably the same genome" % top_hit_name

		else:		
			print "Top hit is to %s (%s, %s) percent span length hit: %s" % (top_hit_name, top_cluster_assignment[0], top_cluster_assignment[1], top_hit_percent)

			if percent_hit_to_query > 50:
				print "Your phage has a span length match greater than 50%% to %s and likely also belongs to Cluster %s" % (top_hit_name, top_cluster_assignment[0])
				if percent_hit_to_query > 66 and (top_cluster_assignment[1] != "NONE"):
					print "And with a span length match greater than 65%% is also in Subcluster %s" % top_cluster_assignment[1]
				
				if percent_hit_to_query <= 65 and (top_cluster_assignment[1] != "NONE"):
					print "But your phage matches %s with a span length match less than 65%% and is a candidate for a new Sublcuster" % top_hit_name
			else:
				print "Your phage does not have a span length greater than 50% and is a candidate for a new Singleton"
				if percent_hit_to_query > 40:
					print "But your phage matches %s with a span length greater than 40%% and thus may be in the same cluster (%s)" % (top_hit_name, top_cluster_assignment[0])
				
