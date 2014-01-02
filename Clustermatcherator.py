#! /usr/bin/env python

import re
import sys
import subprocess

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

DEFAULT_BLAST_PARAMS = {
	"db": "mycobacteriophages471",
	"evalue": 10,
	"outfmt": 5,
	"num_threads": 1
}

def blast_records_fasta(record):
	"""This function takes a fasta file format and converts to blast records.
	"""
	blast_params = DEFAULT_BLAST_PARAMS
	blastn_cline = NcbiblastnCommandline(**blast_params)
	child = subprocess.Popen(str(blastn_cline),
		stdin=subprocess.PIPE,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE,
		shell=(sys.platform!="win32"))
	SeqIO.write(record, child.stdin, "fasta")
	child.stdin.close()
	blast_record = NCBIXML.read(child.stdout)
	return blast_record
	
def getblastrecords(InputFile): 
	"""This function takes an input phagename.fasta file, searches against a database
	(mycobacteriophages471 in this version), and returns the blast records
	"""
	OutputFile = InputFile[:-6] + ".xml"
	
	blast_params = DEFAULT_BLAST_PARAMS
	blast_params["query"] = InputFile
	blast_params["out"] = OutputFile
	blastn_cline = NcbiblastnCommandline(**blast_params)
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
	cluster = None
	subcluster = None
	for line in lines:
		line = line.strip('\n').upper()
		SearchStr="^%s,(\w+),(\w+)" % blasthit.upper()
		Result = re.match(SearchStr, line)
		if Result is not None:
			cluster = Result.group(1)
			subcluster = Result.group(2)
			InFile.close()
			break
	clusterassignment = cluster,subcluster
	return clusterassignment

if __name__ == "__main__":

	if len(sys.argv)<2:
		print Usage
		sys.exit(1)

	InputFile = sys.argv[1]

	for seq_record in SeqIO.parse(InputFile, "fasta"):
		blast_record = blast_records_fasta(seq_record)
		queryname = blast_record.query
	
		if len(blast_record.alignments) == 0:
			#This was included because the blast output has funky alignments with multiple
			#iterations of the same sequence, but with no blast hits recorded.
			continue 		

		tophitname = blast_record.alignments[0].hit_def[21:]
		top_blast_hit_name = blast_record.alignments[0].hit_def[21:]
		Query_cluster = clusterlookup(tophitname)
		Query_cluster_des = Query_cluster[0]
		Query_subcluster_des = Query_cluster[1]

		top_hit_cluster = clusterlookup(top_blast_hit_name)
		top_hit_cluster_des = top_hit_cluster[0]
		top_hit_subcluster_des = top_hit_cluster[1]

		second_blast_hit_name = blast_record.alignments[1].hit_def[21:]
		second_alignment_data = blast_record.alignments[1]
		second_alignment_total = sumhsps(second_alignment_data)
		querylength = blast_record.query_length
		second_hit_name = blast_record.alignments[1].hit_def[21:]
		second_blast_hit_name = second_hit_name[21:]
		secondhitlength = second_alignment_data.length
		percent_second_hit_to_query = ((float(second_alignment_total) / int(querylength)) *100)	
		second_hit_cluster = clusterlookup(second_hit_name)
		
		if len(blast_record.alignments) == 1:
		
			cluster_des = top_hit_cluster_des
			subcluster_des = top_hit_subcluster_des
			match = tophitname
			percent = "100"
			
	#		print "The query phage (%s) is the same as %s which is in Cluster %s, Subcluster %s" \
	#			% (queryname[20:], top_blast_hit_name[21:], cluster_des, subcluster_des)	
		
		elif percent_second_hit_to_query < 50:
			cluster_des = top_hit_cluster[0]
			subcluster_des = top_hit_cluster[1]
			match = second_blast_hit_name
			percent = "100"

#			print "Testing: %s second top hit is to %s but matches less than 70% %.2f" % (queryname[20:], second_hit_name[21:], percent_second_hit_to_query)

#			print "%s second top hit is to %s, but matches less than 70% (spanmatch is only %.2f).\
#Cluster designation is thus the same as the query (and the top hit; %s), i.e. Cluster: %s, Subcluster: %s" \
#			% (queryname[20:], second_hit_name[21:], percent_second_hit_to_query, queryname[20:], Query_cluster_des, Query_subcluster_des)

		else:
			cluster_des = second_hit_cluster[0]
			subcluster_des = second_hit_cluster[1]
			match = second_blast_hit_name
			percent = str(percent_second_hit_to_query)
	
#			print "The second top hit to %s is %s (Cluster: %s, Subcluster: %s) and matches \
#more than 70%. %s is thus assigned to %s and %s" \
#				% (queryname[20:], second_hit_name[21:], second_hit_cluster[0], second_hit_cluster[1], 
#				queryname[20:], second_hit_cluster[0], second_hit_cluster[1])
	
	print ",".join([queryname[20:], match, percent, cluster_des, subcluster_des])







