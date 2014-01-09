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
Mother_of_all_matcherators.py - version 1.2
Created by Graham Hatfull and Fred Hatfull
Generates as matrix of the span lengths of all phages in a database 
Usage: Mother_of_all_matcherators.py Yourfileofphagefastafiles.fasta > outputfilename.csv
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

#def sumhsps(alignment_data):
#	"""This function takes an alignment record, identifies the lengths of all 
#	hits (hsps), and sums them.
#	"""
#	alignmentList = []
#	for hsp in alignment_data.hsps:
#		align_length = hsp.align_length
#		alignmentList.append(align_length)
#	return sum(alignmentList)

def clusterlookup(blasthit):
	"""This function takes a genome and looks up its cluster assignment"""
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

def get_alignment_dict(blast_record):
	"""This function looks up blast records and places alignment span lengths into a dictionary"""
	alignment_length_dictionary = {}
	for each_alignment in blast_record.alignments:
		alignment_total = sumhsps(each_alignment)		
		alignment_name = each_alignment.hit_def[21:]
		querylength = blast_record.query_length
		percent_hit_to_query = ((float(alignment_total) / int(querylength)) *100)		
		alignment_length_dictionary[alignment_name] = percent_hit_to_query
	return alignment_length_dictionary
	
#	sorted_dictionary = sorted(alignment_length_dictionary.items(), reverse=True, key=lambda x:x[1])
	
if __name__ == "__main__":

	if len(sys.argv)<2:
		print Usage
		sys.exit(1)

	InputFile = sys.argv[1]

	alignment_matrix = {}

	for seq_record in SeqIO.parse(InputFile, "fasta"):
		blast_record = blast_records_fasta(seq_record)
		queryname = blast_record.query[20:]
			
		alignments = get_alignment_dict(blast_record)
		alignment_matrix[queryname] = alignments
	
	sorted_genomes = sorted(alignment_matrix.keys(), key=clusterlookup)
	header = "," + ",".join(sorted_genomes)
	print header
	for genome in sorted_genomes:
		row = []
		row.append(genome)
		for alignment in sorted_genomes:
			value = alignment_matrix[genome].get(alignment, 0)
			row.append("%.2f" % value)
		print ",".join(row)		






#		for alignment_name, percent_hit_to_query in sorted_list_of_alignments:
#			percent = "%.2f" % percent_hit_to_query
		
		
		
		
		
		
		
		
			
		#	print ",".join([alignment_name, percent])



"""
		tophitname = blast_record.alignments[0].hit_def[21:]
		blast_hit_name = blast_record.alignments[0].hit_def[21:]
		Query_cluster = clusterlookup(tophitname)
		Query_cluster_des = Query_cluster[0]
		Query_subcluster_des = Query_cluster[1]

		top_hit_cluster = clusterlookup(blast_hit_name)
		top_hit_cluster_des = top_hit_cluster[0]
		top_hit_subcluster_des = top_hit_cluster[1]

		if len(blast_record.alignments) == 1:
		
			cluster_des = top_hit_cluster_des
			subcluster_des = top_hit_subcluster_des
			match = blast_hit_name
			percent = "100"
			
		else:
			blast_hit_name = blast_record.alignments[1].hit_def[21:]
			second_alignment_data = blast_record.alignments[1]
			second_alignment_total = sumhsps(second_alignment_data)
			querylength = blast_record.query_length
			second_hit_name = blast_record.alignments[1].hit_def[21:]
			second_blast_hit_name = second_hit_name[21:]
			secondhitlength = second_alignment_data.length
			percent_second_hit_to_query = ((float(second_alignment_total) / int(querylength)) *100)	
			second_hit_cluster = clusterlookup(second_hit_name)
		
			if percent_second_hit_to_query < 70:
				cluster_des = top_hit_cluster[0]
				subcluster_des = top_hit_cluster[1]
				match = blast_hit_name
				percent = "%.2f" % percent_second_hit_to_query

			else:
				cluster_des = second_hit_cluster[0]
				subcluster_des = second_hit_cluster[1]
				match = second_hit_name
				percent = "%.2f" % percent_second_hit_to_query
	
		print ",".join([queryname[20:], match, percent, cluster_des, subcluster_des])
"""






