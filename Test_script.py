#! /usr/bin/env python

import re
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

def getblastrecords(InputFile): 
	"""This function takes an input phagename.fasta file, searches against a database
	(mycobacteriophages471 in this version), and returns the blast records
	"""
	OutputFile = InputFile[:-6] + ".xml"
	blastn_cline = NcbiblastnCommandline(query=InputFile, db="mycobacteriophages471", evalue=10, outfmt=5, out=OutputFile, culling_limit=2, gapopen=2, gapextend=2, reward=1, penalty=-1)
	stdout, stderr = blastn_cline()

	result_handle = open(OutputFile)
	blast_records = NCBIXML.parse(result_handle)
	return blast_records

def sumhsps(alignment_data):
	"""This function takes an alignment record, identifies the lengths of all non-overlapping 
	hits (hsps), and sums them.
	"""
	alignment_list = []
	initial_hsp_start = blast_record.alignments[0].hsps[0].query_start # this is OK
	initial_hsp_end = blast_record.alignments[0].hsps[0].query_end # this is OK
	for hsp in alignment_data.hsps:
		hsp_start = hsp.query_start
		hsp_end = hsp.query_end		
		if hsp_start >= initial_hsp_start and hsp_end <= initial_hsp_end:
			continue
		if hsp_end <= initial_hsp_start:
			align_length = int(hsp_end) - int(hsp_start)
			alignment_list.append(align_length)
			continue
		if hsp_start >= initial_hsp_end:
			align_length = int(hsp_end) - int(hsp_start)
			alignment_list.append(align_length)
			continue
		if hsp_start < initial_hsp_start:
			initial_hsp_start = hsp_start
		else:
			initial_hsp_end = hsp_end
	align_length = int(initial_hsp_end) - int(initial_hsp_start)
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

	InputFile = sys.argv[1]
	blastdata = getblastrecords(InputFile)
	blast_record = next(blastdata)

	alignments = get_alignment_dict(blast_record)	

	print alignments







