#! /usr/bin/env python

import re
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML


def clusterlookup(phagename):
	InFile = open('clustertable599.csv', 'r')
	lines = InFile.readlines()
	for line in lines:
		line = line.strip('\n').upper()
		SearchStr="^%s,(\w+),(\w+)" % item
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


if __name__ == '__main__':

#	if len(sys.argv)<2:
#		print Usage
#		sys.exit(1)

	InFile = open('testresults.csv', 'r')
	lines = InFile.readlines()
	header = lines[0]
	for item in header.split(','):
		item = item.upper()
		cluster_assignment = clusterlookup(item)
		cluster = cluster_assignment[0]
		subcluster = cluster_assignment[1]
		print cluster, subcluster
	

#	print header
#	for line in lines:
#		line = line.strip('\n').upper()
#		SearchStr="(\w).*" 
#		result = re.match(SearchStr, line)
#		print result
#			phage_name = result
#			cluster_assignment = clusterlookup(phage_name)
#			cluster = cluster_assignment[0]
#			subcluster = cluster_assignment[1]
#			replacement_string = result, cluster, subcluster
#			print replacement_string

#		replacement = re.sub"(^(\w.)*, replacement_string")
#		InputFile.close()
#	print InputFile
