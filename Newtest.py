#! /usr/bin/env python

import re
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

InfileName = sys.argv[1]

for seq_record in SeqIO.parse(InfileName, "fasta"):
	print seq_record.id
	print seq_record.seq
    
"""
InFile = open(InfileName, 'r')
lines = InFile.readlines()

for line in lines:
	line = line.strip('\n')
	SearchStr = "(>gi\w+)" 
	Result = re.search(SearchStr, line)
	print Result()
"""	




#	>gi|96475|emb|A19112|Wildcat