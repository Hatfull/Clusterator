#! /usr/bin/env python
"""This program fetches fasta files from GenBank that can be used to construct a database
"""
import re
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW

Entrez.email = "gfh@pitt.edu"

#handle = Entrez.esearch(db="nucleotide", term="(Mycobacteriophage[All Fields] AND Hatfull[All Fields]) NOT RefSeq[All Fields]", name="Mycobacteriophage")
#record = Entrez.read(handle)

#print record["Count"]




handle = Entrez.efetch(db="nucleotide", term="(Mycobacteriophage[All Fields] AND Hatfull[All Fields]) NOT RefSeq[All Fields]", name="Mycobacteriophage", rettype="fasta", retmode="text")
seq_record = SeqIO.read(handle, "fasta")
handle.close()
print "%s with %i features" % (seq_record.id, len(seq_record.features))
