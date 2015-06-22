#!/usr/bin/env python
# encoding: utf-8
"""
mgrmdb2fasta.py

Input is tab-delimited file with following database columns:
report_id, sub_type, aliquot_d, region_name, nt_sequence

Created by Mark on 2012-02-02.
Copyright (c) 2012 Monogram Biosciences. All rights reserved.
"""

import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna

filename = raw_input('Enter name of file to convert to .seq: ')
f1 = open(filename,'r')
f2 = open(filename+".ns3.dna.seq",'w')
f3 = open(filename+".ns4a.dna.seq",'w')

for record in f1:
	a = record.rstrip().split('\t')
	if a[3].upper() == 'NS3':
		f2.write('>'+a[2]+'|'+a[3]+'|'+a[1]+'|'+a[0]+'\n'+a[4].upper()+'\n')
	elif a[3].upper() == 'NS4A':
		f3.write('>'+a[2]+'|'+a[3]+'|'+a[1]+'|'+a[0]+'\n'+a[4].upper()+'\n')
