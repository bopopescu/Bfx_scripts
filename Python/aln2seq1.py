#!/usr/bin/env python
# encoding: utf-8
"""
aln2seq1.py

Created by Mark on 2010-11-10.
Copyright (c) 2010 Monogram Biosciences. All rights reserved.
"""

import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna

filename = raw_input('Enter name of fasta file to convert from gapped to ungapped: ')
f1 = open(filename,'r')
f2 = open(filename+".clean",'w')

for record in SeqIO.parse(f1,'fasta'):
   f2.write(">"+record.id+"\n"+record.seq.tostring().replace('-','')+"\n")
   




