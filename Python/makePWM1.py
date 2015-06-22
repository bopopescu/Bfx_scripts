#!/usr/bin/env python
# encoding: utf-8
"""
makePWM.py
Create a log-likelihood position weight matrix from a simple occurance matrix

Created by Mark Evans on 2010-12-01.
Copyright (c) 2010 __Monogram Biosciences__. All rights reserved.
"""

import sys
import os
from Bio.Alphabet import generic_dna
#from Bio import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
#from Bio.Seq import Seq
#from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio import AlignIO

def main():
   aln_file = raw_input("What is the name of the .aln file to convert? ")
   align = AlignIO.read(aln_file,'clustal')
   summary_align = AlignInfo.SummaryInfo(align)
   consensus = summary_align.dumb_consensus()
   pssm1 = summary_align.pos_specific_score_matrix(consensus, chars_to_ignore = ['X'])
   

if __name__ == '__main__':
   main()

