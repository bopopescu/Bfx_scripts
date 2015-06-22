#!/usr/bin/env python
# encoding: utf-8
"""
makePWM.py
Create a log-likelihood position weight matrix from a simple occurance matrix

Created by Mark Evans on 2010-12-01.
Copyright (c) 2010 __Monogram Biosciences__. All rights reserved.
"""

import sys
import os, math
from Bio.Alphabet import generic_dna
#from Bio import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
#from Bio.Seq import Seq
#from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio import AlignIO

def main():
   aln = AlignIO.read("subD_pep.aln",'clustal')
   summ_align = AlignInfo.SummaryInfo(aln)
   con = summ_align.dumb_consensus(0.2)  # use 0.2 to avoid getting X's in the consensus
   seq_length = 404  # change this for each set of seqs
   # create absolute frequency count matrix
   abs_pssm = summ_align.pos_specific_score_matrix(con)
   f = open("abs_frq_matrix.txt","w")
   f.write(str(abs_pssm))
   f.close()
   
   # determine number of sequences in the matrix (N) for future calculations
   n = len(aln)
   freq = []
   pssm = []
   aa_keys = abs_pssm[0].keys()
   aa_keys.sort()
   ids = {'A':0, 'C':0, 'E':0, 'D':0, 'G':0, 'F':0, 'I':0, 'H':0, 'K':0, 'M':0, 'L':0, 'N':0, 'Q':0, 'P':0, 'S':0, 'R':0, 'T':0, 'W':0, 'V':0, 'Y':0}
   
   # Calculate individual frequencies with background correction of 1 to prevent future div0! errors
   # Could use smaller number than 1 if desired, like sqrt(N)
   for x in range(0,seq_length):
      freq.append(ids.copy())
      pssm.append(ids.copy())
      for k in aa_keys:
         freq[x][k] = (float(abs_pssm[x][k]) + float(1))/(float(n)+float(20))
         pssm[x][k]    = math.log(float(freq[x][k]))         
   f = open("python_pssm","w")
   f2 = open("pct_matrix","w")
   f.write("# Python generated PSSM matrix\n# According to the equations and methods from Jensen et.al. J Virol Vol(77)24 2003\n# Mark Evans 02-03-2011\t\n\t")
   for i in range(1,seq_length):
      f.write(str(i)+"\t")
      f2.write(str(i)+"\t")
   f.write(str(seq_length)+"\n")
   f2.write(str(seq_length)+"\n")
   
   # Write out new PSSM to a text file
   for k in aa_keys:
      f.write(k+"\t")
      f2.write(k+"\t")
      for x in range(0,seq_length):
         if x != seq_length - 1: 
            f.write(str(pssm[x][k])+"\t")
            f2.write(str(freq[x][k])+"\t")
         else: 
            f.write(str(pssm[x][k])+"\n")
            f2.write(str(freq[x][k])+"\n")
   f.close()
   f2.close()
   

if __name__ == '__main__':
   main()

