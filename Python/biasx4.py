#!/usr/bin/env python
# encoding: utf-8
"""
biasx4.py

Created by Mark Evans on 2011-01-19.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys, os, itertools, math
#from decimal import Decimal
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Data import IUPACData 
from Bio.Align import AlignInfo
from Bio import AlignIO


# Takes Bio.Seq.Seq object as input
# Returns list of all possible proteins
def generateProtFromAmbiguousDNA(s):
   std_nt = CodonTable.unambiguous_dna_by_name["Standard"]  # create normal codon table object
   nonstd = IUPACData.ambiguous_dna_values                  # create list of all ambiguous DNA values, includes normal bases
   aa_trans = []
   for i in range(0,len(s),3):
      codon = s.tostring()[i:i+3]
      # For each ambiguous (or not) codon, returns list of all posible translations
      aa = CodonTable.list_possible_proteins(codon,std_nt.forward_table,nonstd) 
      aa_trans.append(aa)
      
   # Now have a list of format [ [a], [b,c], [d], [e,f,g], etc ]
   # this function creates a list of tuples containing all possible ordered combinations like
   # [(a,b,d,e), (a,b,d,f), (a,b,d,g), (a,c,d,e), (a,c,d,f), (a,c,d,g)]
   proteins = list(itertools.product(*aa_trans))
   possible_proteins = []
   for x in proteins:
      possible_proteins.append("".join(x))
   return possible_proteins


def generateX4R5_PSSM(x4alignment,r5alignment,format):
   x4aln = AlignIO.read(x4alignment,format)
   r5aln = AlignIO.read(r5alignment,format)
   x4summ_align = AlignInfo.SummaryInfo(x4aln)
   r5summ_align = AlignInfo.SummaryInfo(r5aln)
   x4con = x4summ_align.dumb_consensus(0.2)  # use 0.2 to avoid getting X's in the consensus
   r5con = r5summ_align.dumb_consensus(0.2)
   
   # create absolute frequency count matrix
   x4_abs_pssm = x4summ_align.pos_specific_score_matrix(x4con)
   r5_abs_pssm = r5summ_align.pos_specific_score_matrix(r5con)
   
   # determine number of sequences in the matrix (N) for future calculations
   x4n = len(x4aln)
   r5n = len(r5aln)
   x4_freq = []
   r5_freq = []
   pssm = []
   aa_keys = x4_abs_pssm[0].keys()
   aa_keys.sort()
   ids = {'A':0, 'C':0, 'E':0, 'D':0, 'G':0, 'F':0, 'I':0, 'H':0, 'K':0, 'M':0, 'L':0, 'N':0, 'Q':0, 'P':0, 'S':0, 'R':0, 'T':0, 'W':0, 'V':0, 'Y':0}
   
   # Calculate individual frequencies with background correction of 1 to prevent future div0! errors
   # Could use smaller number than 1 if desired, like sqrt(N)
   for x in range(0,35):
      x4_freq.append(ids.copy())
      r5_freq.append(ids.copy())
      pssm.append(ids.copy())
      for k in aa_keys:
         x4_freq[x][k] = (float(x4_abs_pssm[x][k]) + float(1))/(float(x4n)+float(20))
         r5_freq[x][k] = (float(r5_abs_pssm[x][k]) + float(1))/(float(r5n)+float(20))
         pssm[x][k]    = math.log(float(x4_freq[x][k]) / float(r5_freq[x][k]))         
   f = open("python_pssm","w")
   f.write("# Python generated PSSM matrix from toro2DX and bothR5.TN\n# According to the equations and methods from Jensen et.al. J Virol Vol(77)24 2003\n# Mark Evans 02-03-2011\t\n\t")
   for i in range(1,35):
      f.write(str(i)+"\t")
   f.write("35\n")
   
   # Write out new PSSM to a text file
   for k in aa_keys:
      f.write(k+"\t")
      for x in range(0,35):
         if x != 34: f.write(str(pssm[x][k])+"\t")
         else: f.write(str(pssm[x][k])+"\n")
   f.close()
   return pssm



def main():
   a = Seq('ATGGCARTTGTAHAC')
   print "DNA: ",a.tostring()
   print "Proteins:"
   foo = generateProtFromAmbiguousDNA(a)
   pssm = generateX4R5_PSSM("toro2DX.aln","bothR5.TN.aln","stockholm")
   for s in foo: 
      c = 0
      score = 0
      for i in list(s):
         print "\t",i," = ",pssm[c][i]
         score = score + pssm[c][i]
         c = c + 1
      print s," total score=",score
   
   
   

   # read in clustal alignments. These are stockholm format because they were used for HMMscan
#   x4aln = AlignIO.read("toro2DX.aln","stockholm")
#   r5aln = AlignIO.read("bothR5.TN.aln","stockholm")
#   print "x4aln: ",x4aln
#   print "r5aln: ",r5aln
#   x4summ_align = AlignInfo.SummaryInfo(x4aln)
#   r5summ_align = AlignInfo.SummaryInfo(r5aln)
#   x4con = x4summ_align.dumb_consensus(0.2)  # use 0.2 to avoid getting X's in the consensus
#   r5con = r5summ_align.dumb_consensus(0.2)
   # create absolute frequency count matrix
#   x4_abs_pssm = x4summ_align.pos_specific_score_matrix(x4con)
#   r5_abs_pssm = r5summ_align.pos_specific_score_matrix(r5con)
#   print "x4_abs_pssm:",x4_abs_pssm
#   print "r5_abs_pssm:",r5_abs_pssm
   # determine number of sequences in the matrix (N) for future calculations
#   x4n = len(x4aln)
#   r5n = len(r5aln)
#   print "X4 N=",x4n," R5 N=",r5n
#   x4_freq = []
#   r5_freq = []
#   pssm = []
#   aa_keys = x4_abs_pssm[0].keys()
#   aa_keys.sort()
#   ids = {'A':0, 'C':0, 'E':0, 'D':0, 'G':0, 'F':0, 'I':0, 'H':0, 'K':0, 'M':0, 'L':0, 'N':0, 'Q':0, 'P':0, 'S':0, 'R':0, 'T':0, 'W':0, 'V':0, 'Y':0}
#   x4freq_matrix = open("x4_freq.txt","w")
#   x4freq_matrix.write("\nx4_freq\n"+"\t".join(ids.keys())+"\n")
#   r5freq_matrix = open("r5_freq.txt","w")
#   r5freq_matrix.write("\nr5_freq\n"+"\t".join(ids.keys())+"\n")
#   pssm_matrix = open("x4r5_matrix.txt","w")
#   pssm_matrix.write("\nx4r5_pssm\n"+"\t".join(ids.keys())+"\n")
   
   # Calculate individual frequencies with background correction of 1 to prevent future div0! errors
   # Could use smaller number than 1 if desired, like sqrt(N)
#   for x in range(0,35):
#      x4_freq.append(ids.copy())
#      r5_freq.append(ids.copy())
#      pssm.append(ids.copy())
#      for k in aa_keys:
         #x4f = (Decimal(str(x4_abs_pssm[x][k])) + Decimal(1))/(Decimal(x4n)+Decimal(20))
#         x4f = (float(x4_abs_pssm[x][k]) + float(1))/(float(x4n)+float(20))
         #r5f = (Decimal(str(r5_abs_pssm[x][k])) + Decimal(1))/(Decimal(r5n)+Decimal(20))
#         r5f = (float(r5_abs_pssm[x][k]) + float(1))/(float(r5n)+float(20))
#         print "x= ",x," k= ",k," x4_abs_pssm=[x][k]",x4_abs_pssm[x][k]," r5_abs_pssm[x][k]=",r5_abs_pssm[x][k], " x4f=",x4f," r5f=",r5f, " log p=",math.log(float(x4f) / float(r5f))
#         x4_freq[x][k] = x4f
#         r5_freq[x][k] = r5f
#         pssm[x][k]    = math.log(float(x4f) / float(r5f))
#         print "====> x4_freq[x][k]:",x4_freq[x][k]," r5_freq[x][k]:",r5_freq[x][k]," pssm[x][k]:",pssm[x][k]
#         x4freq_matrix.write(str(x4_freq[x][k])+"\t")
#         r5freq_matrix.write(str(r5_freq[x][k])+"\t")
#         pssm_matrix.write(str(pssm[x][k])+"\t")
#      x4freq_matrix.write("\n")
#      r5freq_matrix.write("\n")
#      pssm_matrix.write("\n")

#   x4freq_matrix.close()
#   r5freq_matrix.close()
#   pssm_matrix.close()
#   print "PSSM:\n",pssm
   


if __name__ == '__main__':
   main()

