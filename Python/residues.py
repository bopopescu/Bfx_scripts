#!/usr/bin/env python
# encoding: utf-8
"""
residues.py

Created by Mark Evans on 2010-11-15.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys,random
import os
from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
#from Bio.SeqUtils.ProtParam import ProteinAnalysis as pa

aa = {'A':'Asp','C':'Cys','D':'Asp','E':'Glu',
      'F':'Phe','G':'Gly','H':'His','I':'Ile',
      'K':'Lys','L':'Leu','M':'Met','N':'Asn',
      'P':'Pro','Q':'Gln','R':'Arg','S':'Ser',
      'T':'Thr','V':'Val','W':'Trp','Y':'Tyr',}
nonstd_aa = {'B':'Asp or Asn','U':'Sec','Z':'Gln or Glu','X':'Unk'}
nt = {'A':'Ade','C':'Cyt','G':'Gua','T':'Thy'}
nonstd_nt = {'I':'Ino','X':'Xan','R':'A|G','U':'Uri',
             'Q':'pseudoU','Y':'C|T|U','N':'Unk',
             'W':'A|T','S':'G|C','M':'A|C','K':'G|T',
             'H':'A|C|T','B':'G|C|T','V':'G|A|C','D':'G|A|T'}


def main():
   d = raw_input("DNA filename: ")
   p = raw_input("matching Protein filename: ")
   a_seq = SeqIO.to_dict(SeqIO.parse(p,'fasta'))
   d_seq = SeqIO.to_dict(SeqIO.parse(d,'fasta'))
   t = "\t".join(["Seqid","Cohort","NT_len","AA_len","NT_N","NT_R","NT_Y","NT_W","NT_S","NT_M","NT_K","NT_H","NT_B","NT_V","NT_D","AA_X","AA_B","AA_Z"])+"\n"
   o = open(d+".freq","w")
   o.write(t)
   for sid in a_seq.keys():
      a = sid.split('|')
      NT_len = str(len(d_seq[sid].seq))
      AA_len = str(len(a_seq[sid].seq))
      NT_N = str(d_seq[sid].seq.count('N'))
      NT_R = str(d_seq[sid].seq.count('R'))
      NT_Y = str(d_seq[sid].seq.count('Y'))
      NT_W = str(d_seq[sid].seq.count('W'))
      NT_S = str(d_seq[sid].seq.count('S'))
      NT_M = str(d_seq[sid].seq.count('M'))
      NT_K = str(d_seq[sid].seq.count('K'))
      NT_H = str(d_seq[sid].seq.count('H'))
      NT_B = str(d_seq[sid].seq.count('B'))
      NT_V = str(d_seq[sid].seq.count('V'))
      NT_D = str(d_seq[sid].seq.count('D'))
      AA_X = str(a_seq[sid].seq.count('X'))
      AA_B = str(a_seq[sid].seq.count('B'))
      AA_Z = str(a_seq[sid].seq.count('Z'))
      o.write("\t".join([a[0],a[1],NT_len,AA_len,NT_N,NT_R,NT_Y,NT_W,NT_S,NT_M,NT_K,NT_H,NT_B,NT_V,NT_D,AA_X,AA_B,AA_Z])+"\n")
   o.close()

if __name__ == '__main__':
   main()

