#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Mark Evans on 2010-11-11.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

def main():
   seqfile = raw_input("What is the sequence filename: ")
   hmmfile = raw_input("What is the HMM table filename: ")
   
   h = open(hmmfile,'r')
   d = {}
   c = 1
   for line in h:
      if c > 3:
         val = line.split()
         v3id = val[3] # accession
         hfrom=int(val[15]) # hmm start pos
         hto = int(val[16]) # hmm end pos
         sfrom = int(val[17]) # sequence match start pos
         sto = int(val[18]) # sequence match end pos
         acc = float(val[21]) # accuracy
         if d.has_key(v3id):
            if d[v3id]['acc'] < acc:
               d[v3id] = {'hfrom':hfrom,'hto':hto,'sfrom':sfrom,'sto':sto,'acc':acc}
         else: d[v3id] = {'hfrom':hfrom,'hto':hto,'sfrom':sfrom,'sto':sto,'acc':acc}
      else: c+=1
   o = open(seqfile+".trim","w")
   log = open(seqfile+".log","w")
   log.write("mgrm_acc\tparsed_seq\tbuffer_seq\th_from\th_to\ts_from\ts_to\tacc\n")
   for record in SeqIO.parse(seqfile,'fasta'):
      try:
         seq = record.seq.tostring()[d[record.id]['sfrom']-1:d[record.id]['sto']]
         buff = record.seq.tostring()[d[record.id]['sfrom']-11:d[record.id]['sto']+10]
         o.write(">"+record.id+"\n"+seq+"\n")
         log.write(record.id+"\t"+seq+"\t"+buff+"\t"+str(d[record.id]['hfrom'])+"\t"+str(d[record.id]['hto'])+"\t"+str(d[record.id]['sfrom'])+"\t"+str(d[record.id]['sto'])+"\t"+str(d[record.id]['acc'])+"\n")
      except:
         log.write(record.id+"\tNo Match\n")
   print "\n Process complete\n"
   
if __name__ == '__main__':
   main()

