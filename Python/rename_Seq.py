#!/usr/bin/env python
# encoding: utf-8
"""
rename_Seq.py

Created by Mark Evans on 2010-11-29.
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
   seqfile = raw_input("Name of sequence file to rename: ")
   lookupfile = raw_input("Name of lookup file that has new and old names: ")
   oldcol = raw_input("Which zero-based column holds old id? ")
   newcol = raw_input("Which zero-based column holds new id? ")
   seqs = SeqIO.to_dict(SeqIO.parse(seqfile,'fasta'))
   f = open(seqfile+".renamed","w")
   lf = open(lookupfile,'r')
   
   ids ={}
   for line in lf:
      b = line.split()
      ids[b[int(oldcol)]] = b[int(newcol)]
      
   print str(ids)
   keys = seqs.keys()
   for k in keys:
      a = k.split('|')
      print str(a)
      if seqs.has_key(k):
         try:
            print "Replacing ",a[0]," with ",ids[a[0]]
            f.write(">"+ids[a[0]]+"\n"+seqs[k].seq.tostring()+"\n") # +"|"+a[1]+"|"+a[2]+"|"+a[3]
         except:
            print "####### WARNING ########: ID ",a[0]," was not found in id list"
      else:
         print "******** Warning ******: Could not replace ",a[0]," key did not exist in seqs"
         
   lf.close()
   f.close()
   
if __name__ == '__main__':
   main()

