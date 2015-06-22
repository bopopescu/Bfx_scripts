#!/usr/bin/env python
# encoding: utf-8
"""
parseLANL.py

Created by Mark Evans on 2011-08-05.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Bio import SeqIO

def main():
   infile = raw_input("what is name of LANL DNA file: ")
   dna = SeqIO.to_dict(SeqIO.parse(infile,'fasta'))
   out = open(infile+".table",'w')
   for x in dna:
      a = x.split('|')
      tropism = a[-1]
      del(a[-1])
      b = a[0].split(".")
      subtype = b[0]
      out.write("\t".join([x,tropism,subtype,dna[x].seq.tostring(),dna[x].seq.translate().tostring(),a[0]])+"\n")
   out.close()

if __name__ == '__main__':
   main()

