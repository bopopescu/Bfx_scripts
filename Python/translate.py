#!/usr/bin/env python
# encoding: utf-8
"""
translate.py
Translates all sequences in a multi-fasta file to protein
Assumes reading frame 1

Created by Mark Evans on 2011-08-11.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Bio import SeqIO

def main():
   infile = raw_input("Enter name of file to translate: ")
   f = open(infile+".aa.seq","w")
   for record in SeqIO.parse(infile,'fasta'):
#      record.seq = record.seq.translate()
#      SeqIO.write(record,f,'fasta')      # works, but leaves in stop codon *
      f.write(">"+record.id+"\n"+record.seq.translate().tostring()+"\n")
   f.close()

if __name__ == '__main__':
   main()

