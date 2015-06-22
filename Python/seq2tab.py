#!/usr/bin/env python
# encoding: utf-8
"""
seq2tab.py
Converts fasta file to tab-text file
Created by Mark Evans on 2011-01-07.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys,getopt
import os
from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO


##################################
# Define usage syntax for script #
##################################   
def usage():
   code =  "\n\n#############\nSeq2tab utility to convert fasta sequences out of a multiple sequence file into a tab-delimited text file\n"
   code += "Usage: $ python seq2tab.py [-i] input_fasta_filename [-o] output_tab_filename\n"
   code += "-i [input_sequence_filename] multi-sequence fasta format\n"
   code += "-o [output_tab_filename] name of tab file that will be created for the sequences\n"
   code += "-h help\n#############\n\n"
   print code
   return

   
########################################
# Read through seq fasta file once and #
# if id is in array, print to stdout   #
########################################   
def processFiles(seqfile,tabfile):
   f = open(tabfile,"w")
   for record in SeqIO.parse(seqfile,'fasta',generic_dna):
      if record.id.find('|') != -1:
         defs = record.id.split()
         f.write("\t".join(defs)+"\t"+record.seq.tostring()+"\n")
      else:
         f.write(record.id+"\t"+record.seq.tostring()+"\n")
   f.close()
   sys.exit()
   
   
def main(argv):
   tabfile=''
   seqfile=''
   
   # Read command line arguments
   try:
      opts, args = getopt.getopt(argv,"i:o:h",["input,output,help"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt in ("-h","--help"):
         usage()
         sys.exit()
      elif opt in ("-i","--input"): seqfile = arg
      elif opt in ("-o","-output"): tabfile = arg
   
   if seqfile !='' and tabfile !='':
      processFiles(seqfile,tabfile)
   else:  # if cmd line flags were not used, check for right args
      if len(sys.argv) == 3 and sys.argv[1] !='' and sys.argv[2] !='':
         seqfile = sys.argv[1]
         tabfile = sys.argv[2]
         processFiles(seqfile,tabfile)
      else:
         usage()
         sys.exit()


if __name__ == '__main__':
   main(sys.argv[1:])

