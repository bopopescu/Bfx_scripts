#!/usr/bin/env python
# encoding: utf-8
"""
yank.py

Created by Mark Evans on 2010-09-13.
Copyright (c) 2010 __Monogram Biosciences__. All rights reserved.
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
   code =  "\n\n#############\nYank utility to extract fasta sequences out of a multiple sequence file that correspond to sequence ids found in an id list\n"
   code += "Usage: $ python yank.py [-i] idlist_filename [-s] sequence_filename\n"
   code += "-i [idlist_filename] list of ids, one id per line\n"
   code += "-s [sequence_filename] name of file containing multiple fasta sequences\n"
   code += "-h help\n#############\n\n"
   print code
   return

########################################
# Read all ids from id list into array #
# and pass back to calling fxn         #
########################################
def loadIDs(filename):
   idlist = []
   try:
      f = open(filename,'r')
   except IOError:
      print "*****\nCan\'t open id file "+filename+" for reading.\n*****\n"
      sys.exit(0)
   for line in f.readlines():
      idlist.append(line.rstrip())
   f.close()
   return idlist
   
########################################
# Read through seq fasta file once and #
# if id is in array, print to stdout   #
########################################   
def processFiles(idfile,seqfile):
   ids = loadIDs(idfile)
   for record in SeqIO.parse(seqfile,'fasta',generic_dna):
      a = record.id.split('|')
      if a[0] in ids:
         sys.stdout.write(">"+record.id+'\n'+record.seq.tostring()+'\n')
   sys.exit()
   

#################################################
# Take a file of sequence id's, one per row and #
# yank out corresponding fasta sequence from a  #
# multi sequence fasta file                     #
#################################################   
def main(argv):
   idfile=''
   seqfile=''
   
   # Read command line arguments
   try:
      opts, args = getopt.getopt(argv,"i:s:h",["input,seq,help"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt in ("-h","--help"):
         usage()
         sys.exit()
      elif opt in ("-i","--input"): idfile = arg
      elif opt in ("-s","-seq"): seqfile = arg
   
   if idfile !='' and seqfile !='':
      processFiles(idfile,seqfile)
   else:  # if cmd line flags were not used, check for right args
      if len(sys.argv) == 3 and sys.argv[1] !='' and sys.argv[2] !='':
         idfile = sys.argv[1]
         seqfile = sys.argv[2]
         processFiles(idfile,seqfile)
      else:
         usage()
         sys.exit()

#### Called by cmd line, launches script and passes arguemnts       
if __name__ == "__main__":
   main(sys.argv[1:])
   

