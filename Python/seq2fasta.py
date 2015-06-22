#!/usr/bin/env python
# encoding: utf-8
"""
seq2fasta.py
Converts ABI sequence txt files into multi-fasta format

Created by Mark Evans on 2011-08-11.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os, glob


def main():
   path =""
   outfile = open("allseqs.txt","w")
   
   for infile in glob.glob (os.path.join(path,"*")):            # Read all filenames in directory
       if infile.find(".TXT")!=-1 : # Take all TXT files
         print "looking in ",infile
         f = open(infile,'r')
         fn = infile.split(".")
         seq = ''
         for line in f:
            if line[0] != '>':
               seq = seq +line.strip()
         outfile.write(">"+fn[0]+"\n"+seq+"\n")
         f.close()
   outfile.close()


if __name__ == '__main__':
	main()

