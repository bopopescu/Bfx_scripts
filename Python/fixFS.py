#!/usr/bin/env python
# encoding: utf-8
"""
fixFS.py
Adds single 'A' after CTAATTTTTT sequence between GAG and TF
domains so that they will artificially translate together
as a single ORF.

Created by Mark Evans on 2011-04-29.
Copyright (c) 2011 __Monogram_Biosciences__. All rights reserved.
"""

import sys
import os,glob
from time import localtime, strftime
from Bio import SeqIO
from Bio.Seq import Seq


def main():
   path =""
   timestamp = strftime("%m%d%y_%H%M%S", localtime())
   outfile = open("GagPrt_FSfixNT_"+timestamp+".txt","w")
   outfile2 = open("GagPrt_FS_rejected_"+timestamp+".txt","w")
   outfile2.write("Could not find the sequence CTAATTTTTT in the following sequences:\n")
   outfile3 = open("GagPrt_FSfixAA_"+timestamp+".txt","w")
   
   for infile in glob.glob (os.path.join(path,"*")):            # Read all filenames in directory
       if infile.find(".exe") == -1 and infile.find(".py")==-1 and infile.find('GagPrt_FS')==-1 and infile.find('.dll')==-1: # Take everything except program files
         print "looking in ",infile
         for nt in SeqIO.parse(infile,'fasta'):
            z = nt.seq.tostring().upper()
            if z.find('CTAATTTTTT') == -1:
               outfile2.write(nt.id+"\n")
            else:
               outfile.write(">"+nt.id+"\n"+z[:z.find('CTAATTTTTT')+10]+"A"+z[z.find('CTAATTTTTT')+10:]+"\n")
               aa = Seq(z[:z.find('CTAATTTTTT')+10]+"A"+z[z.find('CTAATTTTTT')+10:])
               outfile3.write(">"+nt.id+"\n"+aa.translate().tostring()+"\n")
   outfile.close()
   outfile2.close()
   outfile3.close()
   print "\n\nProcess complete"

if __name__ == '__main__':
   main()
