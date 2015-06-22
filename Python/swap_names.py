#!/usr/bin/env python
# encoding: utf-8
"""
swap_names.py

Created by Mark Evans on 2011-08-01.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import AlignInfo
from Bio import AlignIO



def main():
   indexfile = open('indexfile.txt','r')
   for line in indexfile:
      files = line.split()
#      print "Seqfile name= ",files[0]," and aln file= ",files[1]
      seqs = SeqIO.to_dict(SeqIO.parse(files[0],'fasta'))
#      print "seqs = "+str(seqs)
      align = AlignIO.read(files[1],'clustal')
#      print "align= "+str(align)
      seqnames = seqs.keys()
#      print "seqnames = "+str(seqnames)
      name_idx ={}
      for s in seqnames:
#         n = s.split()
#         print "s = ",s," and full desc= ",seqs[s].description
         name_idx[s] = seqs[s].description
#      print "name_idx = "+str(name_idx)
      aln_dict = {}
      for x in range(0,len(align)):
         aln_dict[align[x].id] = x
#      print "aln_dict = "+str(aln_dict)
      for sname in name_idx:
#         print "sname = ",sname
         if aln_dict.has_key(sname): align[aln_dict[sname]].id = name_idx[sname]
      
#      print "new align should be "+str(align)
      newalign = open('new_'+files[1],"w")
      AlignIO.write(align,newalign,'clustal')
      newalign.close()
      

if __name__ == '__main__':
   main()

