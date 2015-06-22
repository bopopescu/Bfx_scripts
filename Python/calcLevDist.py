#!/usr/bin/env python
# encoding: utf-8
"""
calcLevDist.py

Calculates Levenshtein distance between two groups of sequences
Created by Mark Evans on 2011-02-09.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import Levenshtein as L
from Bio import SeqIO


def main():
   x4s = SeqIO.to_dict(SeqIO.parse("seqX4.seq","fasta"))
   r5s = SeqIO.to_dict(SeqIO.parse("seqR5.seq","fasta"))
   dms = SeqIO.to_dict(SeqIO.parse("seqDM.seq","fasta"))
   best = {}
   
   for ids in dms.keys():
      for r5 in r5s.keys():
         lev = L.distance(dms[ids].seq.tostring(),r5s[r5].seq.tostring())
         if not best.has_key(lev):
            best[lev] = {'x4':ids,'r5':r5}
   r = best.keys()
   r.sort()
   print "The distances ranged across ",len(r)
   print "Lowest distance=",r[0]," farthest distance=",r[len(r)-1]
#   print str(best)
   print "X4:",dms[best[r[0]]['x4']]," R5:",r5s[best[r[0]]['r5']]
   


if __name__ == '__main__':
   main()

