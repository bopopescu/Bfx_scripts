#!/usr/bin/env python
# encoding: utf-8
"""
codon_generator.py

Created by Mark Evans on 2011-02-11.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os,itertools
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Data import IUPACData



def main():
   bases = ['T','C','A','G']
   ambig_bases = ['R','Y','S','W','K','M','D','H','B','V','N']
   codons = [a+b+c for a in bases for b in bases for c in bases]
   amino_acids = "F F L L S S S S Y Y stop stop C C stop W L L L L P P P P H H Q Q R R R R I I I M T T T T N N K K S S R R V V V V A A A A D D E E G G G G".split(' ')
   codon_table = dict(zip(codons, amino_acids))
   nonstd = IUPACData.ambiguous_dna_values                  # create list of all ambiguous DNA values, includes normal bases
   std_nt = CodonTable.unambiguous_dna_by_name["Standard"]  # create normal codon table object
   
   all_ambig_trip = []
   for s in list(itertools.product(*[ambig_bases,ambig_bases,ambig_bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[bases,ambig_bases,ambig_bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[ambig_bases,bases,ambig_bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[ambig_bases,ambig_bases,bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[ambig_bases,bases,bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[bases,ambig_bases,bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[bases,bases,ambig_bases])): all_ambig_trip.append(s)
   
   k = codon_table.keys()
   k.sort()
   for b in k:
      print "codon_table:",b," aa:",codon_table[b]
      
   stop_set ={}  # store final results of ambig codon and all possible translations minus the potential stop
   for trip in all_ambig_trip:
      stop_seen = False
      ambig = "".join(trip)
#      print "Ambig = ",ambig
      stop_set[ambig] = []
      for codon in list(itertools.product(*[ list(nonstd[trip[0]]), list(nonstd[trip[1]]), list(nonstd[trip[2]]) ])):  # translate ambig into all possible real combinations
         codon = "".join(codon)
         aa = codon_table[codon]
#         print "\tcodon: ",codon," aa: ",aa
         if aa == 'stop': stop_seen = True
         else: stop_set[ambig].append(aa)
      if stop_seen is False: 
         del(stop_set[ambig])
         try:
            CodonTable.list_possible_proteins(ambig,std_nt.forward_table,nonstd)
         except KeyError, err:
            print "Sanity check failed for: ",ambig," with Key error. stop_seen is ",str(stop_seen)," ",err
         except CodonTable.TranslationError, TransError:
            print "Sanity check failed for: ",ambig," with Translation error. stop_seen is ",str(stop_seen)," ",TransError
   for x in stop_set:
      print "'"+x+"':",stop_set[x],","
   print "length of stop set is ",len(stop_set)
   



if __name__ == '__main__':
   main()

