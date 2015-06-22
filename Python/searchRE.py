#!/usr/bin/env python
# encoding: utf-8
"""
searchRE.py
Intended to search a folder of DNA sequence files for files containing
the restriction enzymes in res.  Reports the name of the files that have
the enzymes, the number of times they occur and where and if they occur
in a location containing mixtures.

Created by Mark Evans on 2010-12-20.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os, glob

nonstd_nt = {'R':['A','G'],'Y':['C','T'],'W':['A','T'],'S':['G','C'],'M':['A','C'],'K':['G','T'],
             'H':['A','C','T'],'B':['G','C','T'],'V':['G','A','C'],'D':['G','A','T']}
res = {'SbfI':"CCTGCAGG", "XbaI":"TCTAGA", "SpeI":"ACTAGT", "AflII":"CTTAAG", "XhoI":"CTCGAG", "ScaI":"AGTACT", "AseI":"ATTAAT"}


def main():
   path = ""
   outfile = open("enzyme_search_result.txt","w")
   outfile.write("Filename\tSites found\tNumber Sites found\tNumber in mixtures\n")
   for infile in glob.glob (os.path.join(path,"*")):                                            # Read all filenames in directory
       if infile.find(".exe") == -1 and infile.find("enzyme") == -1 and infile.find(".pl")==-1: # Take everything except program files
          f = open(infile,"r")
          tempseqs={}
          positions={}
          orig_pos={}
          mix_pos={}
          occurances={'SbfI':0, "XbaI":0, "SpeI":0, "AflII":0, "XhoI":0, "ScaI":0, "AseI":0}
          
          # Load original and all created sequence variants into a dictionary for processing
          ##################################################################################
          for seq in f:
             seq = seq.upper()
             tempseqs['orig'] = seq
             c = 1
             for nt in nonstd_nt:
                if seq.find(nt) != -1:
                   for base in nonstd_nt[nt]:
                      tempseqs['seq'+str(c)] = seq.replace(nt,base)
                      c = c + 1
          f.close()
          seqkeys = tempseqs.keys()
          seqkeys.sort()                    # Put 'orig' first, all seq1,seq2,etc later

          # Go through all enzymes and check for occurance
          ################################################
          for enz in res:
             # Go through all sequence variants 
             for k in seqkeys:
                pos = 0
                while (1):
                   pos = tempseqs[k].find(res[enz],pos)
                   if pos == -1 : break
                   if k =='orig': 
                      orig_pos[pos] = enz
                      positions[pos] = enz
                      occurances[enz] = occurances[enz] + 1
                   else:
                      if not orig_pos.has_key(pos):
                         mix_pos[pos]=enz
                         positions[pos]=enz
                         occurances[enz] = occurances[enz] + 1
                   pos = pos+1
                   
          # Process the search results for printing
          # Reorder things by enzyme rather than by position
          ##################################################
          if len(positions) != 0:
             enz = {}
             for pos in positions:
                if enz.has_key(positions[pos]):
                   enz[positions[pos]].append(pos)
                else: enz[positions[pos]] = [pos]
             k = enz.keys()
             k.sort()
             enzymes = ",".join(k)
             sites =""
             mix = ""
             for x in k: sites = sites+str(occurances[x])+","
             sites =sites[:-1]
             if len(mix_pos) != 1:
                for p in mix_pos: mix = mix+mix_pos[p]+":"+str(p)+","
             mix = mix[:-1]
             outfile.write(infile+"\t"+enzymes+"\t"+sites+"\t"+mix+"\n")
          
   outfile.close()

if __name__ == '__main__':
   main()

