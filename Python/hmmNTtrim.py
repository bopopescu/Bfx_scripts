#!/usr/bin/env python
# encoding: utf-8
"""
hmmNTtrim.py

Created by Mark Evans on 2010-11-11.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys,random
import os
from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

def getReadingFrames(seq):
   s = {}
   for i in range(0,3):
      s[i] = seq[i:].translate().tostring()
      s[i+3] = seq.reverse_complement()[i:].translate().tostring()
   return s

def verifyCorrectFrame(seq):
   bin_path="/usr/local/bin/hmmscan"  # change depending on system being run on
   hmm_tmp = str(random.random())[2:20]  # generate random temp filename
   f = open(hmm_tmp+".seq.tmp","w")
   for i in range(0,3):
      f.write(">frame"+str(i)+"\n"+seq[i]+"\n")
      f.write(">frame"+str(i+3)+"\n"+seq[i+3]+"\n")
   f.close()
   os.system(bin_path+" -o "+hmm_tmp+".log.tmp --noali --domtblout "+hmm_tmp+".out.tmp mgrmv3.hmm "+hmm_tmp+".seq.tmp")
   
   h = open(hmm_tmp+".out.tmp",'r')
   d = {}
   c = 1
   for line in h:
      if c > 3:
         val = line.split()
         seqid = val[3]       # accession
         hfrom=int(val[15])   # hmm start pos
         hto = int(val[16])   # hmm end pos
         sfrom = int(val[17]) # sequence match start pos
         sto = int(val[18])   # sequence match end pos
         acc = float(val[21]) # accuracy
         if d.has_key(seqid):
            if d[seqid]['acc'] < acc:
               d[seqid] = {'hfrom':hfrom,'hto':hto,'sfrom':sfrom,'sto':sto,'acc':acc}
         else: d[seqid] = {'hfrom':hfrom,'hto':hto,'sfrom':sfrom,'sto':sto,'acc':acc}
      else: c+=1
   h.close()
   os.system("rm "+hmm_tmp+".*")
   return d
   
   
def parseDNA(hmm,record):
   err = ""
   f = hmm.keys()
   if len(f) != 0:
      frame = int(f[0][5])
      if frame < 3:  # forward reading frames
         try:
            dna  = record.seq[frame+(int(hmm[f[0]]['sfrom'])*3)-3:frame+(int(hmm[f[0]]['sto'])*3)].tostring()
            prot = record.seq[frame+(int(hmm[f[0]]['sfrom'])*3)-3:frame+(int(hmm[f[0]]['sto'])*3)].translate().tostring()
            print prot
         
         except:
            dna = "No Match"
            prot = "No Match"
            err += "No Match for "+record.id+"\t"+record.seq.tostring()+"\n"
      elif frame >2:  # reverse reading frames
         try:
            dna = record.seq.reverse_complement()[frame-3+int(hmm[f[0]]['sfrom'])-3:frame-3+hmm[f[0]]['sto']].tostring()
            prot = record.seq.reverse_complement()[frame-3+int(hmm[f[0]]['sfrom'])-3:frame-3+hmm[f[0]]['sto']].translate().tostring()
         except:
            dna = "No Match"
            prot = "No Match"
            err += "No Match for "+record.id+"\t"+record.seq.tostring()+"\n"
   else:
      err += "FAILURE: No HMM Keys for "+record.id+"\t"+record.seq.tostring()+"\n"
      print "FAILURE: No HMM Keys for "+record.id+"\t"+record.seq.tostring()+"\n"
      dna = "No Match"
      prot = "No Match"
   return (dna,prot,err)
   
   
   
def main():
   seqfile = raw_input("What is the NT sequence filename: ")
   o = open(seqfile+".trim.dna","w")
   p = open(seqfile+".trim.prot","w")
   log = open(seqfile+".dnalog","w")
   for record in SeqIO.parse(seqfile,'fasta'):
      frames = getReadingFrames(record.seq)
      hmm = verifyCorrectFrame(frames)
      dna,prot,err = parseDNA(hmm,record)
#      print "dna: ",dna," Prot: ",prot
      if prot != "No Match":
         o.write(">"+record.id+"\n"+dna+"\n")
         p.write(">"+record.id+"\n"+prot+"\n")
      log.write(err)
   o.close()
   p.close()
   log.close()
   
   
if __name__ == '__main__':
   main()

