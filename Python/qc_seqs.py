#!/usr/bin/env python
# encoding: utf-8
"""
qc_seqs.py
Check length, duplicates, non-amino acid characters

Modified to work with LANL data 08.05.11

Created by Mark on 2010-09-10.
Copyright (c) 2010 Monogram Biosciences. All rights reserved.
"""
from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import os


filename = raw_input('What is the name of file to evaluate? ')
size = raw_input('Select all peptides of length: ')
dna = raw_input('Are these dna seqs? (y/n) ')
int(size)
f1 = open(filename+'.clean.nodups','w')
f2 = open(filename+'.dups','w')
f3 = open(filename+'.garbage','w')
size_range={}
nodups = {}

# search for duplicate sequences and remove them
for record in SeqIO.parse(filename,"fasta",generic_dna):
   a = record.id.split('.')
   b = record.seq.tostring().upper()
   tropism =''
   if 'CCR5' in record.id.upper(): tropism = 'R5'
   elif 'CXCR4' in record.id.upper(): tropism = 'X4'
   if 'X' in b or ' ' in b or 'B' in b or 'J' in b or 'U' in b or 'Z' in b:
      f3.write('>'+'.'.join(a)+'|'+tropism+'\n'+b+"\n")
   elif nodups.has_key(record.seq.tostring()):
      f2.write(">"+".".join(a)+"|"+tropism+"\n"+record.seq.tostring()+"\n")
      f2.write(">"+".".join(nodups[record.seq.tostring()]['id'])+"|"+tropism+"\n"+record.seq.tostring()+"\n")
      del nodups[record.seq.tostring()]
   else: nodups[record.seq.tostring()] = {'id':a,'record':record,'tropism':tropism}

f2.close()
f3.close()

seqs = nodups.keys()
seqs.sort()
c = 0

# check for name dups, modify them and let user know
namedups={}
dups=[]
for v in seqs:
   if namedups.has_key(nodups[v]['record'].id): 
      dups.append(nodups[v]['record'].id)
      old = nodups[v]['record'].id
      nodups[v]['id'].append(str(c))
      nodups[v]['record'].id = nodups[v]['record'].id+"."+str(c)
      print "old= ",old," new= ",str(nodups[v]['id'])," record= ",nodups[v]['record'].id
      c += 1
   else: namedups[nodups[v]['record'].id] = v
   
   
   
# search for sequences of length 'size'
if dna == 'y':
   for v in seqs:
      v2 = nodups[v]['record'].seq.translate().tostring()
      if size_range.has_key(len(v2)): size_range[len(v2)] += 1
      else: size_range[len(v2)] = 1 
      if len(v2) == int(size):
         f1.write(">"+".".join(nodups[v]['id'])+"|"+nodups[v]['tropism']+"\n"+ v +"\n")
else:
   for v in seqs:
      if size_range.has_key(len(v)): size_range[len(v)] += 1
      else: size_range[len(v)] = 1 
      if len(v) == int(size):
         f1.write(">"+".".join(nodups[v]['id'])+"|"+nodups[v]['tropism']+"\n"+ v +"\n")
   
f1.close()

# check for name dups and let user know
#namedups={}
#print str(nodups)
#dups=[]
#for v in seqs:
 #  print str(nodups[v][0])
#   if namedups.has_key(nodups[v][0]): dups.append(nodups[v][0])
#   else: namedups[nodups[v][0]]=v     
print "Looking for size="+str(size)+", length distribution in file given [length:qty] is "+str(size_range)+'\n'
if len(dups) == 0: print "No name duplicates in file\n"
elif len(dups) != 0: print "Name duplicates found....\n"+str(dups)+'\nnFinished\n'

      

