#!/usr/bin/env python
# encoding: utf-8
"""
compress_align.py

This script is designed to take aligned sequences and compress or expand them to be exactly
35 amino acids in length.  Clustalw2 is used to generate a pairwaise alignment of any sequence
that is not already 35 amino acids long.  If the sequence is short, Clustalw2 will pad gaps
with the '-' character.  If the sequence is longer, the insertion region (based on alignment 
to a 35-mer consensus) will be removed and a 'Z' inserted in it's place.

Created by Mark Evans on 2010-11-19.
Copyright (c) 2010 __Monogram Biosciences__. All rights reserved.
"""
import sys,os,subprocess,random
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

consensus = 'CTRPNNNTRKSIXIGPGRAFYATGXIIGDIRQAHC'  # Consensus of all mgrm 35mers @ 50% occurance

def main():
   filename = raw_input("Enter the name of the mutiple sequence file (.seq)? ")
   a_seq = SeqIO.to_dict(SeqIO.parse(filename,'fasta'))  # read all sequences directly to dictionary
   newfile = open(filename+".len35","w")
   for seqid in a_seq.keys():
      tmp = str(random.random())[2:20]  # generate random temp filename
      f = open(tmp+".seq",'w')
      f.write(">C\n"+consensus+"\n>Q\n"+a_seq[seqid].seq.tostring())
      f.close()
      
      # Create command line request for clustalw
      ##########################################
      cline = ClustalwCommandline("clustalw2",infile=tmp+".seq", outfile=tmp+".aln")
      child = subprocess.call(str(cline),shell=(sys.platform!="win32"))  # shell part of this statement is required to run properly. Go figure
      align = AlignIO.read(tmp+".aln",'clustal')
      summary_align = AlignInfo.SummaryInfo(align)
      seqs = {align[0].id:align[0].seq.tostring(),align[1].id:align[1].seq.tostring()}
      
      # Look for insertions
      ###############################
      if seqs['C'].find('-'):                                                 # Look for insertions in query "Q" by looking for '-' in consensus "C"
         gapnum = seqs['C'].count('-')                                          # count num of gap chars
         length = len(seqs['C'][seqs['C'].find('-'):seqs['C'].rfind('-')+1])    # measure length from first gap char to last gap char
         if gapnum != 0:                                                        # if we have an insertion event, gapnum will not be 0
            if gapnum == length:                                                  # if gapnum == length, then we only have a single insertion = simple case
               start = seqs['C'].find('-') - 1                                      # Add one to grab the left side "true" residue as well, since we will replace with Z
               end = seqs['C'].rfind('-') + 1                                       # Add one because zero based
               seqs['Q'] = seqs['Q'][:start]+"Z"+seqs['Q'][end:]                    # create substituted query sequence
            else:                                                                 # multiple insertion events 
               s = list(seqs['C'])                                                  # convert "C" sequence into a list
               indices = []                                                         # store list indices of '-' occurance
               flag='no'                                                            # are we in the middle of a gap?
               for i,v in enumerate(s):                                             # i = list index, v = list value
                  if v =='-':
                     if flag == 'no': indices.append(i)
                     flag = 'yes'
                  elif v != '-' and flag =='yes':
                     indices.append(i-1)
                     flag = 'no'
               z = list(seqs['Q'])                                                # convert "Q" sequence into a list
               for i in range(0,len(indices),2):                                  # loop through indices of gaps
                  z[indices[i]-1]='Z'                                               # set residue next to insertion = 'Z'
                  for x in range(indices[i],indices[i+1]+1):                        # loop through length of insertion
                     z[x]="|"                                                         # set each value in the insertion range to a 'throw-away' string
               newseq = "".join(z)                                                # Join the list back into a string
               newseq = newseq.replace("|","")                                    # remove all throw-away chars
               seqs['Q'] = newseq
      os.remove(tmp+".seq")                                                   # remove all temp files
      os.remove(tmp+".aln")
      os.remove(tmp+".dnd")
      newfile.write(">"+seqid+"\n"+seqs['Q']+"\n")
   
   newfile.close()
   
if __name__ == '__main__':
   main()