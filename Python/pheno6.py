#!/usr/bin/env python
# encoding: utf-8
"""
pheno5.py

Created by Mark Evans on 2010-09-13.
Copyright (c) 2010 __Monogram Biosciences__. All rights reserved.

This program is intended to only work with output from HMMER 3.0,
using hmmscan.  DO NOT USE WITH HMMPFAM OUTPUT
"""

import sys
import os,math
from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO


#############################################################################################################
# assocCoeff                                                                                                #
# Function to calculate association coefficient, true pos percent, false pos percent, and accuracy          #
# from a dictionary of data                                                                                 #
# @input: dictionary, index of actual call, index of data to be evaluated, pheno to be predicted (DM or R5) #
# @return dictionary of cutoff:[tn,fn,tp,fp,phi,tpp,fpp,acc]                                                #
#############################################################################################################
def assocCoeff(data,ac,sc,predictor):
   # ac = actual call index, sc = stat calculated call index, predictor = pheno a positive score predicts
   # set variables
   acstats={}
   pheno = {'DM':['DM','R5'],
            'R5':['R5','DM']}
   for x in range(50,100):
      fp=0.0
      fn=0.0
      tp=0.0
      tn=0.0
      pid=''
      for v3id in data:
         # if predictor matches mgrm call, this is the positive direction (tp,fp)
         if predictor == 'R5':
            # check if score is below current cutoff x and call pheno for this round
            if float(data[v3id][predictor]) <= x: pid=pheno[predictor][1]
            else: pid=pheno[predictor][0]
         elif predictor == 'DM':
            if float(data[v3id][predictor]) >= x: pid=pheno[predictor][0]
            else: pid=pheno[predictor][1]

         if   pid == 'DM' and pid == data[v3id][ac]: tp = tp + 1     # true positive call
         elif pid == 'DM' and pid != data[v3id][ac]: fp = fp + 1     # false positive
         elif pid == 'R5' and pid == data[v3id][ac]: tn = tn + 1     # true negative
         elif pid == 'R5' and pid != data[v3id][ac]: fn = fn + 1     # false negative

         pid=''
      try:
         # phi = sqrt(((tp*tn - fp*fn)^2/(tp+fp)(fn+tn)(tp+fn)(fp+tn)))  association coefficient 
         float(phi)
     #    phi = math.sqrt((math.pow(tp*tn - fp*fn,2)/(tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
         phi = math.sqrt(math.pow(tp*tn - fp*fn,2)/((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
      except:
         phi = 0.0
      try: tpp = tp/(tp+fn)  # true positive percentage
      except: tpp = 0.0
      try: fpp = fp/(fp+tn)  # false positive percentage
      except: fpp = 0.0
      try: acc = (tp+tn)/(tp+tn+fp+fn)  # Total accuracy for prediction versus reality
      except: acc = 0.0
      try: spec = tn/(tn+fp)  # specificity for X4
      except: spec = 0.0
      try: sens = tp/(tp+fn)  # sensitivity for X4
      except: sens = 0.0
      acstats[str(x)] = {'tn':tn,'fn':fn,'tp':tp,'fp':fp,'phi':phi,'tpp':tpp,'fpp':fpp,'acc':acc,'spec':spec,'sens':sens}
   return (acstats)


def doHMMScan(filename):
   bin_path="/usr/local/bin/hmmscan"  # change depending on system being run on
#   os.system(bin_path+" --noali --tblout "+filename+".hmmresult all.hmm "+filename)
   os.system(bin_path+" --noali --tblout "+filename+".hmmresult all.hmm "+filename)
   return filename+".hmmresult"
   
   
def processHMMresult(filename):   
   f = open(filename,'r')
   hmmR5score = {}
   hmmDXscore = {}
   hmlookup = {'bothDX':'DM','toro2DX':'DM','toro2R5':'R5','bothR5.all':'R5'}
   for x in f.readlines():
      x.rstrip()
      a = x.split()
      if '#' not in x:
         b = a[2].split('|')
         if hmlookup[a[0]] == 'R5':
            hmmR5score[b[0]] = {hmlookup[a[0]]:a[5],'cohort':b[1],'trofile':b[2],'es':b[3]}
         elif hmlookup[a[0]] == 'DM':
            hmmDXscore[b[0]] = {hmlookup[a[0]]:a[5],'cohort':b[1],'trofile':b[2],'es':b[3]}
   if len(hmmR5score) == len(hmmDXscore):
      print "Success!! number of sequences in each array is the same\n"
      keys = hmmR5score.keys()
      merged = {}
      for z in keys:
         call=''
         if hmmR5score[z]['R5'] != hmmDXscore[z]['DM'] and hmmR5score[z]['R5'] > hmmDXscore[z]['DM']: call ='R5'
         else: call = 'DM'
         mgrm =''
         if hmmR5score[z]['trofile'] == 'DM' or hmmR5score[z]['trofile']=='X4' or hmmR5score[z]['es']=='DM' or hmmR5score[z]['es']=='X4': mgrm = 'DM'
         else: mgrm = 'R5'
         merged[z] = {'call'   :call,
                      'mgrm'   :mgrm,
                      'R5'     :hmmR5score[z]['R5'],
                      'DM'     :hmmDXscore[z]['DM'],
                      'cohort' :hmmR5score[z]['cohort']}
   return merged
   
   
   
def main():
   filename = raw_input("Enter the name of the sequence file (.seq)? ")
   result_file = doHMMScan(filename)
   hmmresults = processHMMresult(result_file)
   
   # Statistically Analyze the HMM output 
   ###########################################
   v = [('bothR5',1,'R5'),('bothDX',2,'DM')]
       
   for a,b,c in v:
      f = open(filename+"."+str(a)+".stats","w")
      stats = assocCoeff(hmmresults,'mgrm','call',str(c))
      f.write("Cutoff\tTn\tFn\tTp\tFp\tPhi\tTPP\tFPP\tAccuracy\tSpecificity\tSensitivity\n")
      cutoffs = stats.keys()
      cutoffs.sort()
      for c in cutoffs:
         f.write(str(c)+"\t"+"\t".join([str(stats[c]['tn']),str(stats[c]['fn']),str(stats[c]['tp']),str(stats[c]['fp']),str(stats[c]['phi']),str(stats[c]['tpp']),str(stats[c]['fpp']),str(stats[c]['acc']),str(stats[c]['spec']),str(stats[c]['sens'])])+"\n")
      f.close()
   print "Finished analyzing HMM data\n"


if __name__ == '__main__':
	main()

