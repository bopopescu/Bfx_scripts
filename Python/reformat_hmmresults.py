#!/usr/bin/env python
# encoding: utf-8
"""
reformat_hmmresults.py

Created by Mark Evans on 2010-10-13.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os,math
from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np


def main():
   filename = raw_input("Enter the name of the hmm result file (.hmmresult)? ")
   f = open(filename,'r')
   outfile = open(filename+".summ","w")
   
   # tnscore = bothR5.TN, r5score = bothR5.all, toscore = toro2R5,tnxscore = bothDX.TN,dxscore= bothDX,toxscore = toro2DX
#   tnscore = {}
   r5score = {}
   toscore = {}
#   tnxscore = {}
   dxscore ={}
   toxscore ={}
   # load raw scores for each HMM
   for line in f:
      a = line.split()
#      if a[0] == 'bothR5.TN': tnscore[a[2]] = a[5]
      if a[0] == 'bothR5.all': r5score[a[2]] = a[5]
      elif a[0] == 'toro2R5': toscore[a[2]] = a[5]
      elif a[0] == 'toro2DX': toxscore[a[2]] = a[5]
      elif a[0] == 'bothDX': dxscore[a[2]] = a[5]
#      elif a[0] == 'bothDX.TN': tnxscore[a[2]] = a[5]
      
   keys = r5score.keys()
   keys.sort()
   
   ##############################################################################
   # Do the summary without using cutoff, just comparing R5 HMM vs DX HMM score #
   ##############################################################################
#   R1stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"r5tn_vs_dxtn"}
#   R2stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"r5tn_vs_dx"}
#   R3stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"r5tn_vs_tdx"}
#   R4stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"r5_vs_dxtn"}
   R5stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"r5_vs_dx"}
   R6stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"r5_vs_tdx"}
#   R7stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"tr5_vs_dxtn"}
   R8stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"tr5_vs_dx"}
   R9stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"tr5_vs_tdx"}   
#   stats = [R1stats,R2stats,R3stats,R4stats,R5stats,R6stats,R7stats,R8stats,R9stats]
   stats = [R5stats,R6stats,R8stats,R9stats]
   predict={} # hold the vector of all the various calls by sequence if we want them again
   for k in keys:
      ids = k.split('|')
#      r1 = "R5" # r5tn_vs_dxtn
#      r2 = "R5" # r5tn_vs_dx
#      r3 = "R5" # r5tn_vs_tdx
#      r4 = "R5" # r5_vs_dxtn
      r5 = "R5" # r5_vs_dx
      r6 = "R5" # r5_vs_tdx
#      r7 = "R5" # tr5_vs_dxtn
      r8 = "R5" # tr5_vs_dx
      r9 = "R5" # tr5_vs_tdx
      mgrmcall="R5"
      if ids[2]=='DM' or ids[2]=='X4' or ids[3]=='DM' or ids[3]=='X4': mgrmcall = 'DM'
#      if tnscore[k] <= tnxscore[k]: r1 = "DM"
#      if tnscore[k] <= dxscore[k] : r2 = "DM"
#      if tnscore[k] <= toxscore[k]: r3 = "DM"
#      if r5score[k] <= tnxscore[k]: r4 = "DM"
      if r5score[k] <= dxscore[k] : r5 = "DM"
      if r5score[k] <= toxscore[k]: r6 = "DM"
#      if toscore[k] <= tnxscore[k]: r7 = "DM"
      if toscore[k] <= dxscore[k] : r8 = "DM"
      if toscore[k] <= toxscore[k]: r9 = "DM"
#      predict[k] = [mgrmcall,r1,r2,r3,r4,r5,r6,r7,r8,r9]
      predict[k] = [mgrmcall,r5,r6,r8,r9]
   
   for k in keys:
      y = 1
      for x in stats:
         if   predict[k][0]==predict[k][y] and predict[k][0]=='DM':  x['tp'] += 1
         elif predict[k][0]==predict[k][y] and predict[k][0]=='R5':  x['tn'] += 1
         elif predict[k][0]!=predict[k][y] and predict[k][0]=='DM':  x['fn'] += 1
         elif predict[k][0]!=predict[k][y] and predict[k][0]=='R5':  x['fp'] += 1
         y += 1
   
   outfile.write("\n\n==========================\nCombo\tTP\tTN\tFP\tFN\tPHI\tTPP\tFPP\tACC\tSpec\tSens\n")
   for x in stats:
      phi=0.0
      tp = float(x['tp'])
      tn = float(x['tn'])
      fp = float(x['fp'])
      fn = float(x['fn'])
      try:
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
      outfile.write(x['name']+"\t"+str(x['tp'])+"\t"+str(x['tn'])+"\t"+str(x['fp'])+"\t"+str(x['fn'])+"\t"+str(phi)+"\t"+str(tpp)+"\t"+str(fpp)+"\t"+str(acc)+"\t"+str(spec)+"\t"+str(sens)+"\n")
   outfile.write("=========================\n\n")
   
   ##################################################
   # Now do summary per HMM, using a sliding cutoff #
   ##################################################
#   r1stats = {}  # hold all the stats for each cutoff value
   r2stats = {}
   r3stats = {}
#   r4stats = {}
   r5stats = {}
   r6stats = {}
#   r7stats = {}
#   r8stats = {}
#   r9stats = {}
   for i in np.arange(50,100,0.1):
#      r1stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      r2stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      r3stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
#      r4stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      r5stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      r6stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
#      r7stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
#      r8stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
#      r9stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
#  stats2 = [r1stats,r2stats,r3stats]
#   stats3 = [r4stats,r5stats,r6stats]
   stats2 = [r2stats,r3stats]
   stats3 = [r5stats,r6stats]
#   scores = [tnscore,r5score,toscore]
#   scores2 =[tnxscore,dxscore,toxscore]
   scores = [r5score,toscore]
   scores2 =[dxscore,toxscore]
   # Generate test data for R5 HMM's
   for i in np.arange(50,100,0.1):  # The cutoff being tested
      y = 0
      for x in stats2:
         for k in keys:
            if   float(scores[y][k])  > i and predict[k][0] =='DM':  x[i]['fn'] += 1
            elif float(scores[y][k])  > i and predict[k][0] == 'R5': x[i]['tn'] += 1
            elif float(scores[y][k]) <= i and predict[k][0] == 'DM': x[i]['tp'] += 1
            elif float(scores[y][k]) <= i and predict[k][0] == 'R5': x[i]['fp'] += 1
         y += 1
   # Generate test data for X4 HMM's
   for i in np.arange(50,100,0.1):  # The cutoff being tested
      y = 0
      for x in stats3:
         for k in keys:
            if   float(scores2[y][k])  > i and predict[k][0] =='DM':  x[i]['tp'] += 1
            elif float(scores2[y][k])  > i and predict[k][0] == 'R5': x[i]['fp'] += 1
            elif float(scores2[y][k]) <= i and predict[k][0] == 'DM': x[i]['fn'] += 1
            elif float(scores2[y][k]) <= i and predict[k][0] == 'R5': x[i]['tn'] += 1
         y += 1
         
#   outfile.write("\n\n\nStats for bothR5.TN\nCutoff\tTn\tFn\tTp\tFp\tPhi\tTPR\tFPR\tAccuracy\tSpecificity\tSensitivity\n")
#   for i in np.arange(50,100,0.1):
#      tp = float(r1stats[i]['tp'])
#      tn = float(r1stats[i]['tn'])
#      fp = float(r1stats[i]['fp'])
#      fn = float(r1stats[i]['fn'])
#      phi = 0.0
#      try:
#         phi = math.sqrt(math.pow(tp*tn - fp*fn,2)/((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
#      except:
#         phi = 0.0
#      try: tpr = tp/(tp+fn)  # true positive rate
#      except: tpr = 0.0
#      try: fpr = fp/(fp+tn)  # false positive rate
#      except: fpr = 0.0
#      try: acc = (tp+tn)/(tp+tn+fp+fn)  # Total accuracy for prediction versus reality
#      except: acc = 0.0
#      try: spec = tn/(tn+fp)  # specificity for X4
#      except: spec = 0.0
#      try: sens = tp/(tp+fn)  # sensitivity for X4
#      except: sens = 0.0
#      outfile.write(str(i)+"\t"+"\t".join([str(tn),str(fn),str(tp),str(fp),str(phi),str(tpr),str(fpr),str(acc),str(spec),str(sens)])+"\n")
   
   outfile.write("\n\n\nStats for bothR5.all\nCutoff\tTn\tFn\tTp\tFp\tPhi\tTPR\tFPR\tAccuracy\tSpecificity\tSensitivity\n")
   for i in np.arange(50,100,0.1):
      tp = float(r2stats[i]['tp'])
      tn = float(r2stats[i]['tn'])
      fp = float(r2stats[i]['fp'])
      fn = float(r2stats[i]['fn'])
      phi = 0.0
      try:
         phi = math.sqrt(math.pow(tp*tn - fp*fn,2)/((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
      except:
         phi = 0.0
      try: tpr = tp/(tp+fn)  # true positive rate
      except: tpr = 0.0
      try: fpr = fp/(fp+tn)  # false positive rate
      except: fpr = 0.0
      try: acc = (tp+tn)/(tp+tn+fp+fn)  # Total accuracy for prediction versus reality
      except: acc = 0.0
      try: spec = tn/(tn+fp)  # specificity for X4
      except: spec = 0.0
      try: sens = tp/(tp+fn)  # sensitivity for X4
      except: sens = 0.0
      outfile.write(str(i)+"\t"+"\t".join([str(tn),str(fn),str(tp),str(fp),str(phi),str(tpr),str(fpr),str(acc),str(spec),str(sens)])+"\n")
         
   outfile.write("\n\n\nStats for toro2R5\nCutoff\tTn\tFn\tTp\tFp\tPhi\tTPR\tFPR\tAccuracy\tSpecificity\tSensitivity\n")
   for i in np.arange(50,100,0.1):
      tp = float(r3stats[i]['tp'])
      tn = float(r3stats[i]['tn'])
      fp = float(r3stats[i]['fp'])
      fn = float(r3stats[i]['fn'])
      phi = 0.0
      try:
         phi = math.sqrt(math.pow(tp*tn - fp*fn,2)/((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
      except:
         phi = 0.0
      try: tpr = tp/(tp+fn)  # true positive rate
      except: tpr = 0.0
      try: fpr = fp/(fp+tn)  # false positive rate
      except: fpr = 0.0
      try: acc = (tp+tn)/(tp+tn+fp+fn)  # Total accuracy for prediction versus reality
      except: acc = 0.0
      try: spec = tn/(tn+fp)  # specificity for X4
      except: spec = 0.0
      try: sens = tp/(tp+fn)  # sensitivity for X4
      except: sens = 0.0
      outfile.write(str(i)+"\t"+"\t".join([str(tn),str(fn),str(tp),str(fp),str(phi),str(tpr),str(fpr),str(acc),str(spec),str(sens)])+"\n")

#   outfile.write("\n\n\nStats for bothDX.TN\nCutoff\tTn\tFn\tTp\tFp\tPhi\tTPR\tFPR\tAccuracy\tSpecificity\tSensitivity\n")
#   for i in np.arange(50,100,0.1):
#      tp = float(r4stats[i]['tp'])
#      tn = float(r4stats[i]['tn'])
#      fp = float(r4stats[i]['fp'])
#      fn = float(r4stats[i]['fn'])
#      phi = 0.0
#      try:
#         phi = math.sqrt(math.pow(tp*tn - fp*fn,2)/((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
#      except:
#         phi = 0.0
#      try: tpr = tp/(tp+fn)  # true positive rate
#      except: tpr = 0.0
#      try: fpr = fp/(fp+tn)  # false positive rate
#      except: fpr = 0.0
#      try: acc = (tp+tn)/(tp+tn+fp+fn)  # Total accuracy for prediction versus reality
#      except: acc = 0.0
#      try: spec = tn/(tn+fp)  # specificity for X4
#      except: spec = 0.0
#      try: sens = tp/(tp+fn)  # sensitivity for X4
#      except: sens = 0.0
#      outfile.write(str(i)+"\t"+"\t".join([str(tn),str(fn),str(tp),str(fp),str(phi),str(tpr),str(fpr),str(acc),str(spec),str(sens)])+"\n")

   outfile.write("\n\n\nStats for bothDX\nCutoff\tTn\tFn\tTp\tFp\tPhi\tTPR\tFPR\tAccuracy\tSpecificity\tSensitivity\n")
   for i in np.arange(50,100,0.1):
      tp = float(r5stats[i]['tp'])
      tn = float(r5stats[i]['tn'])
      fp = float(r5stats[i]['fp'])
      fn = float(r5stats[i]['fn'])
      phi = 0.0
      try:
         phi = math.sqrt(math.pow(tp*tn - fp*fn,2)/((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
      except:
         phi = 0.0
      try: tpr = tp/(tp+fn)  # true positive rate
      except: tpr = 0.0
      try: fpr = fp/(fp+tn)  # false positive rate
      except: fpr = 0.0
      try: acc = (tp+tn)/(tp+tn+fp+fn)  # Total accuracy for prediction versus reality
      except: acc = 0.0
      try: spec = tn/(tn+fp)  # specificity for X4
      except: spec = 0.0
      try: sens = tp/(tp+fn)  # sensitivity for X4
      except: sens = 0.0
      outfile.write(str(i)+"\t"+"\t".join([str(tn),str(fn),str(tp),str(fp),str(phi),str(tpr),str(fpr),str(acc),str(spec),str(sens)])+"\n")
   
   outfile.write("\n\n\nStats for toro2DX\nCutoff\tTn\tFn\tTp\tFp\tPhi\tTPR\tFPR\tAccuracy\tSpecificity\tSensitivity\n")
   for i in np.arange(50,100,0.1):
      tp = float(r6stats[i]['tp'])
      tn = float(r6stats[i]['tn'])
      fp = float(r6stats[i]['fp'])
      fn = float(r6stats[i]['fn'])
      phi = 0.0
      try:
         phi = math.sqrt(math.pow(tp*tn - fp*fn,2)/((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
      except:
         phi = 0.0
      try: tpr = tp/(tp+fn)  # true positive rate
      except: tpr = 0.0
      try: fpr = fp/(fp+tn)  # false positive rate
      except: fpr = 0.0
      try: acc = (tp+tn)/(tp+tn+fp+fn)  # Total accuracy for prediction versus reality
      except: acc = 0.0
      try: spec = tn/(tn+fp)  # specificity for X4
      except: spec = 0.0
      try: sens = tp/(tp+fn)  # sensitivity for X4
      except: sens = 0.0
      outfile.write(str(i)+"\t"+"\t".join([str(tn),str(fn),str(tp),str(fp),str(phi),str(tpr),str(fpr),str(acc),str(spec),str(sens)])+"\n")

   
   
   f.close()
   outfile.close()
if __name__ == '__main__':
   main()

