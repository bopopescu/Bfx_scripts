#!/usr/bin/env python
# encoding: utf-8
"""
reformat_hmmresults.py

Created by Mark Evans on 2010-10-26.
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
from rpy2.robjects.packages import importr
from rpy2 import robjects
from rpy2.robjects.lib import grid
from rpy2.robjects.vectors import DataFrame

# create basic R functions
rprint = robjects.globalenv.get("print")
stats = importr('stats')
grdevices = importr('grDevices')
base = importr('base')
palette = grdevices.palette()
ro = robjects
r = robjects.r
plot = robjects.r.plot

def main():
   filename = raw_input("Enter the name of the hmm result file (.hmmresult)? ")
   f = open(filename,'r')
   outfile = open(filename+".summ","w")
   
   # tnscore = bothR5.TN, r5score = bothR5.all, toscore = toro2R5,tnxscore = bothDX.TN,dxscore= bothDX,toxscore = toro2DX
   r5score = {}
   toscore = {}
   dxscore ={}
   toxscore ={}
   # load raw scores for each HMM
   for line in f:
      a = line.split()
      if a[0] == 'bothR5.all': r5score[a[2]] = a[5]
      elif a[0] == 'toro2R5': toscore[a[2]] = a[5]
      elif a[0] == 'toro2DX': toxscore[a[2]] = a[5]
      elif a[0] == 'bothDX': dxscore[a[2]] = a[5]
      
   keys = r5score.keys()
   keys.sort()
   
   ##############################################################################
   # Do the summary without using cutoff, just comparing R5 HMM vs DX HMM score #
   ##############################################################################

   R5stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"r5_vs_dx"}
   R6stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"r5_vs_tdx"}

   R8stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"tr5_vs_dx"}
   R9stats = {'tp':0,'tn':0,'fp':0,'fn':0,'name':"tr5_vs_tdx"}   

   stats = [R5stats,R6stats,R8stats,R9stats]
   predict={} # hold the vector of all the various calls by sequence if we want them again
   for k in keys:
      ids = k.split('|')
      r5 = "R5" # r5_vs_dx
      r6 = "R5" # r5_vs_tdx
      r8 = "R5" # tr5_vs_dx
      r9 = "R5" # tr5_vs_tdx
      mgrmcall="R5"
      if ids[2]=='DM' or ids[2]=='X4' or ids[3]=='DM' or ids[3]=='X4': mgrmcall = 'DM'
      if r5score[k] <= dxscore[k] : r5 = "DM"
      if r5score[k] <= toxscore[k]: r6 = "DM"
      if toscore[k] <= dxscore[k] : r8 = "DM"
      if toscore[k] <= toxscore[k]: r9 = "DM"
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
  # hold all the stats for each cutoff value
   r2stats = {}
   r3stats = {}
   r5stats = {}
   r6stats = {}
   for i in np.arange(50,100,0.1):
      r2stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      r3stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      r5stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      r6stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      
   stats2 = [r2stats,r3stats]
   stats3 = [r5stats,r6stats]
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
    
   dicts = {'TnR5':r2stats,'TxR5':r3stats,'TnDX':r5stats,'TxDX':r6stats}      
   
   rdict = {}
   for hmm in dicts:
      outfile.write("\n\n\nStats for "+hmm+"\nCutoff\tTn\tFn\tTp\tFp\tPhi\tTPR\tFPR\tAccuracy\tSpecificity\tSensitivity\n")
      st = {'tp':[],'tn':[],'fp':[],'fn':[],'phi':[],'tpr':[],'fpr':[],'acc':[],'spec':[],'sens':[],'cut':[]}
      for i in np.arange(50,100,0.1):
         st['cut'].append(i)
         tp = float(dicts[hmm][i]['tp'])
         tn = float(dicts[hmm][i]['tn'])
         fp = float(dicts[hmm][i]['fp'])
         fn = float(dicts[hmm][i]['fn'])
         st['tp'].append(tp)
         st['tn'].append(tn)
         st['fp'].append(fp)
         st['fn'].append(fn)
         phi = 0.0
         try:
            phi = math.sqrt(math.pow(tp*tn - fp*fn,2)/((tp+fp)*(fn+tn)*(tp+fn)*(fp+tn)))
         except:
            phi = 0.0
         st['phi'].append(phi)
         try: tpr = tp/(tp+fn)  # true positive rate
         except: tpr = 0.0
         st['tpr'].append(tpr)
         try: fpr = fp/(fp+tn)  # false positive rate
         except: fpr = 0.0
         st['fpr'].append(fpr)
         try: acc = (tp+tn)/(tp+tn+fp+fn)  # Total accuracy for prediction versus reality
         except: acc = 0.0
         st['acc'].append(acc)
         try: spec = tn/(tn+fp)  # specificity for X4
         except: spec = 0.0
         st['spec'].append(spec)
         try: sens = tp/(tp+fn)  # sensitivity for X4
         except: sens = 0.0
         st['sens'].append(sens)
         outfile.write(str(i)+"\t"+"\t".join([str(tn),str(fn),str(tp),str(fp),str(phi),str(tpr),str(fpr),str(acc),str(spec),str(sens)])+"\n")
      rdict[hmm] = ro.DataFrame({'CUT':ro.FloatVector(st['cut']),'TP':ro.FloatVector(st['tp']),'FP':ro.FloatVector(st['fp']),'TN':ro.FloatVector(st['tn']),'FN':ro.FloatVector(st['fn']),'PHI':ro.FloatVector(st['phi']),'TPR':ro.FloatVector(st['tpr']),'FPR':ro.FloatVector(st['fpr']),'ACC':ro.FloatVector(st['acc']),'SPEC':ro.FloatVector(st['spec']),'SENS':ro.FloatVector(st['sens'])})
   
   f.close()
   outfile.close()
   
   
   #########################
   # Generate pdf graphs
   #########################
   hm = ['TnR5','TnDX','TxR5','TxDX']
   
   # Begin witing graphs to pdf
   grdevices.pdf(file=filename+".stats.pdf",width=7,height=7)
   
   # ROC curve for TnR5/DX and TxR5/DX
   plot(rdict['TnR5'].rx2("FPR"),rdict['TnR5'].rx2("TPR"),type='l',main='ROC of all MGRM non-35aa seq',xlab='FPR',ylab='TPR', col='blue')
   r.lines(rdict['TxR5'].rx2("FPR"),rdict['TxR5'].rx2("TPR"),type='l',col='blue',lty=2)
   r.lines(rdict['TnDX'].rx2("FPR"),rdict['TnDX'].rx2("TPR"),type='l',col='red',lty=1)
   r.lines(rdict['TxDX'].rx2("FPR"),rdict['TxDX'].rx2("TPR"),type='l',col='red',lty=2)
   r.legend('bottomright',base.c('TnR5','TxR5','TnDX','TxDX'),col=base.c('blue','blue','red','red'),lty=base.c(1,2,1,2),bg='white')
   
   # Sensitivity/Specificity plot for R5 HMMs
   plot(rdict['TnR5'].rx2("CUT"),rdict['TnR5'].rx2("SPEC"),type='l',main='R5 Sensitivity & Specificity vs Cutoff for all MGRM non-35aa seq',xlab='Cutoff Score',ylab='', col='blue',xaxs='i',tck=1,ylim=base.c(0,1),xlim=base.c(40,100),xaxp=base.c(40,100,30))
   r.lines(rdict['TnR5'].rx2("CUT"),rdict['TxR5'].rx2("SPEC"),type='l',col='blue',lty=2)
   r.lines(rdict['TnR5'].rx2("CUT"),rdict['TnR5'].rx2("SENS"),type='l',col='red',lty=1)
   r.lines(rdict['TnR5'].rx2("CUT"),rdict['TxR5'].rx2("SENS"),type='l',col='red',lty=2)
   r.legend('right',base.c('TnR5 Specificity','TxR5 Specificity','TnR5 Sensitivity','TxR5 Sensitivity'),col=base.c('blue','blue','red','red'),lty=base.c(1,2,1,2),bg='white')
   
   # Sensitivity/Specificity plot for X4 HMMs
   plot(rdict['TnDX'].rx2("CUT"),rdict['TnDX'].rx2("SPEC"),type='l',main='DX Sensitivity & Specificity vs Cutoff for all MGRM non-35aa seq',xlab='Cutoff Score',ylab='', col='blue',xaxs='i',tck=1,ylim=base.c(0,1),xlim=base.c(40,100),xaxp=base.c(40,100,30))
   r.lines(rdict['TnDX'].rx2("CUT"),rdict['TxDX'].rx2("SPEC"),type='l',col='blue',lty=2)
   r.lines(rdict['TnDX'].rx2("CUT"),rdict['TnDX'].rx2("SENS"),type='l',col='red',lty=1)
   r.lines(rdict['TnDX'].rx2("CUT"),rdict['TxDX'].rx2("SENS"),type='l',col='red',lty=2)
   r.legend('right',base.c('TnDX Specificity','TxDX Specificity','TnDX Sensitivity','TxDX Sensitivity'),col=base.c('blue','blue','red','red'),lty=base.c(1,2,1,2),bg='white')
   
   
   # close pdf file
   grdevices.dev_off()
   
if __name__ == '__main__':
   main()

