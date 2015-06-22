#!/usr/bin/env python
# encoding: utf-8
"""
hmm_score_diff.py

Created by Mark Evans on 2010-10-29.
Copyright (c) 2010 __Monogram Biosciences__. All rights reserved.
"""

import sys
import os,math
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
   outfile = open(filename+".R5X4diff","w")
   
   # r5score = bothR5.all, toxscore = toro2DX
   r5score = {}
   toxscore ={}
   # load raw scores for each HMM
   for line in f:
      a = line.split()
      if a[0] == 'bothR5.all': r5score[a[2]] = a[5]
      elif a[0] == 'toro2DX': toxscore[a[2]] = a[5]
      
   keys = r5score.keys()
   keys.sort()
   
   ########################################################################################
   # Do summary per HMM, using a sliding value for the difference between R5 and DX score #
   ########################################################################################
  # hold all the stats for each cutoff value
   r1stats = {}
   r2stats = {}

   for i in np.arange(0,15,0.1):
      r1stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      r2stats[i] = {'tp':0,'tn':0,'fp':0,'fn':0}
      
#   stats = [r1stats,r2stats]
#   scores = [r5score,toxscore]
   # Generate test data for HMM's
   for i in np.arange(0,15,0.1):  # The HMM score difference being tested
      for k in keys:              # sequence id string
         ids = k.split('|')       # get original tropism calls
         c = 'R5'
         if ids[2] != c or ids[3] != c: c = 'DM' # merge STF and ESTA calls into one DM call
         if   (float(r5score[k]) - float(toxscore[k]) > i) and c == 'DM': r1stats[i]['fn'] += 1
         elif (float(r5score[k]) - float(toxscore[k]) > i) and c == 'R5': r1stats[i]['tn'] += 1
         elif (float(r5score[k]) - float(toxscore[k])<= i) and c == 'DM': r1stats[i]['tp'] += 1
         elif (float(r5score[k]) - float(toxscore[k])<= i) and c == 'R5': r1stats[i]['fp'] += 1
   
   
   outfile.write("\n\n\nCutoff\tTn\tFn\tTp\tFp\tPhi\tTPR\tFPR\tAccuracy\tSpecificity\tSensitivity\n")
   st = {'tp':[],'tn':[],'fp':[],'fn':[],'phi':[],'tpr':[],'fpr':[],'acc':[],'spec':[],'sens':[],'cut':[]}
   rdata={}
   for i in np.arange(0,15,0.1):
      st['cut'].append(i)
      tp = float(r1stats[i]['tp'])
      tn = float(r1stats[i]['tn'])
      fp = float(r1stats[i]['fp'])
      fn = float(r1stats[i]['fn'])
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
      rdata['TnR5'] = ro.DataFrame({'CUT':ro.FloatVector(st['cut']),'TP':ro.FloatVector(st['tp']),'FP':ro.FloatVector(st['fp']),'TN':ro.FloatVector(st['tn']),'FN':ro.FloatVector(st['fn']),'PHI':ro.FloatVector(st['phi']),'TPR':ro.FloatVector(st['tpr']),'FPR':ro.FloatVector(st['fpr']),'ACC':ro.FloatVector(st['acc']),'SPEC':ro.FloatVector(st['spec']),'SENS':ro.FloatVector(st['sens'])})
   
   f.close()
   outfile.close()
   
   
   #########################
   # Generate pdf graphs
   #########################
   
   # Begin witing graphs to pdf
   grdevices.pdf(file=filename+".stats.pdf",width=8,height=7)
   
   # ROC curve for TnR5-TxDX
   plot(rdata['TnR5'].rx2("FPR"),rdata['TnR5'].rx2("TPR"),type='l',main='ROC of TnR5-TxDX Cutoff',xlab='FPR',ylab='TPR', col='blue')
#   r.lines(rdict['TxR5'].rx2("FPR"),rdict['TxR5'].rx2("TPR"),type='l',col='blue',lty=2)
#   r.lines(rdict['TnDX'].rx2("FPR"),rdict['TnDX'].rx2("TPR"),type='l',col='red',lty=1)
#   r.lines(rdict['TxDX'].rx2("FPR"),rdict['TxDX'].rx2("TPR"),type='l',col='red',lty=2)
   r.legend('bottomright',base.c('TnR5-TxDX'),col=base.c('blue'),lty=base.c(1),bg='white')
   
   # Sensitivity/Specificity plot for R5 HMMs
   plot(rdata['TnR5'].rx2("CUT"),rdata['TnR5'].rx2("SPEC"),type='l',main='R5-DX Sensitivity & Specificity vs Difference Cutoff for all public seq',xlab='Difference Cutoff Score',ylab='', col='blue',xaxs='i',tck=1,ylim=base.c(0,1),xlim=base.c(0,20),xaxp=base.c(0,20,20))
#   r.lines(rdata['TnR5'].rx2("CUT"),rdata['TxR5'].rx2("SPEC"),type='l',col='blue',lty=2)
   r.lines(rdata['TnR5'].rx2("CUT"),rdata['TnR5'].rx2("SENS"),type='l',col='red',lty=1)
#   r.lines(rdict['TnR5'].rx2("CUT"),rdict['TxR5'].rx2("SENS"),type='l',col='red',lty=2)
   r.legend('right',base.c('TnR5-TxDX Specificity','TxR5-TxDX Specificity'),col=base.c('blue','red'),lty=base.c(1,1),bg='white')
   
   # Sensitivity/Specificity plot for X4 HMMs
#   plot(rdict['TnDX'].rx2("CUT"),rdict['TnDX'].rx2("SPEC"),type='l',main='DX Sensitivity & Specificity vs Cutoff for all MGRM non-35aa seq',xlab='Cutoff Score',ylab='', col='blue',xaxs='i',tck=1,ylim=base.c(0,1),xlim=base.c(40,100),xaxp=base.c(40,100,30))
#   r.lines(rdict['TnDX'].rx2("CUT"),rdict['TxDX'].rx2("SPEC"),type='l',col='blue',lty=2)
#   r.lines(rdict['TnDX'].rx2("CUT"),rdict['TnDX'].rx2("SENS"),type='l',col='red',lty=1)
#   r.lines(rdict['TnDX'].rx2("CUT"),rdict['TxDX'].rx2("SENS"),type='l',col='red',lty=2)
#   r.legend('right',base.c('TnDX Specificity','TxDX Specificity','TnDX Sensitivity','TxDX Sensitivity'),col=base.c('blue','blue','red','red'),lty=base.c(1,2,1,2),bg='white')
   
   
   # close pdf file
   grdevices.dev_off()
   
if __name__ == '__main__':
   main()

