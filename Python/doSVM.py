#!/usr/bin/env python
# encoding: utf-8
"""
doSVM.py

Created by Mark Evans on 2011-02-04.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from libsvm.svm import *
from libsvm.svmutil import *


def main():
   # Make 0/1 vector
   vect = {}
   mut = {"nt_r":0, "nt_y":3, "nt_w":2, "nt_k":1, "aa_x":6}
#   seq = list('CTRPNNNTRKGIRIGPGRAVYTAEKIIGNIRKAHC')
   seq = list('CTRPNNNTRKGIHIGPGRAFYATGEIIGDIRQAHC') # X4
   ids = {'-':(0,1),'A':(0,1), 'C':(0,1), 'E':(0,1), 'D':(0,1), 'G':(0,1), 'F':(0,1), 'I':(0,1), 'H':(0,1), 'K':(0,1), 'M':(0,1), 'L':(0,1), 'N':(0,1), 'Q':(0,1), 'P':(0,1), 'S':(0,1), 'R':(0,1), 'T':(0,1), 'W':(0,1), 'V':(0,1), 'Y':(0,1),'Z':(0,1)}
   labels = []
   alpha = ids.keys()
   alpha.sort()
   for n in range(0,770):
      for s in seq:
         for a in alpha:
            if s == a: 
               vect[n] = ids[a][1]
            else: 
               vect[n] = ids[a][0]
   vect[770] = (float(mut["nt_r"]))
   vect[771] = (float(mut["nt_y"]))
   vect[772] = (float(mut["nt_w"]))
   vect[773] = (float(mut["nt_k"]))
   vect[774] = (float(mut["aa_x"]))
   
   print vect
   print "\n\nLength of vect is ",len(vect)
   x = [vect]
   y = [0]*len(x)
   # Begin SVM
   m = svm_load_model('python_svm_model.svm')
   p_label, p_acc, p_val = svm_predict(y, x, m)
   print "p_label: ",p_label
   print "p_acc: ",p_acc
   print "p_val: ",p_val
   
if __name__ == '__main__':
	main()

