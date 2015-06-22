#!/usr/bin/env python
# encoding: utf-8
"""
assignTestNum.py

Created by Mark Evans on 2011-05-11.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from time import localtime, strftime

def main():
   print "\nstarted at ",strftime("%H:%M:%S", localtime()),"\n\n"
   f1 = open('patient_test_foo.txt','r')
   f2 = open('patient_test.txt','w')
   c = 1
   cid = ''
   
   for line in f1:
      a = line.split()
      a[-1].rstrip()
      if a[0] != cid:
         cid = a[0]
         c = 1
         f2.write('\t'.join(a)+"\t1\n")
      elif a[0]==cid:
         c = c + 1
         f2.write('\t'.join(a)+"\t"+str(c)+"\n")
   f2.close()
   f1.close()
   print "\nfinished at ",strftime("%H:%M:%S", localtime()),"\n\n"

if __name__ == '__main__':
   main()

