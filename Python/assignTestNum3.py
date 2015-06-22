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

# Returns current time
#########################
def gt():
   return strftime("%d %b %Y %H:%M:%S",localtime())



def main():
   filename = raw_input("What is name of patient_test_export file? ")
   f1 = open(filename,'r')
   f2 = open(filename+'.processed.txt','w')
   print gt()+"\tstarted... \n"
   
   c = 1
   cid = ''
   
   for line in f1:
      a = line.split('\t')
      a[-1].rstrip()  # if did not export with trailing null column, use this
 #     del a[-1]  # if exported with last column as null, use this
      if a[1] != cid:
         cid = a[1]
         c = 1
         f2.write(a[0]+"\t"+str(c)+"\n")
      elif a[1]==cid:
         c = c + 1
         f2.write(a[0]+"\t"+str(c)+"\n")
   f2.close()
   f1.close()
   print gt()+"\tfinished....\n\n"

if __name__ == '__main__':
   main()

