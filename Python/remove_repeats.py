#!/usr/bin/env python
# encoding: utf-8
"""
remove_repeats.py

Created by Mark Evans on 2011-06-08.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import hashlib


def main():
   patients={}
   exclude_acc =[]
   keep_acc =[]
   rpt_patients=[]
   y09={}
   y10={}
   f = open("patient_export3.txt","r")
   for x in f:
      a = x.upper().split('\t')
      if a[0][0:2] in ('03','04','05','06','07','08','09','10'):
         ptid = hashlib.md5(a[1]+a[2]+a[3]+a[4]).hexdigest()
         if ptid:
            if patients.has_key(ptid):
               patients[ptid].append(a[0])
            else:
               patients[ptid] = [a[0]]
   f.close()
   print "\nThere are "+str(len(patients))+" unique patients\n"
   
   for ptid in patients:
      if len(patients[ptid]) == 1:
         keep_acc.append(patients[ptid][0])
      else:
         patients[ptid].sort()
         keep_acc.append(patients[ptid][0])
         del(patients[ptid][0])
         for x in patients[ptid]:
            if x not in exclude_acc: 
               exclude_acc.append(x)
               if x[0:2]=='09' and not y09.has_key(x): y09[x]=''
               if x[0:2]=='10' and not y10.has_key(x): y10[x]=''
         rpt_patients.append(ptid)
         
   print "There are "+str(len(rpt_patients))+" patients that have been seen more than once\n"
   print "Resulting in "+str(len(exclude_acc))+" accessions to be excluded\n"
   print "There are "+str(len(y09))+" accessions in 2009\n"
   print "There are "+str(len(y10))+" accessions in 2010\n"
   f2 = open("acc2exclude_total.txt","w")
   f3 = open("acc2exclude_y0910.txt","w")
   for val in exclude_acc:
      f2.write(val+"\n")
   f2.close()
   for val in y09.keys():
      f3.write(val+"\n")
   for val in y10.keys():
      f3.write(val+"\n")
   f3.close()
   print "Finished...\n"

if __name__ == '__main__':
	main()

