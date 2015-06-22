#!/usr/bin/env python
# encoding: utf-8
"""
convert2fasta.py.py

Created by Mark Evans on 2010-10-25.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os


def main():
   infile = raw_input("Name of text file to convert: ")
   outfile = infile+".seq"
   f = open(infile,'r')
   o = open(outfile,'w')
   
   for x in f:
      a = x.split()
      print "Len="+str(len(a))+" "+str(a)
      if len(a) == 5:
         o.write(">"+"|".join([a[0],a[1],a[2],a[3]])+"\n"+a[4]+"\n")
      else: o.write(">"+"|".join([a[0],a[1],a[2]," "])+"\n"+a[3]+"\n")
   f.close()
   o.close()
   sys.exit()

if __name__ == '__main__':
   main()

