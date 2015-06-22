#!/usr/bin/env python
# encoding: utf-8
"""
format4svm.py

Created by Mark Evans on 2011-02-07.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os


def main():
   print "\n\nThis script reformats a data file from the R format of 775-length vector of 0/1"
   print "to a python format 775-length dictionary"
   filename = raw_input("What is the name of file in R format to convert: ")
   f = open(filename,"r")
   f2 = open(filename+".pythonsvm","w")
   for line in f:
      a = line.split()
      f2.write(a[len(a)-1]+" ")
      for x in range(0,775):
         f2.write(str(x)+":"+a[x]+" ")
      f2.write("\n")
   f.close()
   f2.close()
   


if __name__ == '__main__':
	main()

