#!/usr/bin/env python
# encoding: utf-8
"""
editFileName.py
Changes the file type of all files in a directory from one
specified extension to another.  If files do not have
an extension to begin with, all such files will be changed.

Created by Mark Evans on 2011-05-04.
Copyright (c) 2011 __Monogram_Biosciences__. All rights reserved.
"""

import sys
import os,glob,subprocess

def main():
   old = raw_input("What part of filename to remove (e.g.  -RR )?" )
  # new = raw_input("New file extension to use (e.g.  .txt)?" )
   path=""
   os.mkdir('renamed')
   for infile in glob.glob (os.path.join(path,"*")):   # Read all filenames in directory
       if infile.find(".py")==-1: # Take everything except program files
         if infile.find(old) != -1:
            a = infile.split(old)
            idx = len(a)
            if idx != 1:
               new = "".join(a)
               
               print "changing ",infile," to ",new
               os.rename(infile,new)
               os.system("mv "+new+' renamed/'+new)
         

if __name__ == '__main__':
   main()