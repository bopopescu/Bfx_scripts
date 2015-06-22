#!/usr/bin/env python
# encoding: utf-8
"""
changeFileExt.py
Changes the file type of all files in a directory from one
specified extension to another.  If files do not have
an extension to begin with, all such files will be changed.

Created by Mark Evans on 2011-05-04.
Copyright (c) 2011 __Monogram_Biosciences__. All rights reserved.
"""

import sys
import os,glob

def main():
   old = raw_input("What file extension to change (e.g.  .seq)?" )
   new = raw_input("New file extension to use (e.g.  .txt)?" )
   path=""
   for infile in glob.glob (os.path.join(path,"*")):            # Read all filenames in directory
       if infile.find(".py")==-1: # Take everything except program files
         a = infile.split('.')
         idx = len(a)
         if idx != 1:
            old = old.replace('.','')
            if a[idx-1].upper() == old.upper():
               del(a[idx-1])
               print "changing ",infile," to ",''.join(a)+new
               os.rename(infile,''.join(a)+new)
         else:
            print "changing ",infile," to ",infile+new
            os.rename(infile,infile+new)

if __name__ == '__main__':
   main()