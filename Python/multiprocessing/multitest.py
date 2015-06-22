#!/usr/bin/env python
# encoding: utf-8
"""
multitest.py

Created by Mark on 2010-11-07.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import multiprocessing
from multiprocessing import Process, Pool, Manager
import time


def worker(d, c, s,val):
   name = multiprocessing.current_process().name
   print "starting process ",name
   d[val] = {"s value=":s[val],"c value=":c['foo1']}


def main():
   if __name__=="__main__":
      # Establish communication queues
      tasks = multiprocessing.JoinableQueue()
      results = multiprocessing.Queue()
      
      mgr = multiprocessing.Manager()
      d = mgr.dict()
      s = mgr.dict()
      c = mgr.dict()
      s['seq1'] = {'a':1,'b':2,'c':3}
      s['seq2'] = {'d':1,'e':2,'f':3}
      c['foo1'] = {'m':1,'b':2,'f':3}
      c['foo2'] = {'h':1,'e':2,'c':3}
      k = s.keys()
      jobs = []
      for val in k:
         p = multiprocessing.Process(target=worker, args=(d, c, s,val))
         p.start()
         jobs.append(p)
      for j in jobs: j.join()
      print 'Results:', d
      
	      


if __name__ == '__main__':
	main()

