#!/usr/bin/env python
# encoding: utf-8
"""
processPatches.py
Loop through all patch files and evaluate significance of each one

Created by Mark Evans on 2012-01-11.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""

import sys,os,glob,getopt,subprocess, random
import os.path
from multiprocessing import Process, Manager, Queue
import multiprocessing
from operator import itemgetter
from rpy2.robjects.vectors import DataFrame
from rpy2.robjects.packages import importr
from rpy2 import robjects

#####################################################################################
# glbl                                                                              #
# Creates a class containing global variables that can be accessed across processes # 
#####################################################################################
class glbl:
   jobs =[]
   data_dirs = {}

def analyzeR():
    mypatchData = DataFrame.from_csvfile(inputfile,header=True, sep = "\t")
    
    ici = list(mypatchData.colnames).index('ic50_val')
    starti = ici + 1
    lasti = mypatchData.nrow
    myData = mypatchData[,ici:lasti]


def processData(nid,results):
    path =''
    while 1:
        try: path = nid.get(False)                                                             # Get next query id from Queue generator, 
        except: 
            print " CPU not in use"
            break
        print path

        tmp = str(random.random())[2:20]  # generate random temp filename
        tmp2 = str(random.random())[2:20]
        f = open(tmp,'w')
        f.write(rcode)
        f.close()
        f2 = open(tmp2,'w')
        try:
            process = subprocess.Popen(["R","--vanilla","--subordinate","-f",tmp], stdout=f2)
            process.wait()
        except:
            print "R subprocess failed"
            os.remove(tmp)
            f2.close()
            os.remove(tmp2)
      #      multiprocessing.current_process.terminate()
    
        os.remove(tmp)
        f2.close()
        if os.path.getsize(tmp2) != 0:
          f = open(tmp2,'r')
          results[path] = f.read()
          f.close()
          os.remove(tmp2)
        else:
          results[path] = 'None'
          os.remove(tmp2)
    return




def processFiles(patch_dir):
    root = os.getcwd()
    glbl.data_dirs = {}
    if root != patch_dir: working_path = root+"/"+patch_dir
    else: working_path = root

    for path, dirs, files in os.walk(working_path):
        if len(dirs) == 0: glbl.data_dirs[path] = ''
    

    # Multiprocessing Section
    #########################################
    Qids = glbl.data_dirs.keys()
    manager = Manager()                                      # creates shared memory manager object
    results = manager.dict()                                 # Add dictionary to manager, so it can be accessed across processes
    nextid = Queue()                                         # Create Queue object to serve as shared id generator across processes
    for qid in Qids: nextid.put(qid)                         # Load the ids to be tested into the Queue
    for x in range(0,multiprocessing.cpu_count()):           # Create one process per logical CPU
        p = Process(target=processData, args=(nextid,results)) # Assign process to processCBR function, passing in the Queue and shared dictionary
        glbl.jobs.append(p)                                   # Add the process to a list of running processes
        p.start()                                             # Start process running
    for j in glbl.jobs:
        j.join()                                              # For each process, join them back to main, blocking on each one until finished
    
    # write out results
    c = 1
    sets = results.keys()
    sets.sort()
    for x in sets:
        if results[x] != 'None':
            FINAL = open('result'+str(c)+'.txt','w')
            n = "\n************************************************************************************************\n"
            FINAL.write(n+"* "+x+'    *\n'+n+results[x]+"\n")
            FINAL.close()     
            c += 1
            

def main(argv):
    rlib = importr('e1071')
    base = importr('base')
    r = robjects.r
    patch_dir=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"d:",["patch_dir"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-d": 
            patch_dir = arg
            processFiles(patch_dir)
        else:
            print "There was an error, not right args\n"
            sys.exit()
       
   

if __name__ == '__main__':
    main(sys.argv[1:])