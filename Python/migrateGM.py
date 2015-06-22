#!/usr/bin/env python
# encoding: utf-8
"""
migrateGM.py
Loop through all legacy Genemachine folders to process data
and place in database

Created by Mark Evans on 2012-01-13.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""

import sys,os,glob,getopt,subprocess, random
from multiprocessing import Process, Manager, Queue
import multiprocessing
from operator import itemgetter


#####################################################################################
# glbl                                                                              #
# Creates a class containing global variables that can be accessed across processes # 
#####################################################################################
class glbl:
   jobs =[]
   data_dirs = {}


def processFiles(data_dir):
    root = os.getcwd()
    data_dirs = []
    if root != data_dir: working_path = os.path.join(root,data_dir)
    else: working_path = root

    LOG = open('bad_data.txt','w')
    data_dirs = os.listdir(working_path)							# get project level directory names
    if data_dirs[0] == '.DS_Store': data_dirs = data_dirs[1:]		# if on a Mac, .DS_Store will be there
    for dirname in data_dirs:										# Loop through project directories
        if os.path.isdir(os.path.join(working_path,dirname)):		# verify that it IS a directory
            next_dir = os.path.join(working_path,dirname)			
            filelist = os.listdir(next_dir)							# get list of files and dirs in the project dir
            if 'Report'+dirname+'.csv' in filelist and 'seqs' in filelist:	# make sure that project is 'complete', has .csv and seqs dir
                csv = open(os.path.join(next_dir,'Report'+dirname+'.csv'))	# open csv file for reading
                sample_data = {}
                files = {}
                suffix = ['TXT']
                print "opening CSV file....\n"
                for line in csv:
                    line = line.rstrip()
                    if line[-1]==',': line = line[:-1]
                    a = line.split(',')
                    if a[0] != 'Sample Name':
                        if len(a) == 24:
                            if a[0].find('.') != -1 and a[0][-3:] not in suffix: 
                                print "##### Didn't see ending"
                                suffix.append(a[0][-3:])
                                print suffix
                            sample_data[a[0].upper().replace('-','_')] = {'Genbank':a[1],
                                             'PNGS<br>Total':a[2],
                                             'V1<br>length':a[3],
                                             'V2<br>length':a[4],
                                             'V4<br>length':a[5],
                                             'V5<br>length':a[6],
                                             'HMM Aln':a[7],
                                             'V3 Sequence':a[8],
                                             'Tropism<br>SVM<br>Prediction':a[9],
                                             'Tropism<br>PSSM<br>Prediction':a[10],
                                             'Tropism (11/25)<br>Prediction':a[11],
                                             'Tropism Pillai<br>Prediction':a[12],
                                             'Length (aa)':a[13],
                                             'PNGS V1':a[14],
                                             'PNGS V2':a[15],
                                             'PNGS V3':a[16],
                                             'PNGS V4':a[17],
                                             'PNGS V5':a[18],
                                             'V1-V4<br>Length':a[19],
                                             'V1-V4<br>PNGS':a[20],
                                             'V1-V5<br>Length':a[21],
                                             'V1-V5<br>PNGS':a[22],
                                             'SVM Pred':a[23]}
                        elif len(a) == 25:
                            sample_data[a[0].upper().replace('-','_')] = {'Genbank':a[1],
                                             'PNGS<br>Total':a[2],
                                             'V1<br>length':a[3],
                                             'V2<br>length':a[4],
                                             'V4<br>length':a[5],
                                             'V5<br>length':a[6],
                                             'HMM Aln':a[7],
                                             'V3 Sequence':a[8],
                                             'Tropism<br>SVM<br>Prediction':a[9],
                                             'Tropism<br>PSSM<br>Prediction':a[10],
                                             'Tropism (11/25)<br>Prediction':a[11],
                                             'Tropism Pillai<br>Prediction':a[12],
                                             'Length (aa)':a[13],
                                             'PNGS V1':a[14],
                                             'PNGS V2':a[15],
                                             'PNGS V3':a[16],
                                             'PNGS V4':a[17],
                                             'PNGS V5':a[18],
                                             'V1-V4<br>Length':a[19],
                                             'V1-V4<br>PNGS':a[20],
                                             'V1-V5<br>Length':a[21],
                                             'V1-V5<br>PNGS':a[22],
                                             'SVM Pred':a[23],
                                             'PSSM<br>Score':a[24]}
                        else: 
                            print "###"+dirname+"'s CSV file did not have the right number of columns ###" 
                            LOG.write("###"+dirname+"'s CSV file did not have the right number of columns ###\n")
  #              print "PATH= ",os.path.join(next_dir,'seqs')
   #             print "length of sample_data = ",str(len(sample_data))
                if len(sample_data) != 0:

                    filelist = os.listdir(os.path.join(next_dir,'seqs'))
 #                   print filelist
                    print sample_data.keys()
                    print suffix
                    for fn in filelist:
                        if fn[-3:] in suffix:
                            print fn.upper()
                            if sample_data.has_key(fn.upper()): files[fn.upper()] = os.path.join(next_dir,'seqs/'+fn)
                            elif sample_data.has_key(fn[:-4].upper()): files[fn[:-4].upper()] = os.path.join(next_dir,'seqs/'+fn)
                        elif sample_data.has_key(fn[:-4].upper().replace('-','_')): files[fn[:-4].upper().replace('-','_')] = os.path.join(next_dir,'seqs/'+fn)
                        
                print files.keys()
                print dirname+" has "+str(len(sample_data))," samples and ",str(len(files))," files"
            else: 
                print dirname, " did not contain all required files"
                LOG.write(dirname+" did not contain all required files\n")
                for path,dirs,filez in os.walk(os.path.join(working_path,dirname), topdown=False):
                    print "Removing empty dir ",dirname
                    for name in filez: 
                        os.remove(os.path.join(path,name))
                    for name in dirs: 
                        os.rmdir(os.path.join(path,name))
                    os.rmdir(os.path.join(working_path,dirname))


    			








    #os.path.join(os.path.join(root,data_dir),os.listdir(data_dir)[5])  # <= this is an alternative method, would give more control

# probably not going to use this method this time
 #   for path, dirs, files in os.walk(working_path):
 #   	if len(files) != 0 and len(dirs)==1:			# should be in a project directory with a seq dir
 #   		if files[0].find('.csv') != -1:				# make sure that there is a complete report, i.e. .csv file
 #   			csv = open(os.path.join(path,a[0])) 	# open csv file to read contents


  #      if len(dirs) == 0: glbl.data_dirs[path] = ''
    

    # Multiprocessing Section
    #########################################
#    Qids = glbl.data_dirs.keys()
#    manager = Manager()                                      # creates shared memory manager object
#    results = manager.dict()                                 # Add dictionary to manager, so it can be accessed across processes
#    nextid = Queue()                                         # Create Queue object to serve as shared id generator across processes
#    for qid in Qids: nextid.put(qid)                         # Load the ids to be tested into the Queue
#    for x in range(0,multiprocessing.cpu_count()):           # Create one process per logical CPU
#        p = Process(target=processData, args=(nextid,results)) # Assign process to processCBR function, passing in the Queue and shared dictionary
#        glbl.jobs.append(p)                                   # Add the process to a list of running processes
#        p.start()                                             # Start process running
#    for j in glbl.jobs:
#        j.join()                                              # For each process, join them back to main, blocking on each one until finished
    
    # write out results
#    c = 1
#    sets = results.keys()
#    sets.sort()
#    for x in sets:
#        FINAL = open('result'+str(c)+'.txt','w')
#        n = "\n************************************************************************************************\n"
#        FINAL.write(n+"* "+x+'    *\n'+n+results[x]+"\n")
#        FINAL.close()     
#        c += 1
            

def main(argv):
    data_dir=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"d:",["data_dir"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-d": 
            data_dir = arg
            processFiles(data_dir)
        else:
            print "There was an error, not right args\n"
            sys.exit()
       
   

if __name__ == '__main__':
    main(sys.argv[1:])