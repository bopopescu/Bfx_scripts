#!/usr/bin/env python
# encoding: utf-8
"""
ns5b_fastq_trim.py
Trim fastq files to just ns5b alignment to eliminate primer sequences

Created by Mark Evans on 2012-10-08.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""

import sys,os,getopt,subprocess,glob,random,time
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from multiprocessing import Process, Manager, Queue
import multiprocessing
from operator import itemgetter
from commands import *  # This should be replaced in Python 2.7+ with subprocess, but we need it for Bali


class glbl:
    jobs =[]
    file_chunks = []


#############
# Usage     #
#############
def usage():
    code =  "\n\n#############\nProgram parses fastq files and trims reads to NS5B seq only \n\n"
    code += "Usage: $ python ns5b_fastq_trim.py [-s] fastq_filename \n\n"
    code += "-s [seq] name of fastq file to process\n"
    code += "-h help\n\n#############\n\n"
    print code
    return



#########################################################
# Calls hmmscan from HMMER3.0 and generates result file #
#########################################################
def doHMMScan(filename,local_path,hmmscan_bin):
    tmp = str(random.random())[2:]
    #hmmscan_bin = "/usr/local/bin/hmmscan"
    err = open(local_path+"hmmscan_error_log"+tmp,"w")
    process = subprocess.Popen([hmmscan_bin,"-o",local_path+"hmmscan.log.tmp"+tmp,"--noali","--domT","20","--domtblout",local_path+filename+".hmmresult",local_path+"ns5b_ends.hmm",local_path+filename],stderr=err)
    process.wait()
    err.close()
    os.remove(local_path+"hmmscan_error_log"+tmp)
    os.remove(local_path+"hmmscan.log.tmp"+tmp) 
    os.remove(local_path+filename)
    return filename+".hmmresult"


##################################################################################
# Processes result file generated from hmmscan and returns dictionary of results #
##################################################################################
def processHMMresult(filename,local_path):
    f = open(local_path+filename,'r')
    cregions = {'ns5b_5prime':'','ns5b_3prime':'','dir':''}
    
    for x in f.readlines():
        x.rstrip()
        a = x.split()
        if '#' not in x:                     #[ ali_from_corr, ali_to_corr, (ali_from,  ali_to,   env_from,  env_to) ]  
            if a[0] == 'ns5b_3prime': cregions['ns5b_3prime'] = [int(a[17])-1,int(a[18]),(int(a[15]),int(a[16]),int(a[17]),int(a[18]),int(a[19]),int(a[20]))]    # alreday corrected for python zero-based position count
            elif a[0] == 'ns5b_5prime': cregions['ns5b_5prime'] = [int(a[17])-1,int(a[18]),(int(a[15]),int(a[16]),int(a[17]),int(a[18]),int(a[19]),int(a[20]))]
            cregions['dir'] = a[3]

    f.close()
    os.remove(local_path+filename)
    return cregions

# scanSequences
###################################
def scanSequences(title,sequence,quality):
    tmp = str(random.random())[2:]
    seq = Seq(sequence)
    tempfile = open("hmm.seq"+tmp,"w")
    tempfile.write(">forward\n"+sequence+"\n>reverse\n"+seq.reverse_complement().tostring())   # writing full length seqto file for hmmscan
    tempfile.close()
    
    local_path = os.getcwd()+"/"
    hmmscan_bin = "/usr/local/bin/hmmscan"
    hmmresult_filename = doHMMScan("hmm.seq"+tmp,local_path,hmmscan_bin)
    target_cregions = processHMMresult(hmmresult_filename,local_path)
    s = q =''
    if target_cregions['ns5b_5prime'] =='' and target_cregions['ns5b_3prime'] =='':
        s = sequence
        q = quality
    elif target_cregions['ns5b_5prime'] != '' and target_cregions=='forward':
        x = target_cregions['ns5b_5prime'][0]
        s = sequence[x:]
        q = quality[x:]
    elif target_cregions['ns5b_3prime'] != '' and target_cregions=='forward':
        x = target_cregions['ns5b_3prime'][0]
        s = sequence[:x]
        q = quality[:x]
    elif target_cregions['ns5b_5prime'] and target_cregions=='reverse':
        x = len(sequence)-target_cregions['ns5b_5prime'][0]
        s = sequence[:x]
        q = quality[:x]
    elif target_cregions['ns5b_3prime'] and target_cregions=='reverse':
        x = len(sequence)-target_cregions['ns5b_3prime'][0]+1     # This is a little tricky to compensate for 0-based index of string
        s = sequence[x:]
        q = quality[x:]
    
  #  print "title: ",title," q: ",q," s: ",s
    return (title,s,q)


# assignWork
######################
def assignWork(fn):
    print "we are inside and chunk is ",fn
    local_path = os.getcwd()+'/'
    fq = open(local_path+fn,'rU')
    outfile = open(local_path+'z'+fn+'.out','w')
    before  = open(local_path+'z'+fn+'.before','w')
    after   = open(local_path+'z'+fn+'.after','w')

    for (title, sequence, quality) in FastqGeneralIterator(fq):
        (t,s,q) = scanSequences(title,sequence,quality)
        outfile.write('@'+t+'\n'+s+'\n+\n'+q+'\n')
        if s != sequence:
            before.write('>'+title+'\n'+sequence+'\n')
            after.write('>'+title+'\n'+s+'\n')
    outfile.close()
    before.close()
    after.close()
    return 


################
# processFiles #
################
def processFiles(seqfile):
    print "Beginning multiprocess QC...."
    cpus = multiprocessing.cpu_count()

    print "CPUs = ",cpus
    lc = getoutput('wc -l '+seqfile).split(' ')[1] # getoutput is deprectaed after v2.6. Should replace with subprocess.check_output(['wc','-l',seqfile]).split(' ')[0] on v2.7+ systems
   # lc = subprocess.check_output(['wc','-l',seqfile]).split(' ')[0]
    print "lc is ",lc
    reads = int(lc)/4
    size = int(round(float(reads)/cpus))*4
    print "Lines: ",lc," Reads: ",reads," Chunk size: ",size
    subprocess.call(['split','-l',str(size),seqfile])
    path =""
    file_chunks = []
    for infile in glob.glob (os.path.join(path,"x*")): file_chunks.append(infile)

    jobs = []
    for chunk in file_chunks:
        p = multiprocessing.Process(target=assignWork, args=(chunk,))
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()

    getoutput('cat z*.out > '+seqfile+'.stripped.fastq')
    getoutput('cat z*.before > '+seqfile+'.before.fasta')
    getoutput('cat z*.after > '+seqfile+'.after.fasta')
    getoutput('rm -f z*')
    getoutput('rm -f x*')

    print "QC complete\n\n"


###########
# Main    #
###########
def main(argv):
    seqfile=''
    
    # Read command line arguments
    try:
       opts, args = getopt.getopt(argv,"s:h:",["seq,help"])
    except getopt.GetoptError:
       usage()
       sys.exit(2)
    for opt, arg in opts:
       if opt in ("-h","--help"):
          usage()
          sys.exit()
       elif opt in ("-s","--seq"): seqfile = arg
   
    if seqfile !='':
       processFiles(seqfile)
    else:
       print opts
       print sys.argv
       usage()
       sys.exit()



if __name__ == '__main__':
    main(sys.argv[1:])