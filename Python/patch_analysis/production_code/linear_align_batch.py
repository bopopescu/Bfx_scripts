#!/usr/bin/env python
# encoding: utf-8
"""
linear_align_batch.py
Loop through all result files and perform linear alignment of the epitopes identified by position
Calculate a normalized percent frequency of occurance by position and use this to identify
the final epitope network for each antibody 

Created by Mark Evans on 2012-05-11.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""

import sys,os,getopt,string
from decimal import *

#############
# Usage     #
#############
def usage():
    code =  "\n\n#############\nLoop through all result files and perform linear alignment of the epitopes identified by position\n"
    code += "Calculate a normalized percent frequency of occurance by position and use this to identify the final epitope network for each antibody  \n\n"
    code += "Usage: $ python linear_align_batch.py -d data_directory_path \n\n"
    code += "-d [data_dir] path to data directory containing result files to process\n"
    code += "-h help\n\n#############\n\n"
    print code
    return


def processFiles(datadir):
    files = os.listdir(datadir)     # get all of the result file names
    MASTER_RESULTS = open(datadir+'/'+"main_batch_results.txt","w")  # Main result file, human readable
    EPITOPES = open(datadir+'/'+'chimera_epitope_input.txt',"w")       # used as input for next script, calcFinalEpitope, to generate images and the final epitope call
    MASTER_RESULTS.write("Dir1\tDir2\tDir3\tpatch_size\tmodels_used\tpatch_depth\tmolecule\tAntibody\tResult_file\tPrimary_Epitope\tSecondary_Epitope\tpri_ep_cmd\tsec_ep_cmd\n")
    
    # Loop through each result file and process all of the data for each patch    
    for f in files:
        print "processing ",f
        if f[:6] == 'result' and f[-3:] == 'txt':   # Make sure it is a result file
            linecount = 0
            INFILE = open(datadir+"/"+f,'r')
            OUT = open(datadir+'/'+f[:-4]+".tmp","w")   # create a temp file to hold all the epitope data for the next step.
            mab =''

            # Parse result files, merge epitope data into one tmp file
            for line in INFILE:
                if linecount == 2:
                    row =line.replace('*','').lstrip().rstrip().split('/')
                    del(row[0])
                    MASTER_RESULTS.write('\t'.join(row)+'\t')
                    EPITOPES.write(row[-4:-3][0]+'\t'+row[-3:-2][0]+'\t'+row[-1:][0]+'\t')
                if linecount >= 6:
                    a = line.replace("[1] ","").replace('"','').split(',')
                    if len(a) != 1:
                        if int(a[2]) == 4 and int(a[3]) >= 1:   # ====>> Key parameter: Adjusting will impact results. a[2] = # of models, a[3] = # of sig sites in patch. If the format of the result file changes, this will break
                            OUT.write(a[0]+'\t'+a[1]+'\n')
                linecount += 1
            INFILE.close()
            OUT.close()
            

            # Parse tmp file, store data in data structures
            TMP = open(datadir+'/'+f[:-4]+".tmp","r")               
            RESULT = open(datadir+'/'+f[:-4]+".linalign.txt","w") 
            dataset = {}
            pos = {}    # Keep count of how many times position is seen
            all_aa ={}  # unique set of all amino acids seen, indexed by position
            
            # load top patch result set
            for line in TMP:
                a = line.rstrip().split('\t')
                b = a[1].lstrip().split('-')
                vals = {}   # Holds alphanumeric dictionary of zero-padded positions and amino acids 
                for x in b:
                    p = int(x[1:])  # Get numeric position from alphanumeric string, i.e. 262 from K262 
                    if not all_aa.has_key(p): all_aa[p] = x
                    q =''
                    if len(str(p)) == 2: q = '0'+str(p)
                    elif len(str(p)) == 1: q = '00'+str(p)
                    else: q = str(p)
                    vals[q] = x
                    if pos.has_key(p): pos[p] = pos[p] + 1
                    else: pos[p] = 1
            
                mykey = vals.keys() # get this before adding non-numerical keys
                mykey.sort()
                vals['patch'] = a[0] # patch name
                dataset[string.join(mykey,"|")] = vals  # vals keys become unique key for dataset
            
            print f[:-4]," max length is ",len(pos)
            print f[:-4]," num of patches is ",len(dataset)
            N = len(dataset)
            k = pos.keys()
            k.sort()
            order = dataset.keys()
            order.sort()
    
            # print out all patch values by position
            for patch in order:
                for p in k: # Loop through sorted positions in the patch
                    q =''
                    if len(str(p)) == 2: q = '0'+str(p)
                    elif len(str(p)) == 1: q = '00'+str(p)
                    else: q = str(p)
                    # Write linearly aligned value
                    if dataset[patch].has_key(q): RESULT.write(dataset[patch][q]+"\t")
                    else: RESULT.write("\t")
                RESULT.write(str(dataset[patch]['patch']+"\n"))
            RESULT.write("\n")
            
            # Print out raw frequency line
            for p in k:
                RESULT.write(str(pos[p])+"\t")
            RESULT.write("Raw count\n")
    
            # Print out pct frequency line
            getcontext().prec=2
            maxfreq = 0
            for p in k:
                pctfreq = Decimal(pos[p])/Decimal(N)*Decimal(100)
                if pctfreq > maxfreq : maxfreq = pctfreq
                RESULT.write(str(pctfreq)+"\t")

            RESULT.write("Pct Freq\n")
    
            # Print out Normalized PCT frequency
            pri_epitope_res=[]
            sec_epitope_res=[]
            pri_epitope_pos=[]
            sec_epitope_pos=[]
            for p in k:
                normfreq = (Decimal(pos[p])/Decimal(N)*Decimal(100))/maxfreq*100
                if normfreq >= 70: 
                    pri_epitope_pos.append(str(p))
                    pri_epitope_res.append(all_aa[p])
                elif normfreq >= 50:    # ===>> Key parameter: default 50, user can adjust to impact result
                    sec_epitope_pos.append(str(p))
                    sec_epitope_res.append(all_aa[p])
                RESULT.write(str(normfreq)+"\t")

            RESULT.write("Normalized Pct Freq\n")

            # write epitope strings and Chimera model selection commands
            per = pep = ser = sep =''

            if len(pri_epitope_res) != 0: 
                per = ",".join(pri_epitope_res)
                pep = "select #0 :"+",".join(pri_epitope_pos)
            else: 
                per = pep = "No Result"
            
            if len(sec_epitope_res) != 0:
                ser = ",".join(sec_epitope_res)
                sep = "select #0 :"+",".join(sec_epitope_pos)
            else:
                ser = sep = "No Result"
            RESULT.write("\nPrimary Epitope:\t"+per+"\n")
            RESULT.write("Secondary Epitope:\t"+ser+"\n")
            RESULT.write("\nTo select epitope in Chimera, paste these into the command line:\n")
            RESULT.write("For primary epitope ===>\t"+pep+"\n")
            RESULT.write("For secondary epitope ===>\t"+sep+"\n")
            
            EPITOPES.write(f+'\t'+per+'\t'+ser+'\t'+','.join(pri_epitope_pos)+'\t'+','.join(sec_epitope_pos)+'\n')
            MASTER_RESULTS.write(f+'\t'+per+'\t'+ser+'\t'+pep+'\t'+sep+'\n')

            TMP.close()
            RESULT.close()
            os.remove(datadir+'/'+f[:-4]+".tmp")
    MASTER_RESULTS.close()
    EPITOPES.close()
    print "\n\n******** Analysis complete\n"

def main(argv):
    datadir=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"d:h:",["datadir:"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-d":
            datadir = arg    
            processFiles(datadir)
        elif opt == "-h": usage()
        else:
            print "There was an error, not right args\n"
            usage()
            sys.exit()
       
   

if __name__ == '__main__':
    main(sys.argv[1:])