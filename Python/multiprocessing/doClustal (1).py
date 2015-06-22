"""
v3cbr3.py
Case-based reasoning analysis for V3 amino acid sequences
- currently only handles 35aa sequences, other lengths are ignored

Multiprocess-based code for symmetric multiprocessing
on a single machine with multi-CPU or multi-core CPU configuration.
Will not properly run on cluster machine, in that it will not use the available cluster
CPUs, only those available on the head node.  Needs to be converted to message passing
based code for that to work properly.

Created by Mark Evans on 2010-11-09.
Copyright (c) 2010 __Monogram Biosciences__. All rights reserved.
"""
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from operator import itemgetter
import sys
from multiprocessing import Process, Manager, Queue
import multiprocessing

# Clustal command line
# clustalw2 -INFILE=HA13_protein_seq.txt -ALIGN -QUICKTREE -TYPE=PROTEIN -OUTFILE=new_HA13_protein_seq.aln

#####################################################################################
# glbl                                                                              #
# Creates a class containing global variables that can be accessed across processes # 
#####################################################################################
class glbl:
   mat = {}
   cases = {}
   queries = {}
   jobs =[]
   limit = int()
   qfilename=''
   cfilename=''

################################################################################################
# LoadScoringMatrix                                                                            #
# Loads system defined scoring matrix into a dictionary where the keys                         #
# are amino acids and the values are lists containing a score for each position                #
# The list positions correspond to the sequence string index, so they can be accessed directly #
################################################################################################
def LoadScoringMatrix():
   matrix = 'scoring/11_25.mat'  # Identify scoring matrix to use
   m = open(matrix,'r')          # Read file
   for line in m:                # Loop file by line
      a = line.split()           # split line contents into array by whitespace
      if a[0] != '':             # ignore first line of matrix
         aa = a[0]               # copy amino acid residue to aa 
         del(a[0])               # remove 1st array position
         glbl.mat[aa] = a        # Add dictionary entry, key=amino acid, value = array values
   return

##########################################################################################
# LoadSequences                                                                          #
# Reads a multi-sequence FASTA file from filename and parses each FASTA record.          #
# Parses the defline of the record to get the id and tropism. Calculates sequence length #
# Writes data to seqs dictionary with key=seq id and value=dictionary containing the     #
# sequence, tropism and seqlength                                                        #
# @param: filename, empty dictionary                                                     #
# @return: Dictionary updated in place, not explicitly returned                          #
##########################################################################################
def LoadSequences(filename,seqs):
   for record in SeqIO.parse(filename,'fasta'):
      record_id = record.id.split('|')
      tropism = 'R5'
      if record_id[2] != 'R5' or (record_id[3] != 'R5' and record_id[3] !=''): tropism='DX'
      seqs[record_id[0]] = {'seq':record.seq.tostring(),'trop':tropism,'seqlen':len(record.seq.tostring())}
   return

##################################################################################################
# processCBR                                                                                     #
# Calculates identity matrix score all of the query sequences against each of the case sequences #
# Uses multiprocess queue object to get the next query sequence id to test. Accesses sequence    #
# data from glbl class dictionaries                                                              #
# @param: multiprocess queue object, multiprocess accessible dictionary for results              #
# @return: Dictionary is updated in place, not explicitly returned                               #
##################################################################################################
def processCBR(nid,results):
   cbrscores = {}                                                                           # Hold scores for each case sequence until done testing a query seq
   count=0
   qid = ''
   while 1:
      try: qid = nid.get(False)                                                             # Get next query id from Queue generator, 
      except: break                                                                         # unless it is empty, then exit
      for cid in glbl.cases:                                                                # loop through sequence records in case seq dict
         if glbl.queries[qid]['seqlen'] == glbl.cases[cid]['seqlen'] == 35 and qid != cid:  # Limit comparisons to equal length of 35 aa
            score = 0
            for i in range(0,35):                                                           # Loop through each position in sequence
               try:
                  if glbl.queries[qid]['seq'][i] == glbl.cases[cid]['seq'][i]:              # if amino acids match
                     score += int(glbl.mat[glbl.cases[cid]['seq'][i]][i])                   # add score from scoring matrix to total score
               except:
                  print "+++++ Error +++++\nType: "+str(sys.exc_type)+"\nValue: "+str(sys.exc_value)+"\nTraceback: "+str(sys.exc_traceback)+"\n"
                  # Next line is for debugging, can eventually be removed
                  print "Index = "+str(i)+" cid="+cid+" seq="+glbl.cases[cid]['seq']+" qid:"+qid+" seq="+glbl.queries[qid]['seq']+"\nMat value= "+str(glbl.mat[glbl.cases[cid]['seq'][i]][i])+"\n++++++++++\n"
                  multiprocessing.current_process.terminate()                               # since we are running multiple processes, kill this process. Program will continue
               cbrscores[cid] = score                                                       # store final score in a dictionary keyed to case sequence id
      # creates sorted copy of dictionary that is a list of tuples 
      # [('seq1',score),('seq2',score), etc] and stores in final result dictionary by query id
      results[qid] = sorted(cbrscores.iteritems(),key=itemgetter(1),reverse=True)
      count +=1
   return

######################################################################################
# main                                                                               #
# Home of main control thread. Will only execute once in main process                #
# and will not be called in child processes becuase of __name__="__main__" statement #
######################################################################################                
def main():
   if __name__=="__main__":                                   # This statement prevents this code from executing in child processes
      # Get user variables
      #######################
      glbl.qfilename = raw_input("Enter the name of the query sequence file (multi-FASTA format)? ")
      glbl.cfilename = raw_input("Enter the name of the case sequence file (multi-FASTA format?) ")
      glbl.limit = int(raw_input("Use the top N sequence matches to consider for X4ness. N = "))
      
      # Load scoring matrix and sequence dictionaries
      ################################################
      LoadScoringMatrix()
      LoadSequences(glbl.cfilename,glbl.cases)
      LoadSequences(glbl.qfilename,glbl.queries)
      Qids = glbl.queries.keys()
      
      # Multiprocessing Section
      #########################################
      manager = Manager()                                      # creates shared memory manager object
      results = manager.dict()                                 # Add dictionary to manager, so it can be accessed across processes
      nextid = Queue()                                         # Create Queue object to serve as shared id generator across processes
      for qid in Qids: nextid.put(qid)                         # Load the ids to be tested into the Queue
      for x in range(0,multiprocessing.cpu_count()):           # Create one process per logical CPU
         p = Process(target=processCBR, args=(nextid,results)) # Assign process to processCBR function, passing in the Queue and shared dictionary
         glbl.jobs.append(p)                                   # Add the process to a list of running processes
         p.start()                                             # Start process running
      for j in glbl.jobs:
         j.join()                                              # For each process, join them back to main, blocking on each one until finished
      
      # Generate result file
      ########################################
      f = open(glbl.qfilename+".cbr_results",'w')
      f.write("QueryID\tActual Tropism\tPredicted Tropism\tScore1\tScore2\tScore3\tWorstScore\n")
      totalR5=0
      totalDX=0
      numR5Correct = 0
      numDXCorrect = 0
      for Qid in glbl.queries:                                  # Loop through results for each query
         actualT=glbl.queries[Qid]['trop']                      # Assign actual tropism call for the query sequence
         if actualT == 'R5': totalR5 += 1
         else: totalDX += 1 
         count = 0
         cb = 'R5'
         x=[]    
         try:
            for i in range(0,glbl.limit+1):                     # Loop through results 'limit' times
               if results[Qid][i][0] != Qid: 
                  if count == glbl.limit: break                 # user-defined limit for the number of sequences to check for X4-ness
                  if glbl.cases[results[Qid][i][0]]['trop'] != 'R5': cb='DX'     # if any of the top N are DX, set call to DX
                  x.append(results[Qid][i][1])                  # add actual scores to x for reporting in output file
                  count += 1
         except:
            print "+++++ Error +++++\nType: "+str(sys.exc_type)+"\nValue: "+str(sys.exc_value)+"\nTraceback: "+str(sys.exc_traceback)+"\n"
         if actualT=='R5' and cb =='R5': numR5Correct += 1
         if actualT=='DX' and cb =='DX': numDXCorrect += 1
         f.write("\t".join([Qid,actualT,cb])+"\t"+str(x[0])+"\t"+str(x[1])+"\t"+str(x[2])+"\t"+str(x[len(x)-1])+"\n")  # Write out results as "SeqID" "Actual_tropism" "Predicted_tropism" 
      f.close()
      print "\n\n\n========================\nFinal Results\nR5 Total\tR5 Correct\t% R5 Correct\tDX Total\tDX Correct\t% DX Correct\n"
      print str(totalR5)+"\t"+str(numR5Correct)+"\t"+str((float(numR5Correct)/float(totalR5))*100)+"\t\t"+str(totalDX)+"\t"+str(numDXCorrect)+"\t"+str((float(numDXCorrect)/float(totalDX))*100)+"\n========================\n\n\n"


if __name__ == '__main__':
   main()