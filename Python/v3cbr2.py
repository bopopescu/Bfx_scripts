"""
v3cbr.py

Created by Mark Evans on 2010-11-02.
Copyright (c) 2010 __Monogram Biosciences__. All rights reserved.
"""
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from operator import itemgetter
import sys
from multiprocessing import Pool,Process,Manager, Lock, Queue
import multiprocessing
import pprint


# Declare global scope variables

class glbl:  # global variables
   mat = {}
   cases = {}
   queries = {}
   nthreads = 8
   jobs =[]
   limit = int()
   qfilename=''
   cfilename=''
   
def LoadScoringMatrix():
   matrix = 'scoring/11_25.mat'  # Identify scoring matrix to use
   m = open(matrix,'r')          # Read file
   for line in m:                # Loop file by line
      a = line.split()           # split line contents into array by whitespace
      if a[0] != '':             # ignore first line of matrix
         aa = a[0]               # copy amino acid residue to aa 
         del(a[0])               # remove 1st array position
         glbl.mat[aa] = a             # Add dictionary entry, key=amino acid, value = array values
   return

def LoadSequences(filename,seqs):
   for record in SeqIO.parse(filename,'fasta'):
      record_id = record.id.split('|')
      tropism = 'R5'
      if record_id[2] != 'R5' or (record_id[3] != 'R5' and record_id[3] !=''): tropism='DX'
      seqs[record_id[0]] = {'seq':record.seq.tostring(),'trop':tropism,'seqlen':len(record.seq.tostring())}
   return

def chunks(list,width):
   """ Yield successive width-sized chunks from list """ 
   for i in xrange(0,len(list),width):
      yield list[i:i+width]
      
      
def processCBR(nid,results):
   cbrscores = {}
   count=0
   name = multiprocessing.current_process().name
   print "Starting process ",name
   qid = ''
   while 1:
      try: qid = nid.get(False)
      except: break
#      print "qid=",qid
      for cid in glbl.cases:
#         print "cid=",cid                                                           # get one sequence record from case seq dict
         if glbl.queries[qid]['seqlen'] == glbl.cases[cid]['seqlen'] == 35 and qid != cid:  # Limit comparisons to equal length of 35 aa
            score = 0
            for i in range(0,35):                                                # Loop through each position in sequence
               try:
                  if glbl.queries[qid]['seq'][i] == glbl.cases[cid]['seq'][i]:                # if amino acids match
                     score += int(glbl.mat[glbl.cases[cid]['seq'][i]][i])                        # add score from scoring matrix to total score
               except:
                  print "+++++ Error +++++\nType: "+str(sys.exc_type)+"\nValue: "+str(sys.exc_value)+"\nTraceback: "+str(sys.exc_traceback)+"\n"
                  print "Index = "+str(i)+" cid="+cid+" seq="+glbl.cases[cid]['seq']+" qid:"+qid+" seq="+glbl.queries[qid]['seq']+"\nMat value= "+str(glbl.mat[glbl.cases[cid]['seq'][i]][i])+"\n++++++++++\n"
                  multiprocessing.current_process.terminate()
                  sys.exit()
               cbrscores[cid] = score                                              # store final score in a dictionary keyed to case sequence id
   # creates sorted copy of dictionary that is a list of tuples 
   # [('seq1',score),('seq2',score), etc] and stores in final result dictionary by query id
#   print "QID: "+qid+" CIDscores: "+str(cbrscores)
      results[qid] = sorted(cbrscores.iteritems(),key=itemgetter(1),reverse=True)
#      print "QID:"+qid+" "+str(sorted(cbrscores.iteritems(),key=itemgetter(1),reverse=True))+"\n"
      count +=1
   print "thread ",name, " exiting; processed ",count, "sequence ids"
   return


def main():
   if __name__=="__main__":
      
      glbl.qfilename = raw_input("Enter the name of the query sequence file (multi-FASTA format)? ")
      glbl.cfilename = raw_input("Enter the name of the case sequence file (multi-FASTA format?) ")
      glbl.limit = int(raw_input("Use the top N sequence matches to consider for X4ness. N = "))
      
      manager = Manager()  # creates shared memory manager object
#      mat = manager.dict()  # Dictionary
#      cases = manager.dict()
#      queries = manager.dict()
      results = manager.dict()
      nextid = Queue()
      
      LoadScoringMatrix()
      LoadSequences(glbl.cfilename,glbl.cases)
      LoadSequences(glbl.qfilename,glbl.queries)
      
      Qids = glbl.queries.keys()
#      print str(Qids)
      for qid in Qids: nextid.put(qid)
#      jobs =[]
#   pool = Pool(processes=6)                       # set up pool of potential processes, max = CPU number
      for x in range(0,glbl.nthreads):
#      pool.apply(processCBR, (qid,cases))   # Execute function asynchronously from pool
         p = Process(target=processCBR, args=(nextid,results))
         glbl.jobs.append(p)
         p.start()
      for j in glbl.jobs:
         j.join()
      
      print "================================\nTHE RESULTS\n"
      print len(results)
      print "================================\nTHE END\n"

      
      # Generate result file
      ########################################
      f = open(glbl.qfilename+".cbr_results",'w')
      f.write("QueryID\tActual Tropism\tPredicted Tropism\tScore1\tScore2\tScore3\tWorstScore\n")
      totalR5=0
      totalDX=0
      numR5Correct = 0
      numDXCorrect = 0
      for Qid in glbl.queries:                                   # Loop through results for each query
         actualT=glbl.queries[Qid]['trop']
         if actualT == 'R5': totalR5 += 1
         else: totalDX += 1 
         count = 0
         cb = 'R5'
         x=[]    
         print results[Qid][0][0],", ",results[Qid][0][1]                                           # hold scores for printing
         print "length of results[Qid] is ",len(results[Qid])
         try:
            for i in range(0,glbl.limit+1):
               print "i = ",i                     # for a query id, loop through top 10 scoring sequences
               print "checking result for Qid",Qid," ",results[Qid][i]
               if results[Qid][i][0] != Qid: 
                  print "\n cid= ",results[Qid][i][0]," score is ",results[Qid][i][1]," cases[v]= ", glbl.cases[results[Qid][i][0]]
                  if count == glbl.limit: break                        # user-defined limit for the number of sequences to check for X4-ness
                  if glbl.cases[results[Qid][i][0]]['trop'] != 'R5': cb='DX'          # if any of the top 10 are DX, set call to DX
                  x.append(results[Qid][i][1])
                  count += 1
         except:
            print "+++++ Error +++++\nType: "+str(sys.exc_type)+"\nValue: "+str(sys.exc_value)+"\nTraceback: "+str(sys.exc_traceback)+"\n"
         if actualT=='R5' and cb =='R5': numR5Correct += 1
         if actualT=='DX' and cb =='DX': numDXCorrect += 1
         print "x=",x,"\n"
         f.write("\t".join([Qid,actualT,cb]))
         f.write("\t"+str(x[0]))
         f.write("\t"+str(x[1]))
         f.write("\t"+str(x[2]))
         f.write("\t"+str(x[len(x)-1])+"\n")          # Write out results as "SeqID Actual_tropism Predicted_tropism" 
      f.close()
      print "\n\n\n========================\nFinal Results\nR5 Total\tR5 Correct\t% R5 Correct\tDX Total\tDX Correct\t% DX Correct\n"
#   print str(totalR5)+"\t"+str(numR5Correct)+"\t"+str((float(numR5Correct)/float(totalR5))*100)+"\t\t"+str(totalDX)+"\t"+str(numDXCorrect)+"\t"+str((float(numDXCorrect)/float(totalDX))*100)+"\n========================\n\n\n"
      print str(totalR5)+"\t"+str(numR5Correct)+"\t"+str((float(numR5Correct)/float(totalR5))*100)+"\t\t"+str(totalDX)+"\t"+str(numDXCorrect)+"\t"+str((float(numDXCorrect)/float(totalDX))*100)+"\n========================\n\n\n"


if __name__ == '__main__':
   main()