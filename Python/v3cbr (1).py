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
from multiprocessing import Pool

# Declare global scope variables
mat = {}  # Dictionary

def LoadScoringMatrix():
   matrix = 'scoring/11_25.mat'  # Identify scoring matrix to use
   m = open(matrix,'r')          # Read file
   for line in m:                # Loop file by line
      a = line.split()           # split line contents into array by whitespace
      if a[0] != '':             # ignore first line of matrix
         aa = a[0]               # copy amino acid residue to aa 
         del(a[0])               # remove 1st array position
         mat[aa] = a             # Add dictionary entry, key=amino acid, value = array values
   return

def main():
   qfilename = raw_input("Enter the name of the query sequence file (multi-FASTA format)? ")
   cfilename = raw_input("Enter the name of the case sequence file (multi-FASTA format?) ")
   limit = int(raw_input("Use the top N sequence matches to consider for X4ness. N = "))
   results = {}
   LoadScoringMatrix()
   
   # Process all query sequences and generate CBR scores
   ######################################################
   for Q in SeqIO.parse(qfilename,'fasta'):                        # get one sequence record from query file
      cbrscores = {}
      for C in SeqIO.parse(cfilename,'fasta'):                     # get one sequence record from case file
         if len(Q.seq.tostring()) == len(C.seq.tostring()) == 35 and Q.id != C.id:  # Limit comparisons to equal length of 35 aa
            score = 0
            for i in range(0,35):                                  # Loop through each position in sequence
               try:
                  if Q.seq[i] == C.seq[i]:                         # if amino acids match
                     score += int(mat[C.seq[i]][i])                # add score from scoring matrix to total score
               except:
                  print "Index = "+str(i)+" Cid="+C.id+" C.seq="+C.seq.tostring()+" Qid:"+Q.id+" Q.seq="+Q.seq.tostring()
                  sys.exit()
            cbrscores[C.id] = score                                # store final score in a dictionary keyed to case sequence id
      # creates sorted copy of dictionary that is a list of tuples 
      # [('seq1',score),('seq2',score), etc] and stores in final result dictionary by query id
      results[Q.id] = sorted(cbrscores.iteritems(),key=itemgetter(1),reverse=True)  
   
   # Generate result file
   ########################################
   f = open(qfilename+".cbr_results",'w')
   f.write("QueryID\tActual Tropism\tPredicted Tropism\tScore1\tScore2\Score3\WorstScore\n")
   totalR5=0
   totalDX=0
   numR5Correct = 0
   numDXCorrect = 0
   for Qid in results:                                   # Loop through results for each query
      b = Qid.split('|')                                 # Split ID string to extract actual tropisms
      actualT='R5'
      if b[2] != 'R5' or (b[3] !='R5' and b[3] !=''): 
         actualT = 'DX'
         totalDX += 1
      else: totalR5 += 1
      count = 1
      cb = 'R5'
      x=[]                                               # hold scores for printing
      for Cid,score in results[Qid]:                     # for a query id, loop through top 10 scoring sequences
         if count == limit: break                        # user-defined limit for the number of sequences to check for X4-ness
         c = Cid.split('|')                              # Split ID string for top scoring sequences to check tropism
         if c[2] != 'R5' or (c[3] != 'R5' and c[3] !=''): cb='DX'  # if any of the top 10 are DX, set call to DX
         x.append(str(score))
         count += 1
         
      if actualT=='R5' and cb =='R5': numR5Correct += 1
      if actualT=='DX' and cb =='DX': numDXCorrect += 1
      print str(x)+'\n'
      f.write("\t".join([b[0],actualT,cb])+"\t"+x[0]+"|"+x[1]+"|"+x[2]+x[len(x)-1]+"\n")          # Write out results as "SeqID Actual_tropism Predicted_tropism" 
   f.close()
   print "\n\n\n========================\nFinal Results\nR5 Total\tR5 Correct\t% R5 Correct\tDX Total\tDX Correct\t% DX Correct\n"
#   print str(totalR5)+"\t"+str(numR5Correct)+"\t"+str((float(numR5Correct)/float(totalR5))*100)+"\t\t"+str(totalDX)+"\t"+str(numDXCorrect)+"\t"+str((float(numDXCorrect)/float(totalDX))*100)+"\n========================\n\n\n"
   print str(totalR5)+"\t"+str(numR5Correct)+"\t"+str((float(numR5Correct)/float(totalR5))*100)+"\t\t"+str(totalDX)+"\t"+str(numDXCorrect)+"\t"+str((float(numDXCorrect)/float(totalDX))*100)+"\n========================\n\n\n"
   
if __name__ == '__main__':
   main()