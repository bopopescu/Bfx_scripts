#!/usr/bin/env python
# encoding: utf-8
"""
predv3tropism.py
Version 1.0
Created by Mark Evans on 2011-01-25.
Revised 2011.02.21
Copyright (c) 2011 __Monogram Biosciences__. All rights reserved.
"""

####################################################################################
# Dependancies for Version 1.0                                                     #
#                                                                                  #
# Requires Perl, EMBOSS package (pepstats), HMMER 3.0 package (hmmscan), clustalw2 #
#                                                                                  #
# verifyCorrectFrame calls hmmscan which will require the following files:         #
#     - /usr/local/bin/hmmscan                                                     #
#     - mgrmv3.hmm                                                                 #
#     - mgrmv3.hmm.h3f                                                             #
#     - mgrmv3.hmm.h3i                                                             #
#     - mgrmv3.hmm.h3m                                                             #
#     - mgrmv3.hmm.h3p                                                             #
# runPepstats calls pepstats from EMBOSS package                                   #
#     - /usr/local/bin/pepstats                                                    #
# doHMMScan calls hmmscan which will require the following files:                  #
#     - /usr/local/bin/hmmscan                                                     #
#     - mgrmV3Tropo.hmm                                                            #
#     - mgrmV3Tropo.hmm.h3m                                                        #
#     - mgrmV3Tropo.hmm.h3i                                                        #
#     - mgrmV3Tropo.hmm.h3f                                                        #
#     - mgrmV3Tropo.hmm.h3p                                                        #
# loadPSSM requires the file:                                                      #
#     - x4r5_pssm.matrix                                                           #
# makeLength35 calls clustalw2                                                     #
#     - /usr/local/bin/clustalw2                                                   #
# runSVM requires the model file:                                                  #
#     - pySVM.Model.Feb072011.svm                                                  #
# runCBR requires the following files:                                             #
#     - v3CBREngine.pl                                                             #
#     - Scoring.pm                                                                 #
#     - Sequence.pm                                                                #
#     - Score.txt                                                                  #
#     - CaseLibrary.txt                                                            #
#                                                                                  #
####################################################################################

import sys,getopt,random,itertools,math
import os
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Data import IUPACData
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from libsvm.svm import *
from libsvm.svmutil import *


##################################
# Define usage syntax for script #
##################################   
def usage():
   code =  "\n\n#############\nProgram predicts V3 tropism for a V3 nucleotide sequence and accession number\n"
   code += "Usage: $ python predv3tropism.py [-a] accession_number [-s] sequence\n"
   code += "-a [temp_id] any single sequence identifier or temp id\n"
   code += "-s [sequence] nucleotide sequence string corresponding to V3 region\n"
   code += "-h help\n#############\n\n"
   print code
   return   



############################################################################
# Determine if sequence received is DNA or other (protein, gibberish, etc) #
############################################################################
def guess_if_dna(seq, thresh = 0.90, dna_letters = ['G', 'A', 'T', 'C','R','Y','S','W','K','M','D','H','B','V','N']):
   """"Guess if the given sequence is DNA.
      It's considered DNA if more than 90% of the sequence is GATCs. The threshold
      is configurable via the thresh parameter. dna_letters can be used to configure
      which letters are considered DNA; for instance, adding N might be useful if
      you are expecting data with ambiguous bases.
   """
   if isinstance(seq,Seq):
      seq = seq.data
   elif isinstance(seq, type("")) or isinstance(seq, type(u"")):
      seq = str(seq)
   else:
      raise ValueError("Do not know provided type: %s" % seq)
   seq = seq.upper()
   dna_alpha_count = 0
   for letter in dna_letters:
      dna_alpha_count += seq.count(letter)
   if (len(seq) == 0 or float(dna_alpha_count) / float(len(seq)) >= thresh):
      return True
   else:
      return False



#########################################################################################################
# Generate ambiguous codon translation table because of a translation "bug" in the Biopython module     #
# where ambiguous DNA that codes for 'I' or 'L' gets translated to 'J'. Returns dictionary of ambiguous #
# translations that do not include stops and another dictionary of the ones that do, minus the stop     #
#########################################################################################################
def generateAmbigCodonTables():
   bases = ['T','C','A','G']
   ambig_bases = ['R','Y','S','W','K','M','D','H','B','V','N']
   codons = [a+b+c for a in bases for b in bases for c in bases]
   amino_acids = "F F L L S S S S Y Y stop stop C C stop W L L L L P P P P H H Q Q R R R R I I I M T T T T N N K K S S R R V V V V A A A A D D E E G G G G".split(' ')
   codon_table = dict(zip(codons, amino_acids))
   nonstd = IUPACData.ambiguous_dna_values                  # create list of all ambiguous DNA values, includes normal bases
   std_nt = CodonTable.unambiguous_dna_by_name["Standard"]  # create normal codon table object
   
   all_ambig_trip = []
   for s in list(itertools.product(*[ambig_bases,ambig_bases,ambig_bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[bases,ambig_bases,ambig_bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[ambig_bases,bases,ambig_bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[ambig_bases,ambig_bases,bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[ambig_bases,bases,bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[bases,ambig_bases,bases])): all_ambig_trip.append(s)
   for s in list(itertools.product(*[bases,bases,ambig_bases])): all_ambig_trip.append(s)
   
   k = codon_table.keys()
   k.sort()
   
   stop_set = {}  # store final results of ambig codon and all possible translations minus the potential stop
   ambig_set = {} # store all possible translations of everthing else that does not code for stops
   for trip in all_ambig_trip:
      stop_seen = False
      ambig = "".join(trip)
      stop_set[ambig] = []
      ambig_set[ambig] = []
      for codon in list(itertools.product(*[ list(nonstd[trip[0]]), list(nonstd[trip[1]]), list(nonstd[trip[2]]) ])):  # translate ambig into all possible real combinations
         codon = "".join(codon)
         aa = codon_table[codon]
         if aa == 'stop': 
            stop_seen = True
         else: 
            stop_set[ambig].append(aa)
            ambig_set[ambig].append(aa)
      if stop_seen is False: del(stop_set[ambig])
      if stop_seen is True: del(ambig_set[ambig])
   
   return (ambig_set,stop_set)



###############################################################
# Translate all six reading frames and return in a dictionary #
###############################################################
def getReadingFrames(seq):
   s = {}
   for i in range(0,3):
      s[i] = seq[i:].translate().tostring().replace('J','X')
      s[i+3] = seq.reverse_complement()[i:].translate().tostring().replace('J','X')
   return s   



##################################################################################
# Run HMM against all frame translations and parse output to detect proper frame #
##################################################################################
def verifyCorrectFrame(seq, local_path, seqid):
   # Create multi-sequence FASTA file containing each reading frame translation for testing
   hmmscan_bin = "/usr/local/bin/hmmscan"
   f = open(local_path+seqid+".seq.tmp","w")
   for i in range(0,3):
      f.write(">frame"+str(i)+"\n"+seq[i]+"\n")
      f.write(">frame"+str(i+3)+"\n"+seq[i+3]+"\n")
   f.close()
   
   # Run hmmscan on input file and parse output file to identify correct frame
   err = open(local_path+"hmmscan_error_log","w")
   process = subprocess.Popen([hmmscan_bin,"-o",local_path+seqid+".log.tmp","--noali","--domtblout",local_path+seqid+".out.tmp",local_path+"mgrmv3.hmm",local_path+seqid+".seq.tmp"],stderr=err)
   process.wait()
   err.close()
   h = open(local_path+seqid+".out.tmp",'r')
   d = {}
   c = 1
   for line in h:
      if c > 3:
         val = line.split()
         frame = val[3]       # which frame sequence
         hfrom=int(val[15])   # hmm start pos
         hto = int(val[16])   # hmm end pos
         sfrom = int(val[17]) # sequence match start pos
         sto = int(val[18])   # sequence match end pos
         acc = float(val[21]) # accuracy
         if d.has_key(frame):
            if d[frame]['acc'] < acc:
               d[frame] = {'hfrom':hfrom,'hto':hto,'sfrom':sfrom,'sto':sto,'acc':acc}
         else: d[frame] = {'hfrom':hfrom,'hto':hto,'sfrom':sfrom,'sto':sto,'acc':acc}
      else: c+=1
   h.close()
   os.remove(local_path+seqid+".log.tmp")
   os.remove(local_path+seqid+".out.tmp")
   os.remove(local_path+seqid+".seq.tmp")
   os.remove(local_path+"hmmscan_error_log")
   return d



####################################################################################
# Once proper frame identified, translate and parse amino acids to exact V3 region #
####################################################################################
def parseDNA(hmm,record):
   err = ""
   f = hmm.keys()
   if len(f) != 0:
      frame = int(f[0][5])
      if frame < 3:  # forward reading frames
         try:
            dna  = record.seq[frame+(int(hmm[f[0]]['sfrom'])*3)-3:frame+(int(hmm[f[0]]['sto'])*3)].tostring()
            prot = record.seq[frame+(int(hmm[f[0]]['sfrom'])*3)-3:frame+(int(hmm[f[0]]['sto'])*3)].translate().tostring().replace('J','X')
            
         except:
            dna = "No Match"
            prot = "No Match"
            err += "No Match for "+record.id+"\t"+record.seq.tostring()+"\n"
      elif frame >2:  # reverse reading frames
         try:
            dna = record.seq.reverse_complement()[frame-3+int(hmm[f[0]]['sfrom'])-3:frame-3+hmm[f[0]]['sto']].tostring()
            prot = record.seq.reverse_complement()[frame-3+int(hmm[f[0]]['sfrom'])-3:frame-3+hmm[f[0]]['sto']].translate().tostring().replace('J','X')
         except:
            dna = "No Match"
            prot = "No Match"
            err += "No Match for "+record.id+"\t"+record.seq.tostring()+"\n"
   else:
      err += "FAILURE: No HMM Keys for "+record.id+"\t"+record.seq.tostring()+"\n"
      print "FAILURE: No HMM Keys for "+record.id+"\t"+record.seq.tostring()+"\n"
      dna = "No Match"
      prot = "No Match"
   return (dna,prot,err)



###############################################################################################
# Counts mutation frequencies in both DNA and AA sequence and returns a dictionary of results #
###############################################################################################
def getMutationFreq(nt,aa):
   
   a_seq = {aa.id:aa.seq}
   d_seq = {nt.id:nt.seq}
   result = {"sid":nt.id,"NT_len":0,"AA_len":0,"NT_N":0,"NT_R":0,"NT_Y":0,"NT_W":0,"NT_S":0,"NT_M":0,"NT_K":0,"NT_H":0,"NT_B":0,"NT_V":0,"NT_D":0,"AA_X":0,"AA_B":0,"AA_Z":0}
   
   a = nt.id
   result['NT_len'] = str(len(d_seq[a]))
   result['AA_len'] = str(len(a_seq[a]))
   result['NT_N'] = str(d_seq[a].count('N'))
   result['NT_R'] = str(d_seq[a].count('R'))
   result['NT_Y'] = str(d_seq[a].count('Y'))
   result['NT_W'] = str(d_seq[a].count('W'))
   result['NT_S'] = str(d_seq[a].count('S'))
   result['NT_M'] = str(d_seq[a].count('M'))
   result['NT_K'] = str(d_seq[a].count('K'))
   result['NT_H'] = str(d_seq[a].count('H'))
   result['NT_B'] = str(d_seq[a].count('B'))
   result['NT_V'] = str(d_seq[a].count('V'))
   result['NT_D'] = str(d_seq[a].count('D'))
   result['AA_X'] = str(a_seq[a].count('X'))
   result['AA_B'] = str(a_seq[a].count('B'))
   result['AA_Z'] = str(a_seq[a].count('Z'))
   
   return result



##########################################
# Calls pepstat, returns output filename #
##########################################
def runPepstats(aa_file, local_path):
   pepstats_bin = "/usr/local/bin/pepstats"
   outfilename = local_path+aa_file+".pepstat"
   err = open(local_path+"pepstats.err.log","w")
   process = subprocess.Popen([pepstats_bin,"-sequence",local_path+aa_file,"-outfile",outfilename,"-sformat1","fasta","-sprotein1","-aadata","Eamino.dat","-mwdata","Emolwt.dat","-termini","-nomono","-auto"],stderr=err )
   process.wait()
   err.close()
   os.remove(local_path+"pepstats.err.log")
   return outfilename



###############################################################################################################
# Parses pepstat output into two dictionaries, one for the amino acid DayhoffStat and one for everything else #
###############################################################################################################
def parsePepstats(pepstat_file,local_path):
   p = open(pepstat_file,"r")
   result = {'mw':0,'charge':0,'iep':0,'Tiny':0,'Small':0,'Aliphatic':0,'Aromatic':0,'Non-polar':0,'Polar':0,'Charged':0,'Basic':0,'Acidic':0}
   DayhoffStat = {'A':0,'B':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'J':0,'K':0,'L':0,'M':0,'N':0,'O':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'U':0,'V':0,'W':0,'X':0,'Y':0,'Z':0}
   lines = p.readlines()
   p.close()
   os.remove(pepstat_file)
   
   result['mw']     = lines[2].split()[3]
   result['charge'] = lines[3].split()[7]
   result['iep']    = lines[4].split()[3]
   for x in range(38,47):
      result[lines[x].split()[0]] = lines[x].split()[3]
   for y in range(10,36):
      DayhoffStat[lines[y].split()[0]] = lines[y].split()[5]
   del DayhoffStat['J']   #
   del DayhoffStat['O']   # These are never used, but inserted just to make the code simple
   del DayhoffStat['U']   #
   
   return (result,DayhoffStat)



#########################################################
# Calls hmmscan from HMMER3.0 and generates result file #
#########################################################
def doHMMScan(filename,local_path):
   hmmscan_bin = "/usr/local/bin/hmmscan"
   err = open(local_path+"hmmscan_error_log","w")
   process = subprocess.Popen([hmmscan_bin,"-o",local_path+"hmmscan.log.tmp","--noali","--tblout",local_path+filename+".hmmresult",local_path+"mgrmV3Tropo.hmm",local_path+filename],stderr=err)
   process.wait()
   err.close()
   os.remove(local_path+"hmmscan_error_log")
   os.remove(local_path+"hmmscan.log.tmp")
   return filename+".hmmresult"



##################################################################################
# Processes result file generated from hmmscan and returns dictionary of results #
##################################################################################
def processHMMresult(filename,local_path):
   f = open(local_path+filename,'r')
   hmlookup = {'TxDX':'DM','TnR5':'R5'}
   R5_Score=0
   DX_Score=0
   for x in f.readlines():
      x.rstrip()
      a = x.split()
      if '#' not in x:
         if hmlookup[a[0]] == 'R5': R5_Score = a[5]
         elif hmlookup[a[0]] == 'DM':  DX_Score = a[5]
   f.close()
   call=''
   if R5_Score != DX_Score and R5_Score > DX_Score: call ='R5'
   else: call = 'DM'
   result = {'call'    :call,
             'R5_Score':R5_Score,
             'DX_Score':DX_Score,
             'DX-R5': float(DX_Score)-float(R5_Score)}
   os.remove(local_path+filename)
   return result



###########################################
# Takes Bio.Seq.Seq object as input       #
# Returns list of all possible proteins   #
###########################################
def generateProtFromAmbiguousDNA(s,codon_table,stop_alt):
   std_nt = CodonTable.unambiguous_dna_by_name["Standard"]  # create normal codon table object
   nonstd = IUPACData.ambiguous_dna_values                  # create list of all ambiguous DNA values, includes normal bases
   aa_trans = []
   
   for i in range(0,len(s),3):
      codon = s.tostring()[i:i+3]
      if stop_alt.has_key(codon):   aa_trans.append(stop_alt[codon])
      elif codon == '---':   aa_trans.append('-')
      elif codon_table.has_key(codon):   aa_trans.append(codon_table[codon])
      elif std_nt.forward_table.has_key(codon):   aa_trans.append(list(std_nt.forward_table[codon]))
      
   return aa_trans


###################################################################################################
# Generate the most-X4-like sequence by evaluating the most X4-like codon translation by position #
# instead of translating all peptide permuations before scoring them. Much faster solution and    #
# is independant of the number of potential permutations                                          #
###################################################################################################
def getX4ByCodon(aa_trans,pssm):
   prot =""
   pos = 0
   final_score=0
   for x in aa_trans:
      if len(x) == 1: 
         if x[0] != '-':
            final_score = final_score + pssm[pos][x[0]]
         prot = prot+x[0]
      else:
         c = {}
         for y in x:
            score = pssm[pos][y]
            c[score] = y
         best = c.keys()
         best.sort()
         prot = prot + c[best[len(best)-1]]
         final_score = final_score + best[len(best)-1]
      pos = pos + 1
   return (prot,final_score)



###################################
# Load X4R5 PSSM matrix from file #
###################################
def loadPSSM(local_path):
   f = open(local_path+"x4r5_pssm.matrix","r")
   ids = {'A':0, 'C':0, 'E':0, 'D':0, 'G':0, 'F':0, 'I':0, 'H':0, 'K':0, 'M':0, 'L':0, 'N':0, 'Q':0, 'P':0, 'S':0, 'R':0, 'T':0, 'W':0, 'V':0, 'Y':0}
   pssm = []
   for x in range(0,35):
      pssm.append(ids.copy())
   
   for row in f:
      vals = row.split()
      if vals[0] != "#" and vals[0] != "1":
         for idx in range(0,35):
            pssm[idx][vals[0]] = float(vals[idx+1])
            
   f.close()
   
   return pssm



##########################################################
# Takes aa sequence and compresses or pads it to be 35aa #
##########################################################
def makeLength35(dna,seqid,local_path):
   nt = dna.seq.tostring()
   clustalw2_bin = "/usr/local/bin/clustalw2"
   #consensus = 'CTRPNNNTRKSIHIGPGRAFYATGDIIGDIRQAHC'    # Consensus of all mgrm 35mers @ 30% occurance
   consensus = 'CTRPNNNTRKSIHIGPGRAFYATGEIIGDIRQAHC'    # Subtype B consensus
   f = open(local_path+seqid+".clustal.seq",'w')
   f.write(">C\n" + consensus + "\n>Q\n" + dna.seq.translate().tostring().replace('J','X') )  # Have to add the replace function to compensate for weird biopython 'bug'
   f.close()
   
   # Create command line request for clustalw
   ##########################################
   err = open(local_path+"clustalw2_error_log","w")
   process = subprocess.Popen( [ clustalw2_bin, "-infile="+ local_path+seqid +".clustal.seq", "-outfile="+ local_path+seqid +".aln" ], stderr = err, stdout = err)
   process.wait()
   err.close()
   align = AlignIO.read(local_path+seqid+".aln",'clustal')
   summary_align = AlignInfo.SummaryInfo(align)
   seqs = { align[0].id : align[0].seq.tostring(), align[1].id : align[1].seq.tostring() }
   os.remove(local_path+"clustalw2_error_log")
   
   # Look for insertions
   ###############################
   changes ={'ins':[],'del':[]}
   if seqs['C'].find('-') != -1:                                                # Look for insertions in query "Q" by looking for '-' in consensus "C"
      gapnum = seqs['C'].count('-')                                             # count num of gap chars
      gaplength = len(seqs['C'][seqs['C'].find('-'):seqs['C'].rfind('-')+1])    # measure length from first gap char to last gap char
      if gapnum != 0:                                                           # if we have an insertion event, gapnum will not be 0
         if gapnum == gaplength:                                                # if gapnum == gaplength, then we have an uninterrupted insertion = simple case
            start = seqs['C'].find('-')                                         # Not adding 1 because we are not using the insertion info, so will preserve original aa to left. Also Z insertion was bad choice as Z is already ambig char
            end = seqs['C'].rfind('-') + 1                                      # Add one because zero based
            seqs['Q'] = seqs['Q'][:start]+seqs['Q'][end:]                       
            return_dna = nt[:start*3]+nt[end*3:]
            changes['ins'].append(str(seqs['C'].find('-')+1)+":"+str(seqs['C'].rfind('-')+1))
         else:                                                                 # multiple insertion events 
            s = list(seqs['C'])                                                # convert "C" sequence into a list
            indices = []                                                       # store list indices of '-' occurance
            flag='no'                                                          # are we in the middle of a gap?
            
            for idx,val in enumerate(s):                                       # idx = list index, val = list value
               if val =='-':
                  if flag == 'no': indices.append(idx)
                  flag = 'yes'
               elif val != '-' and flag =='yes':
                  indices.append(idx-1)
                  flag = 'no'
                  
            z = list(seqs['Q'])                                                # convert "Q" sequence into a list
            d = list(nt)                                                       # Convert dna to list
            for idx in range(0,len(indices),2):                                # loop through indices of gaps
               changes['ins'].append(str(indices[idx]+1)+":"+str(indices[idx+1]+1))
               for xidx in range(indices[idx],indices[idx+1]+1):               # loop through length of insertion
                  z[xidx]="|"                                                  # set each value in the insertion range to a 'throw-away' string
                  d[xidx*3]="|"
                  d[xidx*3+1]="|"
                  d[xidx*3+2]="|"                                              
            newseq = "".join(z)                                                # Join the list back into a string
            newseq = newseq.replace("|","")                                    # remove all throw-away chars
            newnt = "".join(d)
            newnt = newnt.replace("|","")
            return_dna = newnt
            seqs['Q'] = newseq
   else:
      gapnum = seqs['Q'].count('-')
      gaplength = len(seqs['Q'][seqs['Q'].find('-'):seqs['Q'].rfind('-')+1])
      if gapnum != 0:                                                            # if we have an deletion event, gapnum will not be 0
         if gapnum == gaplength:                                                 # if gapnum == gaplength, then we have an uninterrupted deletion = simple case
            start = seqs['Q'].find('-')                                          # Not adding 1 because we are not using the insertion info, so will preserve original aa to left. Also Z insertion was bad choice as Z is already ambig char
            end = seqs['Q'].rfind('-')                                           
            return_dna = nt[:start*3]+("---"*gapnum)+nt[end*3:]
            changes['del'].append(str(seqs['Q'].find('-')+1)+":"+str(seqs['Q'].rfind('-')+1))
         else:                                                                 
            s = list(seqs['Q'])                                                
            indices = []                                                       
            flag='no'                                                          
            for idx,val in enumerate(s):                                       
               if val =='-':
                  if flag == 'no': indices.append(idx)
                  flag = 'yes'
               elif val != '-' and flag =='yes':
                  indices.append(idx-1)
                  flag = 'no'
            if flag == 'yes': indices.append(idx)
            
            
            d = list(nt)                                                       # Convert dna to list   
            delcount = 0
            lastidx = range(0,len(indices),2)[len(range(0,len(indices),2))-1]   # get last value of what idx will be in case we need it
            for idx in range(0,len(indices),2):                                 # loop through indices of gaps
               changes['del'].append(str(indices[idx]+1)+":"+str(indices[idx+1]+1))
               for xidx in range(indices[idx],indices[idx+1]+1):                        # loop through length of insertion
                  try:
                     d[(xidx-delcount)*3]= "---"+d[(xidx-delcount)*3]   
                  except:
                     if idx == lastidx:
                        d.append("---")
                  delcount += 1
            newnt = "".join(d)
            return_dna = newnt
      
   os.remove( local_path+seqid +".clustal.seq")   # remove all temp files
   os.remove( local_path+seqid +".aln")
   os.remove( local_path+seqid +".clustal.dnd")
   return (seqs['Q'],return_dna,changes)



#########################################################
# Runs SVM prediction, returns SVM call and SVM score   #
#########################################################
def runSVM(local_path,aa,mut_freq):
   # Make 0/1 vector
   # Libsvm requires data format to be {param_number1 : param_value1, param_number2 : param_value2, etc}
   vect = {}
   seq = list(aa)  # Splits sequence into array of letters
   ids = {'-':(0,1),'A':(0,1), 'C':(0,1), 'E':(0,1), 'D':(0,1), 'G':(0,1), 'F':(0,1), 'I':(0,1), 'H':(0,1), 'K':(0,1), 'M':(0,1), 'L':(0,1), 'N':(0,1), 'Q':(0,1), 'P':(0,1), 'S':(0,1), 'R':(0,1), 'T':(0,1), 'W':(0,1), 'V':(0,1), 'Y':(0,1),'Z':(0,1)}
   alpha = ids.keys()
   alpha.sort()
   n = 0
   for s in seq:
      for a in alpha:
         if s == a: 
            vect[n] = 1
         else: 
            vect[n] = 0
         n += 1
   vect[770] = (float(mut_freq["NT_R"]))
   vect[771] = (float(mut_freq["NT_Y"]))
   vect[772] = (float(mut_freq["NT_W"]))
   vect[773] = (float(mut_freq["NT_K"]))
   vect[774] = (float(mut_freq["AA_X"]))
   
   # Begin SVM
   x = [vect]     # formats x into [{0:1,1:3, etc}], required by libsvm
   y = [0]*len(x) # creates an empty list of zeros of length = len(x)
   m = svm_load_model(local_path+'pySVM.Model.Feb072011.svm')
   p_label, p_acc, p_val = svm_predict(y, x, m)
   
   svm_call = ""
   svm_score = 0
   if p_label[0] == 0.0: svm_call = "X4"
   elif p_label[0] == 1.0: svm_call = "R5"
   svm_score = p_val[0][0]
   
   return (svm_call,svm_score)



#####################################################
# Run CBR in Perl, parses output and cleans up      #
# Returns CBR Call, Score and Case ID of best match #
#####################################################
def runCBR(filename,local_path):
   cbr_out = open(local_path+filename+".cbrout","w")
   process = subprocess.Popen(['perl',local_path+'v3CBREngine.pl','-q',local_path+filename] , stdout=cbr_out)
   process.wait()
   cbr_out.close()
   f = open(local_path+filename+".cbrout","r")
   c = 0
   cbr_call=""
   cbr_score = 0
   caseid =""
   for x in f.readlines():
      if c == 1:
         a = x.split()
         cbr_call = a[1]
         cbr_score = a[2]
         caseid = a[3]
      c += 1
   f.close()
   os.remove(local_path+filename+".cbrout")
   return (cbr_call,cbr_score,caseid)



##########################################################################
# Perform decision tree on svm_score, cbr_score and hmm_difference_score #
# and return the final traopism call                                     #
##########################################################################
def runDecisionTree(svm_score,cbr_score,hmm_score):
   final_call = "UNK"
   svm_score = float(svm_score)
   cbr_score = float(cbr_score)
   hmm_score = float(hmm_score)
   if svm_score > 0: final_call = "X4"
   else:
      if cbr_score < 110:
         if svm_score < -1:
            final_call = "R5"
         elif hmm_score < -2:
            final_call = "R5"
         else: final_call = "X4"
      else: final_call = "X4"
   return final_call



#############################################
# Main control for processing V3 prediction #
#############################################
def processSequence(ntseq,seqid):
   # Set environment variables
   if seqid != '0': os.chdir("v3")
   algo_ver = "Version 1.0"
#   os.chdir("v3")
   local_path = os.getcwd()+"/"
   (codon_table,stop_alt) = generateAmbigCodonTables()
   pssm = loadPSSM(local_path)
   
   # Error Checking to verify DNA was submitted
   flag = guess_if_dna(ntseq)
   if flag is False: 
      f = open(seqid+".tropo","w")
      f.write("<b><font color='darkred'>Error&nbsp;</font></b>:The sequence submitted must be a DNA sequence")
      f.close()
      return  
   
   # get nt seq and identify proper frame by HMM
   nt = SeqRecord(Seq(ntseq))
   nt.id = seqid
   frames = getReadingFrames(nt.seq)
   
   # translate and trim to C-C for V3
   hmm = verifyCorrectFrame(frames,local_path,seqid)
   dna,prot,err = parseDNA(hmm,nt)
   
   # Create new sequence objects and write seq to file for later
   trimmed_nt = SeqRecord(Seq(dna))
   trimmed_nt.id = seqid
   trimmed_aa = SeqRecord(Seq(prot))
   trimmed_aa.id = seqid
   f = open(seqid+".trimmed.seq","w")
   f.write(">"+trimmed_aa.id+"\n"+trimmed_aa.seq.tostring())
   f.close()
   trimmed_aa_file = seqid+".trimmed.seq"
   
   # Run mutation frequency analysis on nt and aa seq
   mut_freq = getMutationFreq(trimmed_nt,trimmed_aa)
   
   # Run pepstats, parse results, catch returned dictionaries with results
   pepstat_filename = runPepstats(trimmed_aa_file,local_path)
   pepstats,dayhoff_stats = parsePepstats(pepstat_filename,local_path)
   
   # run HMM analysis
   result_file = doHMMScan(trimmed_aa_file,local_path)
   hmmresults = processHMMresult(result_file,local_path)
   
   # force to 35aa and bias to X4-like
   final_dna =""
   if len(trimmed_aa) == 35:
      aa_trans = generateProtFromAmbiguousDNA(trimmed_nt.seq,codon_table,stop_alt)
      x4aa,pssm_score = getX4ByCodon(aa_trans,pssm)
      final_dna = trimmed_nt
   else:
      aa35len,x4seqnt,changes = makeLength35(trimmed_nt,seqid,local_path)
      nt35 = SeqRecord(Seq(x4seqnt, IUPAC.IUPACAmbiguousDNA()))
      aa_trans = generateProtFromAmbiguousDNA(nt35.seq,codon_table,stop_alt)
      x4aa,pssm_score = getX4ByCodon(aa_trans,pssm)
      final_dna = nt35
   
   # write X4-biased seq to file for future steps
   x4aa35_file = seqid+".35.seq"
   f = open(x4aa35_file,"w")
   f.write(">"+seqid+"\n"+x4aa)
   f.close()
   
   # Create input file for SVM, CBR
   f = open(seqid+".txt","w")
   f.write("mgrm_acc\tTropism\tseq_len35\tnt_r\tnt_y\tnt_w\tnt_k\taa_x\tCharge\tIEP\tSmall\tBasic\ttnr5_hmm\n")
   f.write(seqid+"\tUNK\t"+x4aa+"\t"+"\t".join([str(mut_freq['NT_R']),str(mut_freq['NT_Y']),str(mut_freq['NT_W']),str(mut_freq['NT_K']),str(mut_freq['AA_X']),str(pepstats['charge']),str(pepstats['iep']),str(pepstats['Small']),str(pepstats['Basic']),str(hmmresults['R5_Score'])])+"\n")
   f.close()
   filename = seqid+".txt"
   
   # run SVM
   (svm_call,svm_score) = runSVM(local_path,x4aa,mut_freq)
   
   # run CBR
   (cbr_call,cbr_score,caseid) = runCBR(filename,local_path)
   
   # Clean up temp files
   os.remove(filename)
   os.remove(x4aa35_file)
   os.remove(trimmed_aa_file)
   
   # perform decision tree and make final call
   final_call = runDecisionTree(svm_score,cbr_score,hmmresults['DX-R5'])
   
   # assemble output and return
   f = open(seqid+".tropo","w")
   f.write("<b>final_tropism_call:"+final_call+"</b><br>\n")
   f.write("Algorithm:"+algo_ver+"<br>\n")
   f.write("raw_trimmed_nt:"+dna+"<br>\n")
   f.write("raw_aa:"+prot+"<br>\n")
   f.write("35len_dna:"+final_dna.seq.tostring()+"<br>\n")
   f.write("x4biased35aa:"+x4aa+"<br>\n")
   f.write("x4hmm_score:"+str(hmmresults['DX_Score'])+"<br>\n")
   f.write("r5hmm_score:"+str(hmmresults['R5_Score'])+"<br>\n")
   f.write("DX-R5:"+str(hmmresults['DX-R5'])+"<br>\n")
   f.write("HMM call:"+hmmresults['call']+"<br>\n")
   f.write("PSSM_Score:"+str(pssm_score)+"<br>\n")
   f.write("SVM_score:"+str(svm_score)+"<br>\n")
   f.write("SVM_call:"+svm_call+"<br>\n")
   f.write("CBR_score:"+str(cbr_score)+"<br>\n")
   f.write("CBR_call:"+cbr_call+"<br>\n")
   f.write("AA_X:"+str(mut_freq['AA_X'])+"<br>\n")
   f.write("NT_R:"+str(mut_freq['NT_R'])+"<br>\n")
   f.write("NT_Y:"+str(mut_freq['NT_Y'])+"<br>\n")
   f.write("NT_W:"+str(mut_freq['NT_W'])+"<br>\n")
   f.write("NT_K:"+str(mut_freq['NT_K'])+"<br>\n")
   f.write("net_charge:"+str(pepstats['charge'])+"<br>\n")
   f.write("iep:"+str(pepstats['iep'])+"<br>\n")
   f.write("basic_m_pct:"+str(pepstats['Basic'])+"<br>\n")
   f.write("small_m_pct:"+str(pepstats['Small'])+"<br>\n")
   f.close()



#################################################
# Main                                          #
#################################################   
def main(argv):
   seqid=''
   ntseq=''
   
   # Read command line arguments
   try:
      opts, args = getopt.getopt(argv,"a:s:h:f",["accession,seq,help,filename"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt in ("-h","--help"):
         usage()
         sys.exit()
      elif opt in ("-a","--accession"): seqid = arg
      elif opt in ("-s","--seq"): ntseq = arg
      elif opt in ("-f","--filename"): seqfilename = arg
   
   if seqid !='' and ntseq !='' and seqfilename =='':
      processSequence(ntseq,seqid)
   else:  # if cmd line flags were not used, check for right args
      if len(sys.argv) == 3 and sys.argv[1] !='' and sys.argv[2] !='':
         seqid = sys.argv[1]
         ntseq = sys.argv[2]
         processSequence(ntseq,seqid)
      else:
         usage()
         sys.exit()


if __name__ == '__main__':
	main(sys.argv[1:])

