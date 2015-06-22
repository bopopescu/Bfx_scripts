#!/usr/bin/env python
# encoding: utf-8
"""
getHIVmutations.py
Take HIV env sequence and generate a mutation list based on HxB2

This version tries to generate mutation list from full length HXB2, rather than just constant regions

Created by Mark Evans on 2012-02-22.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

"""

import sys,getopt
import os, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import AlignInfo
from Bio import AlignIO
import re


#############
# Usage     #
#############
def usage():
   code =  "\n\n#############\nProgram indentifies mutations relative to HXB2 and generates a mutation list \n\n"
   code += "Usage: $ python getHIVmutations.py [-s] seq_filename \n\n"
   code += "-s [seq] nucleotide sequence string corresponding to gp160 region\n"
   code += "-h help\n\n#############\n\n"
   print code
   return


#########################################################
# Calls hmmscan from HMMER3.0 and generates result file #
#########################################################
def doHMMScan(filename,local_path,hmmscan_bin):
   #hmmscan_bin = "/usr/local/bin/hmmscan"
   err = open(local_path+"hmmscan_error_log","w")
   process = subprocess.Popen([hmmscan_bin,"-o",local_path+"hmmscan.log.tmp","--noali","--domT","20","--domtblout",local_path+filename+".hmmresult",local_path+"hxb2_constant_regions.hmm",local_path+filename],stderr=err)
   process.wait()
   err.close()
   os.remove(local_path+"hmmscan_error_log")
   os.remove(local_path+"hmmscan.log.tmp") 
   os.remove(local_path+filename)
   return filename+".hmmresult"


##################################################################################
# Processes result file generated from hmmscan and returns dictionary of results #
##################################################################################
def processHMMresult(filename,local_path):
   f = open(local_path+filename,'r')
   cregions = {'C1':'','C3':'','C4':'','C5':'','C6':''}
   
   for x in f.readlines():
      x.rstrip()
      a = x.split()
      if '#' not in x:                     #[ ali_from_corr, ali_to_corr, (ali_from,  ali_to,   env_from,  env_to) ]      
         if a[0] == 'C1.1': cregions['C1'] = [int(a[17])-1,int(a[18]),(int(a[15]),int(a[16]),int(a[17]),int(a[18]),int(a[19]),int(a[20]))]    # alreday corrected for python zero-based position count
         elif a[0] == 'C3': cregions['C3'] = [int(a[17])-1,int(a[18]),(int(a[15]),int(a[16]),int(a[17]),int(a[18]),int(a[19]),int(a[20]))]
         elif a[0] == 'C4': cregions['C4'] = [int(a[17])-1,int(a[18]),(int(a[15]),int(a[16]),int(a[17]),int(a[18]),int(a[19]),int(a[20]))]
         elif a[0] == 'C5': cregions['C5'] = [int(a[17])-1,int(a[18]),(int(a[15]),int(a[16]),int(a[17]),int(a[18]),int(a[19]),int(a[20]))]
         elif a[0] == 'C6': cregions['C6'] = [int(a[17])-1,int(a[18]),(int(a[15]),int(a[16]),int(a[17]),int(a[18]),int(a[19]),int(a[20]))]
   f.close()
   os.remove(local_path+filename)
   return cregions


#########################################################################################
# processAlignment                                                                      #
# Requires local path and start position. This position is used as an offset            #
# when calculating relative position to be recorded in mut_sum, instead of abs_pos      #
# used in actual calculation.  For GP41, just sending 0 as start, since this is treated #
# as seperate molecule event though it is parsed out of gp160                           #
#########################################################################################
def processAlignment(local_path,start,idxstart,molec):
   ref_aln = ''
   q_aln = ''
   mut_pos = {}
   mut_sum =''
   idx = idxstart - 1
   
   # Call Clustalw2
   subprocess.call(["clustalw2","-align","-infile="+local_path+"foo.seq", "-outfile="+local_path+"foo.aln"])
   
   # Parse the alignment file and get the aligned seq strings.  Should now be equal length
   align = AlignIO.read(local_path+'foo.aln','clustal')
   for seq in align:
      if   seq.id == 'Ref': ref_aln = seq.seq.tostring()
      elif seq.id == 'Query': q_aln = seq.seq.tostring()
   os.remove('foo.seq')
   os.remove('foo.aln')
   os.remove('foo.dnd')
   
   # Loop through each position in the sequence and compare query to ref.  If it is different, record in mut_pos
   insert = 'no'
   for pos in range(0,len(ref_aln)):
      
      if insert != 'no' and ref_aln[pos] != '-' and q_aln[pos] != '-': insert = 'no'   # Turn off insert flag

      # Case 1: residues don't match, no gaps present
      if ref_aln[pos].upper() != q_aln[pos].upper() and ref_aln[pos] != '-' and q_aln[pos] != '-':  # Record simple mutation
         idx += 1
         mut_pos[pos] = ref_aln[pos]+str(idx)+q_aln[pos]
      
      # Case 1.5: gp41 query is longer at start than ref, skip until ref begins
      elif insert=='bad ref start' and ref_aln[pos].upper() == '-':
         continue


      # Case 2: Query contains insertion
      elif ref_aln[pos] == '-' and insert == 'no' and pos != 0:                   
         
         # check that a mutation did not occur just before insertion
         if mut_pos.has_key(pos-1):                                               
            mut_pos[pos-1] = mut_pos[pos-1]+"/B"                                  # if it did, add insertion symbol 'B' to the mutation that is already recorded
         else:
            mut_pos[pos-1] = ref_aln[pos-1]+str(idx)+'B'                  # if it didn't, Step back one position in ref and use 'B' to indicate insertion occurs after this position
         
         # turn insertion flag on. Insertion residues will be skipped until end of insertion is reached
         insert = 'insertion'

      # Case 3: Query contains a deletion
      elif q_aln[pos] == '-' and insert == 'no' and pos != 0:
      	 
         # check that mutation did not occur just before deletion
         if mut_pos.has_key(pos-1):                                               
            mut_pos[pos-1] = mut_pos[pos-1]+"/^"                                  # if it did, add deletion '^' symbol to it
         else:
            mut_pos[pos-1] = ref_aln[pos-1]+str(idx)+"^"                  # if it did not, Step back one position in ref and use '^' to indicate a deletion occurs after this position
         
         # turn insertion flag on. Deletion residues will be skipped until end of deletion is reached
         insert = 'deletion'
         idx += 1

      # Case 3.5: Query contains deletion from start in gp41
      elif insert == 'bad q start' and q_aln[pos].upper()=='-':
         mut_pos[pos-1] = ref_aln[0]+str(idx)+"^"
         idx += 1
         insert = "deletion"

      # Case 4: Continue to be within an insertion or deletion
      elif insert == 'deletion' and (q_aln[pos] == '-' or ref_aln[pos] == '-'):
         idx += 1

      # Case 5: Everything matches, not in insertion or deletion
      elif ref_aln[pos].upper() == q_aln[pos].upper() and ref_aln[pos] != '-' and q_aln[pos] != '-':
         idx += 1

      # Case 6: First position is a '-' in either ref or q aln
      elif insert == 'no' and pos == 0 and (ref_aln[pos] == '-' or q_aln[pos] == '-'):
         if ref_aln[pos].upper()=='-': insert = 'bad ref start'
         elif ref_aln[pos].upper() != '-' and q_aln[pos].upper() == '-':
            idx += 1
            insert = 'bad q start'
         else: continue

   
   # Dump mut_pos contents into string to be written as part of the output
   pos_keys = mut_pos.keys()
   pos_keys.sort()
   for pos in pos_keys:
      mut_sum += mut_pos[pos]+", "
   
   return mut_sum


#########################
# countPNGs             #
#########################
def countPNGS(seq):
   seqclean = seq.strip().replace(' ','').replace('-','').replace(':','').replace('*','')
   pngs_len = len(re.findall(r'([N][^P][ST])',seqclean))
   return pngs_len


###################################################################
# getMutations                                                    #
# Collapse all insertions and deletions to fit in HXB2 numbering  #
###################################################################
def getMutations(x):
   #######################################
   # HXB2 references sequences from LANL #
   #######################################
   gp160 = "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKRAVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL"
   gp120 = "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR"
   gp41  = "AVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL"
   # HXB2 33-130 Constant region 1
   #gp120c1 = "KLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLK"
   # HXB2 196-295 Constant region 3
   #gp120c3 = "CNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEIN"
   # HXB2 331-384 Constant region 4
   #gp120c4 = "CNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFY"
   # HXB2 418-459 Constant region 5
   #gp120c5 = "CRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGG"
   # HXB2 470-511 Constant region 6
   #gp120c6 = "PGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR"

   # HXB2 gp160 constant region coords, adjusted to be zero-based
   # coords = {'C1':[32,130],'C3':[195,295],'C4':[330,384],'C5':[417,459],'C6':[469,511]} 

   # const_seqs = {'C1':gp120c1,'C3':gp120c3,'C4':gp120c4,'C5':gp120c5,'C6':gp120c6}
   # r = ['C1','C3','C4','C5','C6']
   # idx_start_pos = {'C1':{1:33, 2:33, 3:34, 4:35, 5:36},  # C1 HMM is starts 1 before C1 ref, so offset this for idx_start
   #                  'C3':{1:196,2:197,3:198,4:199,5:200},
   #                  'C4':{1:331,2:332,3:333,4:334,5:335},
   #                  'C5':{1:418,2:419,3:420,4:421,5:422},
   #                  'C6':{1:470,2:471,3:472,4:473,5:474}}
   
   #################################
   # Begin process for gp120       #
   #################################
   tempfile = open("hmm.seq","w")
   tempfile.write(">"+x.id+"\n"+x.seq.tostring()+"\n")   # writing full length seqto file for hmmscan
   tempfile.close()
   local_path = os.getcwd()+"/"
   hmmscan_bin = "/usr/local/bin/hmmscan"
   hmmresult_filename = doHMMScan("hmm.seq",local_path,hmmscan_bin)
   target_cregions = processHMMresult(hmmresult_filename,local_path)
   gp120_mut_sum=''
   
   
   # Create temp file for Clustal input
   start = 0
   stop  = target_cregions['C6'][1]
   content = ">Ref\n"+gp120+"\n>Query\n"+x.seq.tostring().strip('*')[start:stop]+"\n"
   print start,stop,content
   tempfile = open('foo.seq','w')
   tempfile.write(content)
   tempfile.close()
   
   gp120_len = str(len(x.seq.tostring().strip('*')[start:stop]))   

   gp120_mut_sum = processAlignment(local_path,start,1,'gp120')
      
   gp120_mut_sum = gp120_mut_sum.strip()  # remove trailing whitespace
   gp120_mut_sum = gp120_mut_sum[0:-1]    # remove trailing ,
   gp120_pngs = str(countPNGS(x.seq.tostring().strip('*')[start:stop]))
   
   ####################
   # Process for gp41 #
   ####################
   gp41_mut_sum =''
   start = target_cregions['C6'][1]
   target_gp41 = x.seq.tostring().strip('*')[start:]
   content = ">Ref\n"+gp41+"\n>Query\n"+target_gp41+"\n"
   tempfile = open('foo.seq','w')
   tempfile.write(content)
   tempfile.close()
   
   gp41_len = str(len(x.seq.tostring().strip('*')[start:]))

   gp41_mut_sum = processAlignment(local_path, 0, 1,'gp41')
   gp41_mut_sum = gp41_mut_sum.strip()  # remove trailing whitespace
   gp41_mut_sum = gp41_mut_sum[0:-1]    # remove trailing ,
   gp41_pngs = str(countPNGS(x.seq.tostring().strip('*')[start:]))

   return (gp120_mut_sum, gp41_mut_sum, gp120_len, gp41_len, gp120_pngs, gp41_pngs)


##########################################################
# makeVector                                             #
# Takes mutation list output and converts to 0//1 vector #
##########################################################
def makeVector(mut_file):
   
   f1 = open(mut_file,'r')
   f2 = open(mut_file+".vector.out",'w')

   positions = {}
   data = {}
   misc = {}  # hold length and PNG number

   for line in f1:
      a = line.replace("'","").strip().split('\t')
      b = a[3].split(', ')
      muts = []
      for x in b:                               # Looping through mutation list
         if x.find('/B') != -1:                 # = mutation + insertion, will be code 2 in vector, so use Z here 
            x = x[:-3]
            x = x+'Z'      
         elif x.find('/^') != -1: 
            x = x[:-3]                          # = mutation + deletion, will be code 5 in vector, so use ! here
            x = x+'!'
         print a[0]," var = ",x                 # since the length of the insertion is not recorded anyway and we are making a binary vector, the outcome is the same, so dropping it is ok
         if not positions.has_key(int(x[1:-1])):
            positions[int(x[1:-1])] = x[:-1]    # make to have the form { 19:'T19'}
         muts.append(x)

      data[a[0]] = muts
      misc[a[0]] = (a[1],a[2])  # length, png#
      print data


   print "POSITIONS*******\n",positions
   print "DATA $$$$$$$$$$$$\n",data
   pos_sorted = positions.keys()
   pos_sorted.sort()
   print "SORTED POSITIONS===============\n",pos_sorted
   head = ""
   for x in pos_sorted:
      head = head + positions[x]+'\t'
   head.strip()

   f2.write('Accession\tLength\tPNG_count\t'+head+'\n')
   for acc in data:
      f2.write(acc+'\t'+misc[acc][0]+"\t"+misc[acc][1]+"\t")
      for pos in pos_sorted:
         for m in data[acc]:
            if positions[pos] == m[:-1] : 
               if m[-1] == 'Z': f2.write('Z')   # mutation followed by insertion (was 2)
               elif m[-1] == 'B': f2.write('B') # insertion only (was 3)
               elif m[-1] == '^': f2.write('^') # Deletion only (was 4)
               elif m[-1] == '!': f2.write('!') # mutation followed by deletion (was 5)
               else: f2.write(m[-1].upper())    # mutation only (was 1)
               break
         else: f2.write('0')                    # No change from HxB2
         f2.write('\t')
      f2.write('\n')

   f1.close()
   f2.close()
   return


################
# processFiles #
################
def processFiles(seqfile):
      
   # Begin processing amino acid sequences to get mutations
   f2 = open(seqfile+"_gp120seq_mutations.txt",'w') # open output file
   f3 = open(seqfile+"_gp41seq_mutations.txt","w")
   for x in SeqIO.parse(seqfile,'fasta'):           # read seq file
      mutations = getMutations(x)                   # generate mutation list
      f2.write(x.id.replace('_','-')+"\t"+mutations[2]+"\t"+mutations[4]+"\t'"+mutations[0]+"'\n")   # write output seqid\tlength\t#pngs\tgp120_mutation_list
      f3.write(x.id.replace('_','-')+"\t"+mutations[3]+"\t"+mutations[5]+"\t'"+mutations[1]+"'\n")   # write output seqid\tlength\t#pngs\tgp41_mutation_list
   f2.close()
   f3.close()
   
   # Make 0/1 mutation vector for gp120
   print "/////////////////// Making vector for gp120 //////////////////"
   makeVector(seqfile+"_gp120seq_mutations.txt")

   # Make 0/1 mutation vector gp41
   print "//////////////////Making vector for GP41 /////////////////////"
   makeVector(seqfile+"_gp41seq_mutations.txt")



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
 