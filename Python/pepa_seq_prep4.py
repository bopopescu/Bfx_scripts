#!/usr/bin/env python
# encoding: utf-8
"""
pepa_seq_prep.py
Generate input file for PEPA (Predictive Epitope Patch Analysis)
PEPA input format is seqid\tIC50\t"mutation list"
This script will take IC50 data file and raw sequence data as input

Created by Mark Evans on 2011-08-11.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.

Revised 09.20.2011
"""

import sys,getopt
import os, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import AlignInfo
from Bio import AlignIO


#############
# Usage     #
#############
def usage():
   code =  "\n\n#############\nProgram indentifies mutations relative to HXB2 and generates a mutation list, \n"
   code += "combined with ic50 data to generate input file for surface patch analysis\n\n"
   code += "Usage: $ python pepa_seq_prep3.py [-i] ic50_filename [-s] seq_filename [-p] patch_filename [-m] mab_name\n\n"
   code += "-i [ic50] name of file containing ic50 data in format 'seqid\\tIC50'\n"
   code += "-s [seq] nucleotide sequence string corresponding to gp160 region\n"
   code += "-p41 [patches41] name of file containing gp41 surface patch data\n"
   code += "-p120 [patches120] name of file containing gp120 surface patch data\n"   
   code += "-m [mab] name of monoclonal antibody used to generate the ic50 data\n"
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
      if '#' not in x:
         if a[0] == 'C1.1': cregions['C1'] = [int(a[17])-1,int(a[18])-1]    # alreday corrected for python zero-based position count
         elif a[0] == 'C3': cregions['C3'] = [int(a[17])-1,int(a[18])-1]
         elif a[0] == 'C4': cregions['C4'] = [int(a[17])-1,int(a[18])-1]
         elif a[0] == 'C5': cregions['C5'] = [int(a[17])-1,int(a[18])-1]
         elif a[0] == 'C6': cregions['C6'] = [int(a[17])-1,int(a[18])-1]
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
def processAlignment(local_path,start):
   ref_aln = ''
   q_aln = ''
   mut_pos = {}
   mut_sum =''
   
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
      
      if insert == 'yes' and ref_aln[pos] != '-' and q_aln[pos] != '-': insert = 'no'   # Turn off insert flag
      
      if ref_aln[pos].upper() != q_aln[pos].upper() and ref_aln[pos] != '-' and q_aln[pos] != '-':  # Record simple mutation
         mut_pos[pos] = ref_aln[pos]+str(pos+start+2)+q_aln[pos]
      elif ref_aln[pos] == '-' and insert == 'no' and pos != 0:                   # If query seq has an insertion....
         if mut_pos.has_key(pos-1):                                               # check that a mutation did not occur just before insertion
            mut_pos[pos-1] = mut_pos[pos-1]+"/B"                                  # if it did, add insertion symbol to it
         else:
            mut_pos[pos-1] = ref_aln[pos-1]+str(pos+start+1)+'B'                  # Step back one position in ref and use B to indicate insertion occurs after this position
            insert = 'yes'
      elif q_aln[pos] == '-' and insert == 'no' and pos != 0:                     # If query has a deletion..
         if mut_pos.has_key(pos-1):                                               # check that mutation did not occur just before deletion
            mut_pos[pos-1] = mut_pos[pos-1]+"/^"                                  # if it did, add deletion symbol to it
         else:
            mut_pos[pos-1] = ref_aln[pos-1]+str(pos+1+start)+"^"                  # Use ^ to indicate a deletion occurs after this position
            insert = 'yes'
   
   # Dump mut_pos contents into string to be written as part of the output
   pos_keys = mut_pos.keys()
   pos_keys.sort()
   for pos in pos_keys:
      mut_sum += mut_pos[pos]+", "
   
   return mut_sum


#################
# getMutations  #
#################
def getMutations(x):
   #######################################
   # HXB2 references sequences from LANL #
   #######################################
   gp160 = "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKRAVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL"
   gp120 = "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR"
   gp41  = "AVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL"
   gp120c1 = "KLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLK"
   gp120c3 = "CNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEIN"
   gp120c4 = "CNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFY"
   gp120c5 = "CRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGG"
   gp120c6 = "PGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR"
   const_seqs = {'C1':gp120c1,'C3':gp120c3,'C4':gp120c4,'C5':gp120c5,'C6':gp120c6}
   coords = {'C1':[32,129],'C3':[195,294],'C4':[330,383],'C5':[417,458],'C6':[469,510]} # HXB2 gp160 constant region coords, adjusted to be zero-based
   r = ['C1','C3','C4','C5','C6']
   
   #################################
   # Begin process for gp120       #
   #################################
   tempfile = open("hmm.seq","w")
   tempfile.write(">"+x.id+"\n"+x.seq.tostring()+"\n")
   tempfile.close()
   local_path = os.getcwd()+"/"
   hmmscan_bin = "/usr/local/bin/hmmscan"
   hmmresult_filename = doHMMScan("hmm.seq",local_path,hmmscan_bin)
   target_cregions = processHMMresult(hmmresult_filename,local_path)
   gp120_mut_sum=''
   
   for c in r:
      # Create temp file for Clustal input
      start = target_cregions[c][0]
      stop  = target_cregions[c][1]
      content = ">Ref\n"+const_seqs[c]+"\n>Query\n"+x.seq.tostring().strip('*')[start:stop]+"\n"
      tempfile = open('foo.seq','w')
      tempfile.write(content)
      tempfile.close()
      
      gp120_mut_sum += processAlignment(local_path,start)
      
   gp120_mut_sum = gp120_mut_sum.strip()  # remove trailing whitespace
   gp120_mut_sum = gp120_mut_sum[0:-1]    # remove trailing ,
   
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
   
   gp41_mut_sum = processAlignment(local_path, 0)
   gp41_mut_sum = gp41_mut_sum.strip()  # remove trailing whitespace
   gp41_mut_sum = gp41_mut_sum[0:-1]    # remove trailing ,
   
   return (gp120_mut_sum,gp41_mut_sum)


#######################
# makePatchVectors    #
#######################
def makePatchVectors(p,infile,mab,region):
   inputfile = open(infile,"r")
   records = []
   results = ''
   os.mkdir(region+"_"+mab)
   local_path = os.getcwd()
   local_path += "/"+region+"_"+mab+"/"
   
   for line in inputfile:
      records.append(line.split('\t'))
   
   for x in range(0,len(p)):                                      # Loop through patches
      f = open(local_path+"patch"+str(x)+".txt","w")
      header = "seqid\tTropism\tR5_RLU\tX4_RLU\tCurve_tension\tInit_Dil\tInhib\tPCT_Inhib\tIC50_vec\tic50_val\t"+"\t".join(p[x])+"\n"
      results += header
      for y in range(0,len(records)):                             # Loop through sequences
         results += "\t".join(records[y][0:-1]) #records[y][0]+"\t"+records[y][1]              # Write id and ic50 for sequence 
         for mut in p[x]:                                         # Loop through each patch position
            if records[y][-1].find(mut) != -1:  results += "\t1"   # Check if sequence has that position
            else: results += "\t0"                                # or does not have that position
         results += "\n"
      f.write(results)
      f.close()
      results = ''
      
      
   
   



################
# processFiles #
################
def processFiles(ic50file,seqfile,gp41patchfile,gp120patchfile,mab):
   
   # Get IC50 data and put in dict with seqid as key
   ic50s = {}
   f1 = open(ic50file,'r')
   for line in f1:
      a = line.strip().split('\t')
      datas = a[1:]        # making the following changes allows multiple columns in addition to IC50 to be used
      ic50s[a[0]] = datas
   f1.close()
   
   # Begin processing amino acid sequences to get mutations and merge with ic50 data
   f2 = open(mab+"_pepa_gp120seq_input.txt",'w')           # open output file
   f3 = open(mab+"_pepa_gp41seq_input.txt","w")
   for x in SeqIO.parse(seqfile,'fasta'):          # read seq file
      mutations = getMutations(x)                  # generate mutation list
      data = ''
      if ic50s.has_key(x.id.replace('_','-')): # find IC50 that matches seqid
         data = ic50s[x.id.replace('_','-')]
#         data[0] = str(data[0]).strip()         # clean up IC50 value
#         if data[0] == 'NEG' or data[0]=='-1.00': data[0] = '0'
 #        else: data[0] = '1'
         f2.write(x.id.replace('_','-')+"\t"+"\t".join(data)+"\t'"+mutations[0]+"'\n")   # write output seqid\tIC50\gp120_mutation_list for knime input
         f3.write(x.id.replace('_','-')+"\t"+"\t".join(data)+"\t'"+mutations[1]+"'\n")   # write output seqid\tIC50\gp41_mutation_list for knime input
   f2.close()
   f3.close()
   
   # Generate patch vector files for analysis
   p41 = []
   p120 =[]
   f4 = open(gp41patchfile,"r")
   f5 = open(gp120patchfile,"r")
   for x in f4:
      p41.append(x.strip().split(','))
   for x in f5:
      p120.append(x.strip().split(','))
   f4.close()
   f5.close()
   makePatchVectors(p41,mab+"_pepa_gp41seq_input.txt",mab,'gp41')
   makePatchVectors(p120,mab+"_pepa_gp120seq_input.txt",mab,'gp120')



###########
# Main    #
###########
def main(argv):
   ic50file=''
   seqfile=''
   gp41patchfile=''
   gp120patchfile=''
   mab=''
   
   # Read command line arguments
   try:
      opts, args = getopt.getopt(argv,"i:s:p:q:h:m:",["ic50,seq,patches41,patches120,help,mab"])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt in ("-h","--help"):
         usage()
         sys.exit()
      elif opt in ("-i","--ic50"): ic50file = arg
      elif opt in ("-s","--seq"): seqfile = arg
      elif opt in ("-p","--patches41"): gp41patchfile = arg
      elif opt in ("-q","--patches120"): gp120patchfile = arg      
      elif opt in ("-m","--mab"): mab = arg
   
   if ic50file !='' and seqfile !='' and gp41patchfile !='' and gp120patchfile !='':
      processFiles(ic50file,seqfile,gp41patchfile,gp120patchfile,mab)
   else:
      print opts
      print sys.argv
      usage()
      sys.exit()



if __name__ == '__main__':
   main(sys.argv[1:])

