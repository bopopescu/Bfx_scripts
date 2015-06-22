#!/usr/bin/env python
# encoding: utf-8
"""
pepa_seq_prep.py
Generate input file for PEPA (Predictive Epitope Patch Analysis)
PEPA input format is seqid\tIC50\t"mutation list"
This script will take IC50 data file and raw sequence data as input

Created by Mark Evans on 2011-08-11.
Copyright (c) 2011 __Monogram Biosciences__. All rights reserved.

Revised 09.20.2011
Revised 02.24.2011
"""

import sys,getopt
import os, subprocess, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import AlignInfo
from Bio import AlignIO


#############
# Usage     #
#############
def usage():
   code =  "\n\n#############\nProgram indentifies mutations relative to HXB2 and generates a mutation list, \n"
   code += "combined with ic50 data to generate input file for surface patch analysis\n"
   code += "All antibody IC50 files must be in a subfolder to this script called ic50input\n"
   code += "Sequence data is protein sequence\n\n"
   code += "Usage: $ python pepa_seq_prep7.py [-p] gp41patch_filename [-q] gp120patch_filename [-s] seq_filename\n\n"
   code += "-p [patches41] name of file containing gp41 surface patch data\n"
   code += "-q [patches120] name of file containing gp120 surface patch data\n"   
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


##########################################################################################
# getMutations                                                                           #
# This method only looks for mutations within the constant regions.                      #
# There is no point trying to identify mutations in variable region because its variable #
##########################################################################################
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
   
   return (gp120_mut_sum, gp41_mut_sum, gp120_len, gp41_len)



##########################################################
# makeVector                                             #
# Takes mutation list output and converts to 0//1 vector #
##########################################################
def makeVector(mut_file):
   
   f1 = open(mut_file,'r')
   f2 = open(mut_file+".vector.out",'w')

   positions = {}
   data = {}

   for line in f1:
      a = line.replace("'","").strip().split('\t')
      b = a[2].split(', ')
      muts = []
      for x in b:
         if x.find('/') != -1: x = x[:-2]
         muts.append(x)

      data[a[0]] = muts
      print data
      for var in b:
         print a[0]," var = ",var
         if var.find('/') != -1: var = var[:-2]
         if not positions.has_key(var[1:-1]):
            positions[int(var[1:-1])] = var


   pos_sorted = positions.keys()
   pos_sorted.sort()
   head = ""
   for x in pos_sorted:
      head = head + positions[x][:-1]+'\t'
   head.strip()

   f2.write('Accession'+'\t'+head+'\n')
   for acc in data:
      f2.write(acc+'\t')
      for pos in pos_sorted:
         if positions[pos] in data[acc][:-1]: f2.write('1')
         else: f2.write('0')
         f2.write('\t')
      f2.write('\n')

   f1.close()
   f2.close()
   return



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
      header = "seqid\tX4_RLU\tR5_RLU\tic50_val\t"+"\t".join(p[x])+"\n"
      results += header
      for y in range(0,len(records)):                             # Loop through sequences
         d = records[y][1:-1]       #
         d.reverse()                # added this to put IC50 on the end for Mojgan.
         results += records[y][0]+"\t"+"\t".join(d)
        #results += "\t".join(records[y][0:-1]) #records[y][0]+"\t"+records[y][1]              # Write id and ic50 for sequence 
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
def processFiles(gp41patchfile,gp120patchfile,seqfile):
   datapath='ic50input/'
   path=""

   mut_set = {'gp120':{},'gp41':{}}

   # Begin processing amino acid sequences to get mutations
   f2 = open(seqfile+"_gp120seq_mutations.txt",'w') # open output file
   f3 = open(seqfile+"_gp41seq_mutations.txt","w")
   for x in SeqIO.parse(seqfile,'fasta'):           # read seq file
      mutations = getMutations(x)                   # generate mutation list
      f2.write(x.id.replace('_','-')+"\t"+mutations[2]+"\t'"+mutations[0]+"'\n")   # write output seqid\tseqlength\tgp120_mutation_list
      f3.write(x.id.replace('_','-')+"\t"+mutations[3]+"\t'"+mutations[1]+"'\n")   # write output seqid\tseqlength\tgp41_mutation_list
      mut_set['gp120'][x.id.replace('_','-')] = mutations[0]
      mut_set['gp41'][x.id.replace('_','-')] = mutations[1]
   f2.close()
   f3.close()
   
   # Make 0/1 mutation vector for gp120 ( not used for rest of process, just to have for future reference )
   makeVector(seqfile+"_gp120seq_mutations.txt")

   # Make 0/1 mutation vector gp41 ( not used for rest of process, just to have for future reference )
   makeVector(seqfile+"_gp41seq_mutations.txt")


   # Begin processing IC50 data for each antibody file in input folder
   for datafile in glob.glob (os.path.join(datapath,"*")):            # Read all filenames in directory
      a = datafile.split('/')[1].split('_')
      mab = a[0]
      ic50s = {}
      # Get IC50 data and put in dict with seqid as key
      f1 = open(datafile,'r')
      for line in f1:
         a = line.strip().split('\t')
         datas = a[1:]        # making the following changes allows multiple columns in addition to IC50 to be used
         ic50s[a[0]] = datas  # added reverse to make IC50 the last column before the 0/1 vector stuff
      f1.close()
      print ic50s
      # Already generated mutations, now merge them with each IC50 files
      f2 = open(mab+"_pepa_gp120seq_input.txt",'w')           # open output file
      f3 = open(mab+"_pepa_gp41seq_input.txt","w")
      seqids = mut_set['gp120'].keys()                     # get list of sequence ids from mutation dict
      seqids.sort()
      print seqids
      for x in seqids:                                     # go through each seqid
         gp41mutations = mut_set['gp41'][x]                # get mutations associated with this id
         gp120mutations = mut_set['gp120'][x]
         print gp120mutations
         data = ''
         if ic50s.has_key(x.replace('_','-')): # find IC50 that matches seqid
            data = ic50s[x.replace('_','-')]
# do not use             data[0] = str(data[0]).strip()         # clean up IC50 value
#                        if data[0] == 'NEG' or data[0]=='-1.00': data[0] = '0'
#                        else: data[0] = '1'
            f2.write(x.replace('_','-')+"\t"+"\t".join(data)+"\t'"+gp120mutations+"'\n")   # write output seqid\tIC50\gp120_mutation_list for knime input
            f3.write(x.replace('_','-')+"\t"+"\t".join(data)+"\t'"+gp41mutations+"'\n")   # write output seqid\tIC50\gp41_mutation_list for knime input
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
   gp41patchfile=''
   gp120patchfile=''
   seqfile=''

   #Read command line arguments
   try:
       opts, args = getopt.getopt(argv,"p:q:s:h:",["gp41patchfilename,gp120patchfilename,seq_filename,help"])
   except getopt.GetoptError:
       usage()
       sys.exit(2)
   for opt, arg in opts:
      if opt in ("-h","--help"):
         usage()
         sys.exit(2)

      elif opt == "-p": gp41patchfile = arg
      elif opt =="-q": gp120patchfile = arg
      elif opt =="-s": seqfile = arg
       
   if gp41patchfile !='' and gp120patchfile !='' and seqfile !='':
       processFiles(gp41patchfile,gp120patchfile,seqfile)
   else:
       print "There was an error\n"
       print opts
       print sys.argv
       usage()
       sys.exit()


if __name__ == '__main__':
   main(sys.argv[1:])

