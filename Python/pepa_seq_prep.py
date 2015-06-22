#!/usr/bin/env python
# encoding: utf-8
"""
pepa_seq_prep.py
Generate input file for PEPA (Predictive Epitope Patch Analysis)
PEPA input format is seqid\tIC50\t"mutation list"
This script will take IC50 data file and raw sequence data as input

Created by Mark Evans on 2011-08-11.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import AlignInfo
from Bio import AlignIO

# From KNIME, Array = IC50, aa1, aa2, etc
##########################################

#setwd("/Users/Bali/Documents/Monogram/R_space/node1");
#R<-R
#mylm=lm(IC50 ~ .,data=R)
#R=mylm

#R<-cbind(RDATA, predict(RMODEL, RDATA));
##R<-predict(RMODEL, RDATA);
#print(RMODEL)

#R[,length(R)+1] = 10^R$"IC50" >=50
#pred=R[,length(R)]
#print(names(R))
#R[,length(R)+1] = 10^R$"predict(RMODEL, RDATA)" >=50

#a=cbind(R[,1],R[,(length(R)-2):length(R)])
#labs=c("IC50","Predicted IC50","Resistant","Predicted Resistant")

#print(names(a))
#print(labs)
#names(a)=labs

#R=a

#setwd("/Users/Bali/Documents/Monogram/R_space/node1");
#R<-R
#df=R
#   cv.lm<-lm(df[,1] ~ df[,2],data=df)
#      psum<-summary(cv.lm)
#      rsq<-c(psum$r.squared)


#      f.stat<-psum$fstatistic
#      if(psum$r.squared == 0){
#         p.value = 1
#      }else{
#         p.value <- 1-pf(f.stat["value"],f.stat["numdf"],f.stat["dendf"])
#      }



#   mse<-sum((df[,1] - df[,2])^2)/length(df[,1])
#   mmse<-median((df[,1] - df[,2])^2)

#   results <- c("$pnum ", "$patch", p.value,psum$r.squared,mse)

#R=t(results)









def getMutations(x):
   # HXB2 references sequences from LANL
   gp160 = "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKRAVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL"
   gp120 = "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR"
   gp41  = "AVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL"

   # Create temp file for Clustal input
   content = ">Ref\n"+gp160+"\n>Query\n"+x.seq.tostring().strip('*')+"\n"
   tempfile = open('foo.seq','w')
   tempfile.write(content)
   tempfile.close()
   local_path="/Users/Bali/Documents/Monogram/Terri/patch_analysis/multiclade_project/dev/"
   
   # Call Clustalw2
   subprocess.call(["clustalw2","-align","-infile="+local_path+"foo.seq", "-outfile="+local_path+"foo.aln"])
   
   ref_aln = ''
   q_aln = ''
   mut_pos = {}
   mut_sum=''
   
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
   for x in range(0,len(ref_aln)):
      
      if insert == 'yes' and ref_aln[x] != '-' and q_aln[x] != '-': insert = 'no'   # Turn off insert flag
      
      if ref_aln[x].upper() != q_aln[x].upper() and ref_aln[x] != '-' and q_aln[x] != '-':  # Record simple mutation
         mut_pos[x] = ref_aln[x]+str(x)+q_aln[x]
      elif ref_aln[x] == '-' and insert == 'no' and x != 0:         # If query seq has an insertion....
         if mut_pos.has_key(x-1):                                   # check that a mutation did not occur just before insertion
            mut_pos[x-1] = mut_pos[x-1]+"/B"                        # if it did, add insertion symbol to it
         else:
            mut_pos[x-1] = ref_aln[x-1]+str(x-1)+'B'                # Step back one position in ref and use B to indicate insertion occurs after this position
            insert = 'yes'
      elif q_aln[x] == '-' and insert == 'no' and x != 0:           # If query has a deletion..
         if mut_pos.has_key(x-1):                                   # check that mutation did not occur just before deletion
            mut_pos[x-1] = mut_pos[x-1]+"/^"                        # if it did, add deletion symbol to it
         else:
            mut_pos[x-1] = ref_aln[x-1]+str(x-1)+"^"                # Use ^ to indicate a deletion occurs after this position
            insert = 'yes'
   
   # Dump mut_pos contents into string to be written as part of the output
   pos_keys = mut_pos.keys()
   pos_keys.sort()
   for pos in pos_keys:
      mut_sum += mut_pos[pos]+", "
   mut_sum = mut_sum.strip()  # remove trailing whitespace
   mut_sum = mut_sum[0:-1]    # remove trailing ,
   
   return mut_sum
   
   
def main():
   infile1 = raw_input('What is the name of the file containing IC50 data (seqid\\tIC50 format): ')
   infile2 = raw_input("what is the name of the matching amino acid seq file: ")
   
   # Get IC50 data and put in dict with seqid as key
   ic50s = {}
   f1 = open(infile1,'r')
   for line in f1:
      a = line.split('\t')
      ic50s[a[0]] = a[1]
   f1.close()
   
   # Begin processing amino acid sequences
   f2 = open("pepa_seq_input.txt",'w')           # open output file
   for x in SeqIO.parse(infile2,'fasta'):        # read seq file
      mut_sum = getMutations(x)                  # generate mutation list
      ic50 = ''
      if ic50s.has_key(x.id.replace('_','-')): ic50 = str(ic50s[x.id.replace('_','-')]).strip()   # find IC50 that matches seqid
      if ic50 == 'NEG': ic50 = '-1'
      f2.write(x.id.replace('_','-')+"\t"+ic50+"\t'"+mut_sum+"'\n")   # write output seqid\tIC50\mutation_list for knime input
   f2.close()


if __name__ == '__main__':
   main()

