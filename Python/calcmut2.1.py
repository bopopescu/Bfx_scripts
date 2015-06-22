#!/usr/bin/env python
# encoding: utf-8
"""
calcmut.py

Created by Mark Evans on 2011-04-29.
Revised 2011-05-06
Copyright (c) 2011 __Monogram_Biosciences__. All rights reserved.
"""

import sys
import os
from decimal import *
from time import localtime, strftime

#aa = {A':'',C':'',D':'',E':'',F':'',G':'',H':'',I':'',K':'',L':'',M':'',N':'',P':'',Q':'',R':'',S':'',T':'',V':'',W':'',Y':''}

# NRTI TAMS: M41L,D67N,K70R,L210W,T215F,T215Y,K219
TAMS = {'M41L':'',
        'D67N':'',
        'K70R':'',
        'L210W':'',
        'T215F':'','T215Y':'',
        'K219A':'','K219C':'','K219D':'','K219E':'','K219F':'','K219G':'','K219H':'','K219I':'','K219L':'','K219M':'','K219N':'','K219P':'','K219Q':'','K219R':'','K219S':'','K219T':'','K219V':'','K219W':'','K219Y':'',
        }
        
# NRTI NAMS: K65R,T69,K70E,L74,V75A,V75M,V75S,V75T,Y115F,Q151M,M184
NAMS = {'K65R':'',
        'T69A':'','T69C':'','T69D':'','T69E':'','T69F':'','T69G':'','T69H':'','T69I':'','T69K':'','T69L':'','T69M':'','T69N':'','T69P':'','T69Q':'','T69R':'','T69S':'','T69V':'','T69W':'','T69Y':'','T69^':'',
        'K70E':'',
        'L74A':'','L74C':'','L74D':'','L74E':'','L74F':'','L74G':'','L74H':'','L74I':'','L74K':'','L74M':'','L74N':'','L74P':'','L74Q':'','L74R':'','L74S':'','L74T':'','L74V':'','L74W':'','L74Y':'',
        'V75A':'','V75M':'','V75S':'','V75T':'',
        'Y115F':'',
        'Q151M':'',
        'M184A':'','M184C':'','M184D':'','M184E':'','M184F':'','M184G':'','M184H':'','M184I':'','M184K':'','M184L':'','M184N':'','M184P':'','M184Q':'','M184R':'','M184S':'','M184T':'','M184V':'','M184W':'','M184Y':''
        }
        
# NRTI RAMS (TAMS + NAMS): M41L,K65R,D67N,T69,K70R,K70E,L74,V75A,V75M,V75S,V75T,Y115F,Q151M,M184,L210W,T215F,T215Y,K219
NRTI_RAMS = {'M41L':'',
             'K65R':'',
             'D67N':'',
             'T69':'',
             'K70R':'','K70E':'',
             'L74A':'','L74C':'','L74D':'','L74E':'','L74F':'','L74G':'','L74H':'','L74I':'','L74K':'','L74M':'','L74N':'','L74P':'','L74Q':'','L74R':'','L74S':'','L74T':'','L74V':'','L74W':'','L74Y':'',
             'V75A':'','V75M':'','V75S':'','V75T':'',
             'L210W':'',
             'Y115F':'',
             'Q151M':'',
             'M184A':'','M184C':'','M184D':'','M184E':'','M184F':'','M184G':'','M184H':'','M184I':'','M184K':'','M184L':'','M184N':'','M184P':'','M184Q':'','M184R':'','M184S':'','M184T':'','M184V':'','M184W':'','M184Y':'',
             'T215F':'','T215Y':'',
             'K219A':'','K219C':'','K219D':'','K219E':'','K219F':'','K219G':'','K219H':'','K219I':'','K219L':'','K219M':'','K219N':'','K219P':'','K219Q':'','K219R':'','K219S':'','K219T':'','K219V':'','K219W':'','K219Y':'',
             }
             
# NNRTI RAMS: A98G, L100I, K101E, K101P, K103N, K103S, V106A, V106M, Y181, Y188, G190, P225, F227, M230L, P236L
NNRTI_RAMS = {'A98G':'', 
              'L100I':'', 
              'K101E':'','K101P':'', 
              'K103N':'','K103S':'', 
              'V106A':'','V106M':'', 
              'Y181A':'','Y181C':'','Y181D':'','Y181E':'','Y181F':'','Y181G':'','Y181H':'','Y181I':'','Y181K':'','Y181L':'','Y181M':'','Y181N':'','Y181P':'','Y181Q':'','Y181R':'','Y181S':'','Y181T':'','Y181V':'','Y181W':'',
              'Y188A':'','Y188C':'','Y188D':'','Y188E':'','Y188F':'','Y188G':'','Y188H':'','Y188I':'','Y188K':'','Y188L':'','Y188M':'','Y188N':'','Y188P':'','Y188Q':'','Y188R':'','Y188S':'','Y188T':'','Y188V':'','Y188W':'',
              'G190A':'','G190C':'','G190D':'','G190E':'','G190F':'','G190H':'','G190I':'','G190K':'','G190L':'','G190M':'','G190N':'','G190P':'','G190Q':'','G190R':'','G190S':'','G190T':'','G190V':'','G190W':'','G190Y':'',
              'P225A':'','P225C':'','P225D':'','P225E':'','P225F':'','P225G':'','P225H':'','P225I':'','P225K':'','P225L':'','P225M':'','P225N':'','P225Q':'','P225R':'','P225S':'','P225T':'','P225V':'','P225W':'','P225Y':'',
              'F227A':'','F227C':'','F227D':'','F227E':'','F227G':'','F227H':'','F227I':'','F227K':'','F227L':'','F227M':'','F227N':'','F227P':'','F227Q':'','F227R':'','F227S':'','F227T':'','F227V':'','F227W':'','F227Y':'',
              'M230L':'', 
              'P236L':''
              }
              
# PI RAMS: L23,L24,D30,V32,M46,I47,G48,I50,I54,V82A,V82F,V82S,V82T,V82C,V82G,V82L,V82M,I84,N88,L90
PI_RAMS = {'D30A':'','D30C':'','D30E':'','D30F':'','D30G':'','D30H':'','D30I':'','D30K':'','D30L':'','D30M':'','D30N':'','D30P':'','D30Q':'','D30R':'','D30S':'','D30T':'','D30V':'','D30W':'','D30Y':'',
           'V32A':'','V32C':'','V32D':'','V32E':'','V32F':'','V32G':'','V32H':'','V32I':'','V32K':'','V32L':'','V32M':'','V32N':'','V32P':'','V32Q':'','V32R':'','V32S':'','V32T':'','V32W':'','V32Y':'',
           'M46A':'','M46C':'','M46D':'','M46E':'','M46F':'','M46G':'','M46H':'','M46I':'','M46K':'','M46L':'','M46N':'','M46P':'','M46Q':'','M46R':'','M46S':'','M46T':'','M46V':'','M46W':'','M46Y':'',
           'I47A':'','I47C':'','I47D':'','I47E':'','I47F':'','I47G':'','I47H':'','I47K':'','I47L':'','I47M':'','I47N':'','I47P':'','I47Q':'','I47R':'','I47S':'','I47T':'','I47V':'','I47W':'','I47Y':'',
           'G48A':'','G48C':'','G48D':'','G48E':'','G48F':'','G48H':'','G48I':'','G48K':'','G48L':'','G48M':'','G48N':'','G48P':'','G48Q':'','G48R':'','G48S':'','G48T':'','G48V':'','G48W':'','G48Y':'',
           'I50A':'','I50C':'','I50D':'','I50E':'','I50F':'','I50G':'','I50H':'','I50K':'','I50L':'','I50M':'','I50N':'','I50P':'','I50Q':'','I50R':'','I50S':'','I50T':'','I50V':'','I50W':'','I50Y':'',
           'I54A':'','I54C':'','I54D':'','I54E':'','I54F':'','I54G':'','I54H':'','I54K':'','I54L':'','I54M':'','I54N':'','I54P':'','I54Q':'','I54R':'','I54S':'','I54T':'','I54V':'','I54W':'','I54Y':'',
           'V82A':'','V82F':'','V82S':'','V82T':'','V82C':'','V82G':'','V82L':'','V82M':'',
           'I84A':'','I84C':'','I84D':'','I84E':'','I84F':'','I84G':'','I84H':'','I84K':'','I84L':'','I84M':'','I84N':'','I84P':'','I84Q':'','I84R':'','I84S':'','I84T':'','I84V':'','I84W':'','I84Y':'',
           'N88A':'','N88C':'','N88D':'','N88E':'','N88F':'','N88G':'','N88H':'','N88I':'','N88K':'','N88L':'','N88M':'','N88P':'','N88Q':'','N88R':'','N88S':'','N88T':'','N88V':'','N88W':'','N88Y':'',
           'L90A':'','L90C':'','L90D':'','L90E':'','L90F':'','L90G':'','L90H':'','L90I':'','L90K':'','L90M':'','L90N':'','L90P':'','L90Q':'','L90R':'','L90S':'','L90T':'','L90V':'','L90W':'','L90Y':'',
           }

# INT RAMS: Y143R,Y143H,Y143C,Q148H,Q148K,Q148R,N155H
INT_RAMS = {'E92Q':'',
            'F121C':'','F121N':'','F121Y':'',
            'Y143R':'','Y143H':'','Y143C':'',
            'Q148H':'','Q148K':'','Q148R':'',
            'N155H':'','N155S':''
            }

# RT2-specific RAMS only, not including overlap with RT
RT2_RAMS = {'N348I':'',
            'T369I':'','T369V':''
            }

def checkRT_RAMs(mutant):
   global tams_sum, nams_sum, nnrti_sum
   score = 0
   if TAMS.has_key(mutant):
      score = 1
      if tams_sum != '': tams_sum = tams_sum+", "+mutant
      else: tams_sum = mutant
   elif NAMS.has_key(mutant):
      score = 1
      if nams_sum != '': nams_sum = nams_sum+", "+mutant
      else: nams_sum = mutant
   elif NNRTI_RAMS.has_key(mutant):
      score = 1
      if nnrti_sum != '': nnrti_sum = nnrti_sum+", "+mutant
      else: nnrti_sum = mutant
   return score


def checkPI_RAMs(mutant):
   global pi_sum
   score = 0
   if PI_RAMS.has_key(mutant):
      score = 1
      if pi_sum != '': pi_sum = pi_sum+", "+mutant
      else: pi_sum = mutant
   return score


def checkINT_RAMs(mutant):
   global integrase_sum
   score = 0
   if INT_RAMS.has_key(mutant):
      score = 1
      if integrase_sum != '': integrase_sum = integrase_sum+", "+mutant
      else: integrase_sum = mutant
   return score


def checkRT2_RAMs(mutant):
   global rt2_sum
   score = 0
   if RT2_RAMS.has_key(mutant):
      score = 1
      if rt2_sum != '': rt2_sum = rt2_sum+", "+mutant
      else: rt2_sum = mutant
   return score


def seenPRI(mutant):
   global pri_seen
   if pri_seen.has_key(mutant): pri_seen[mutant] = pri_seen[mutant] + 1  # Check if mut seen before and increment counter.
   else: pri_seen[mutant] = 1
   return


def seenRT(mutant):
   global rt_seen
   if rt_seen.has_key(mutant): rt_seen[mutant] = rt_seen[mutant] + 1
   else: rt_seen[mutant] = 1
   return


def seenINT(mutant):
   global int_seen
   if int_seen.has_key(mutant): int_seen[mutant] = int_seen[mutant] + 1
   else: int_seen[mutant] = 1
   return


def seenRT2(mutant):
   global rt2_seen
   if rt2_seen.has_key(mutant): rt2_seen[mutant] = rt2_seen[mutant] + 1
   else: rt2_seen[mutant] = 1
   return


def processMuts(mutation_summary,mut_type):
   mut_score = 0
   if mut_type == 'pri':
      seen = seenPRI
      checkRAMS = checkPI_RAMs
   elif mut_type == 'rt':
      seen = seenRT
      checkRAMS = checkRT_RAMs
   elif mut_type == 'int':
      seen = seenINT
      checkRAMS = checkINT_RAMs
   elif mut_type == 'rt2':
      seen = seenRT2
      checkRAMS = checkRT2_RAMs
      
   # Process Protease sequences
   for mut in mutation_summary:
 #     if mut.find("*") != -1: print mut
 #     if mut.find("^") != -1: print mut
      seen(mut)
      # Check for single mutants first
      if mut.find('/') == -1: mut_score = mut_score + checkRAMS(mut)
      
      # Have multiple mutants
      else:
         mut2 = mut.split('/')   # ['I50A','C','H']
         mut_count = len(mut2)
         if mut2[0][0]!=mut2[0][-1]:
            # check first mutation in list
            seen(mut2[0])
            s = checkRAMS(mut2[0])
            if s == 1: mut_score = mut_score + Decimal(1)/Decimal(mut_count)
            
         # Generate the rest of mutations in list and check them.
         for i in range(1,len(mut2)):
            m = mut2[0][:-1]+mut2[i] # automatically create I50C, I50H etc
            seen(m)
            s = checkRAMS(m)
            if s == 1: mut_score = mut_score + Decimal(1)/Decimal(mut_count)
   return mut_score
   
   
def main():
   filename = raw_input("What is name of seq_export file? ")
   print "\nstarted at ",strftime("%H:%M:%S", localtime()),"\n\n"
   f = open(filename,'r')
   f2 = open(filename+'.processed.txt','w')
   getcontext().prec=2     # Sets precision for decimal to 2 digits
   global pri_seen, rt_seen, tams_sum, nams_sum, nnrti_sum, pi_sum, int_seen, integrase_sum,rt2_seen,rt2_sum
   pri_seen = {}     # Record of all protease mutations that have been seen, whether they are RAMs or not
   rt_seen = {}      # Record of all RT mutations that have been seen, whether they are RAMs or not
   int_seen = {}
   rt2_seen = {}
   
   # Begin processing sequence file (Oracle export for building HIVWH - contains pri_summary and rt_summary fields)
   for line in f:
      tams_sum=''
      nams_sum=''
      nnrti_sum=''
      pi_sum=''
      integrase_sum=''
      wt=''
      rt2_sum=''
      pri_mut_score = Decimal(0)
      rt_mut_score = Decimal(0)
      int_mut_score = Decimal(0)
      rt2_mut_score = Decimal(0)
      
   #   a = line.split('"')     # Since summary fields are exported as being enclosed by quotes, we will split on this first
      a = line.split('\t')
      if len(a[5]) > 0:
         #pri_sum = a[1].replace(' ','').split(',')
         pri_sum = a[5].replace(' ','').replace('"','').split(',')
         pri_mut_score = processMuts(pri_sum,'pri')
      if len(a[6]) > 0:
         #rt_sum = a[3].replace(' ','').split(',')
         rt_sum = a[6].replace(' ','').replace('"','').split(',')
         rt_mut_score = processMuts(rt_sum,'rt')
      if len(a[7]) > 0:
         int_sum = a[7].replace(' ','').replace('"','').split(',')
         int_mut_score = processMuts(int_sum,'int')
      if len(a[8]) > 0:
         rt2_summary = a[8].replace(' ','').replace('"','').split(',')
         rt2_mut_score = processMuts(rt2_summary,'rt2')
      
      print "rt_mut_Score=",rt_mut_score," pri_mut_Score=",pri_mut_score
      if rt_mut_score + pri_mut_score + int_mut_score== 0: wt='y'
      else: wt = 'n'
      a[len(a)-1] = a[len(a)-1].rstrip() # remove existing carriage return so original data can be merged with new data for printing
      a[5] = a[5].replace('"','')
      a[6] = a[6].replace('"','')
      a[7] = a[7].replace('"','')
      a[8] = a[8].replace('"','')
      f2.write("\t".join(a)+"\t"+wt+"\t"+tams_sum+"\t"+nams_sum+"\t"+nnrti_sum+"\t"+pi_sum+"\t"+integrase_sum+"\t"+rt2_sum+"\t"+str(pri_mut_score)+"\t"+str(rt_mut_score)+"\t"+str(int_mut_score)+"\t"+str(rt2_mut_score)+"\n")
   f2.close()
   f.close()
   f = open("mut_freq.txt","w")
   f.write("Mutation\tAbsFreq\tType\n")
   for x in pri_seen:
      f.write(x.replace('/','')+"\t"+x+"\t"+str(pri_seen[x])+"\tProtease\n")
   for y in rt_seen:
      f.write(y.replace('/','')+"\t"+y+"\t"+str(rt_seen[y])+"\tRT\n")
   for z in int_seen:
      f.write(z.replace('/','')+"\t"+z+"\t"+str(int_seen[z])+"\tInt\n")
   for w in rt2_seen:
      f.write(w.replace('/','')+"\t"+w+"\t"+str(rt2_seen[w])+"\tRT2\n")
   f.close()
   print "\nfinished at ",strftime("%H:%M:%S", localtime()),"\n\n"
   
if __name__ == '__main__':
   main()

