#!/usr/bin/env python
# encoding: utf-8
"""
calcmut.py

Created by Mark Evans on 2011-04-29.
Copyright (c) 2011 __Monogram_Biosciences__. All rights reserved.
"""

import sys
import os
from decimal import *

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
        'T69':'',
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
PI_RAMS = {'L23A':'','L23C':'','L23D':'','L23E':'','L23F':'','L23G':'','L23H':'','L23I':'','L23K':'','L23M':'','L23N':'','L23P':'','L23Q':'','L23R':'','L23S':'','L23T':'','L23V':'','L23W':'','L23Y':'',
           'L24A':'','L24C':'','L24D':'','L24E':'','L24F':'','L24G':'','L24H':'','L24I':'','L24K':'','L24M':'','L24N':'','L24P':'','L24Q':'','L24R':'','L24S':'','L24T':'','L24V':'','L24W':'','L24Y':'',
           'D30A':'','D30C':'','D30E':'','D30F':'','D30G':'','D30H':'','D30I':'','D30K':'','D30L':'','D30M':'','D30N':'','D30P':'','D30Q':'','D30R':'','D30S':'','D30T':'','D30V':'','D30W':'','D30Y':'',
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



def main():
   f = open('seqtest','r')
   f2 = open('processed_seq.txt','w')
   pri_seen = {}     # Record of all protease mutations that have been seen, whether they are RAMs or not
   rt_seen = {}      # Record of all RT mutations that have been seen, whether they are RAMs or not
   getcontext().prec=2     # Sets precision for decimal to 2 digits
   
   # Begin processing sequence file (Oracle export for building HIVWH - contains pri_summary and rt_summary fields)
   for line in f:
      tams_sum=''
      nams_sum=''
      nnrti_sum=''
      pi_sum=''
      wt=''
      pri_mut_score = Decimal(0)
      rt_mut_score = Decimal(0)
      #nnrti_mut_score = 0
      #pi_mut_score = 0
      
      a = line.split('"')     # Since summary fields are exported as being enclosed by quotes, we will split on this first
      pri_sum = a[1].replace(' ','').split(',')
      rt_sum = a[3].replace(' ','').split(',')
      
      # Process Protease sequences
      for mut in pri_sum:
         if pri_seen.has_key(mut): pri_seen[mut] = pri_seen[mut] + 1  # Check if mut seen before and increment counter.
         else: pri_seen[mut] = 1
         # Check for single mutants first
         if mut.find('/') == -1:
            if PI_RAMS.has_key(mut):
               pri_mut_score = pri_mut_score + 1
               if pi_sum != '': pi_sum = pi_sum+", "+mut
               else: pi_sum = mut
         # Have multiple mutants
         else:
            mut2 = mut.split('/')   # ['I50A','C','H']
            mut_count = len(mut2)
            if mut2[0][0]!=mut2[0][-1]:
               # check first mutation in list
               if pri_seen.has_key(mut2[0]): pri_seen[mut2[0]] = pri_seen[mut2[0]] + 1  # Check if mut seen before and increment counter.
               else: pri_seen[mut2[0]] = 1
               if PI_RAMS.has_key(mut2[0]):
                  pri_mut_score = pri_mut_score + Decimal(1)/Decimal(mut_count)
                  if pi_sum != '': pi_sum = pi_sum+", "+mut2[0]
                  else: pi_sum = mut2[0]
            # Generate the rest of mutations in list and check them.
            for i in range(1,len(mut2)):
               m = mut2[0][:-1]+mut2[i] # automatically create I50C, I50H etc
               if pri_seen.has_key(m): pri_seen[m] = pri_seen[m] + 1  # Check if mut seen before and increment counter.
               else: pri_seen[m] = 1 
               if PI_RAMS.has_key(m):
                  pri_mut_score = pri_mut_score + Decimal(1)/Decimal(mut_count)
                  if pi_sum != '': pi_sum = pi_sum+", "+m
                  else: pi_sum = m
                  
            
      # Process RT sequences
      for mut in rt_sum:
         print "Looking at RT mut ",mut
         if rt_seen.has_key(mut): rt_seen[mut] = rt_seen[mut] + 1
         else: rt_seen[mut] = 1
         # Check for single mutants first
         if mut.find('/') == -1:
            if TAMS.has_key(mut):
               rt_mut_score = rt_mut_score + 1
               if tams_sum != '': tams_sum = tams_sum+", "+mut
               else: tams_sum = mut
            elif NAMS.has_key(mut):
               rt_mut_score = rt_mut_score + 1
               if nams_sum != '': nams_sum = nams_sum+", "+mut
               else: nams_sum = mut
            elif NNRTI_RAMS.has_key(mut):
               rt_mut_score = rt_mut_score + 1
               if nnrti_sum != '': nnrti_sum = nnrti_sum+", "+mut
               else: nnrti_sum = mut
         # Have multiple mutants
         else:
            mut2 = mut.split('/')   # ['I50A','C','H']
            print "   found ",mut2
            mut_count = len(mut2)
            if mut2[0][0]!=mut2[0][-1]:  # Make sure that first element is not wt (e.g. I50I)
               # check first mutation in list
               if rt_seen.has_key(mut2[0]): rt_seen[mut2[0]] = rt_seen[mut2[0]] + 1  # Check if mut seen before and increment counter.
               else: rt_seen[mut2[0]] = 1
               if TAMS.has_key(mut2[0]):
                  rt_mut_score = rt_mut_score + Decimal(1)/Decimal(mut_count)
                  if tams_sum != '': tams_sum = tams_sum+", "+mut2[0]
                  else: tams_sum = mut2[0]
               elif NAMS.has_key(mut2[0]):
                  rt_mut_score = rt_mut_score + Decimal(1)/Decimal(mut_count)
                  if nams_sum != '': nams_sum = nams_sum+", "+mut2[0]
                  else: nams_sum = mut2[0]
               elif NNRTI_RAMS.has_key(mut2[0]):
                     rt_mut_score = rt_mut_score + Decimal(1)/Decimal(mut_count)
                     if nnrti_sum != '': nnrti_sum = nnrti_sum+", "+mut2[0]
                     else: nnrti_sum = mut2[0]
            # Generate the rest of mutations in list and check them.
            for i in range(1,len(mut2)):
               m = mut2[0][:-1]+mut2[i] # automatically create I50C, I50H etc
               print "     looking at alt mut ",m
               if rt_seen.has_key(m): rt_seen[m] = rt_seen[m] + 1  # Check if mut seen before and increment counter.
               else: rt_seen[m] = 1 
               if TAMS.has_key(m):
                  rt_mut_score = rt_mut_score + Decimal(1)/Decimal(mut_count)
                  if tams_sum != '': tams_sum = tams_sum+", "+m
                  else: tams_sum = m
               elif NAMS.has_key(m):
                  print "          I see a nam: ",m
                  rt_mut_score = rt_mut_score + Decimal(1)/Decimal(mut_count)
                  if nams_sum != '': nams_sum = nams_sum+", "+m
                  else: nams_sum = m
               elif NNRTI_RAMS.has_key(m):
                     rt_mut_score = rt_mut_score + Decimal(1)/Decimal(mut_count)
                     if nnrti_sum != '': nnrti_sum = nnrti_sum+", "+m
                     else: nnrti_sum = m
                     
      if rt_mut_score + pri_mut_score == 0: wt='y'
      else: wt = 'n'
      a[len(a)-1] = a[len(a)-1].rstrip()
      f2.write("\t".join(a)+"\t"+wt+"\t"+tams_sum+"\t"+nams_sum+"\t"+nnrti_sum+"\t"+pi_sum+"\t"+str(pri_mut_score)+"\t"+str(rt_mut_score)+"\n")
   f2.close()
   f.close()
   
   
if __name__ == '__main__':
   main()

