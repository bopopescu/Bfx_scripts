#!/usr/bin/env python

from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO


def main():
   f = open("TORO_gp160aa_all.txt","r")
   f2 = open("gp160_full.seq","w")
   
   for line in f:
      a = line.split()
      if len(a[1]) > 520:
         f2.write(">"+a[0]+"\n"+a[1]+"\n")
         
   f2.close()
   f.close()

if __name__ == '__main__':
   main()