

#!/usr/bin/env python
# encoding: utf-8
"""
format4Sequin.py

Created by kaori on 2011-09-21.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO


def main():
   # DNA files
   def1 = "PRRT_"
   def11 = " HIV-1 isolate  from Dominican Republic pol protein (pol) gene, partial cds. [organism=Human Immunodeficiency Virus 1][molecule=DNA][moltype=transcribed RNA][location=virion][country=Dominican Republic][gene=pol]"
   def3 = "[collection-date="
   def4 = "[subtype="
   def5 = "[isolate=HIV_PRRT_"
   note = " encodes full length PR; RT amino acids 1-305"
   p1 = " [gene=pol][protein=pol protein] "
   count = 1
   f = open('PJ01967_nt_ncbi.seq','w')
   f2 = open('PJ01967_aa_ncbi.seq','w')
   f3 = open('submission_clone_lookup.txt','w')
   for x in SeqIO.parse('PJ01967.nt.seq','fasta'):
      vals = x.description.split('|')
      acc = 'PJ01967_'+str(count)   # vals[0]
      c = vals[4].split('/')
      col_date = c[2][:4]
      subtype = 'B'
      if acc == '10-108759': subtype='Complex'
      f.write(">"+def1+acc+def11+def3+str(col_date)+"]"+def5+acc+"]"+def4+subtype+"] "+note+"\n"+x.seq.tostring()+"\n")
      f2.write(">"+def1+acc+p1+note+"\n"+x.seq.translate().tostring()+"\n")
      f3.write(acc+"\t"+vals[0]+"\n")
      count += 1
   
   f.close()
   f2.close()
   f3.close()
   # Protein files
   # 100
#   for y in SeqIO.parse('PJ01967.aa.seq','fasta'):
#      vals = y.description.split('|')
#      acc = vals[0]
#      c = vals[4].split('/')
#      col_date = c[2][:4]
#      subtype='B'
#      if acc = '10-108759': subtype = 'Complex'
#      pr = y.seq.tostring()[:101]
#      rt = y.seq.tostring()[101:]
      


if __name__ == '__main__':
	main()

