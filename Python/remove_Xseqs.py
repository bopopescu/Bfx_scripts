# Removes sequences containing 'X'
########################################################

from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

filename = raw_input("Enter sequence file containing X: ")
f0 = open(filename,'r')
f1 = open(filename+'.clean','w')

for line in f0:
   a = line.split(" ")
   if a[-1].find('X') == -1: f1.write(a[0]+'\t'+a[-1])

f0.close()
f1.close()
