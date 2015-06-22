# converts 2-column tab file to fasta
########################################################

from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

filename = raw_input("Enter tab file to convert to fasta: ")
f0 = open(filename,'r')
f1 = open(filename+'.seq','w')

for line in f0:
   a = line.split("\t")
   f1.write(">"+a[0]+"\n"+a[1])

f0.close()
f1.close()
