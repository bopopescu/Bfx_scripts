# Removes both copies of duplicate sequences
# Writes all of the removed sequences to a seperate file
########################################################

from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

filename = raw_input("Enter sequence file containing duplicates: ")
f1 = open(filename+'.nodups','w')
f2 = open(filename+'.dups','w')
nodups = {}
for record in SeqIO.parse(filename,"fasta",generic_dna):
   ac = record.id.split("|")
   if nodups.has_key(record.seq.tostring()):
      f2.write(">"+"|".join(ac)+"\n"+record.seq.tostring()+"\n")
      f2.write(">"+"|".join(nodups[record.seq.tostring()])+"\n"+record.seq.tostring()+"\n")
      del nodups[record.seq.tostring()]
   else: nodups[record.seq.tostring()] = ac

for seq in nodups:
   f1.write(">"+"|".join(nodups[seq])+"\n"+seq+"\n")

f1.close()
f2.close()