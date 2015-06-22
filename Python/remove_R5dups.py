# Removes all copies except one of duplicate sequences
# Bias is to keep dup sequence that is DM or X4
########################################################

from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

filename = raw_input("Enter sequence file containing duplicates: ")
f1 = open(filename+'.noR5dups','w')

nodups = {}
for record in SeqIO.parse(filename,"fasta",generic_dna):
   ac = record.id.split("|")
   if nodups.has_key(record.seq.tostring()):
      if ac[2]=='DM' or ac[2]=='X4' or ac[3]=='DM' or ac[3]=='X4':
         nodups[record.seq.tostring()] = ac
   else: nodups[record.seq.tostring()] = ac

for seq in nodups:
   f1.write(">"+"|".join(nodups[seq])+"\n"+seq+"\n")

f1.close()